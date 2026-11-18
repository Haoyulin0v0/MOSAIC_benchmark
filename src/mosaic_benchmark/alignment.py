"""PAF-CIGAR, pairwise MAF, SAM, and BLAST aligned-string normalization."""
from dataclasses import dataclass
import re
from .common import rc

CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')


def cigar_ops(text):
    ops = [(int(n), op) for n, op in CIGAR.findall(text)]
    if not ops or ''.join(f'{n}{op}' for n, op in ops) != text or any(n < 1 for n, _ in ops):
        raise ValueError('Malformed CIGAR: '+text)
    return ops


@dataclass
class Alignment:
    query: str
    reference: str
    qstart: int
    qend: int
    tstart: int
    tend: int
    strand: str
    ops: list

    def segments(self):
        q, step, t = (self.qstart, 1, self.tstart) if self.strand == '+' else (self.qend-1, -1, self.tstart)
        for n, op in self.ops:
            if op in 'M=X':
                yield q, t, n, step
                q += step*n; t += n
            elif op == 'I':
                q += step*n
            elif op in 'DN':
                t += n
            else:
                raise ValueError('Unexpected operation inside aligned span: '+op)
        if t != self.tend or q != (self.qend if step == 1 else self.qstart-1):
            raise ValueError('CIGAR endpoint consumption disagrees')

    @property
    def paired(self):
        return sum(n for n, op in self.ops if op in 'M=X')


def string_ops(qa, ta):
    if len(qa) != len(ta):
        raise ValueError('Aligned string lengths differ')
    ops = []
    for q, t in zip(qa.upper(), ta.upper()):
        if q == t == '-':
            raise ValueError('Both-gap alignment column')
        op = 'I' if t == '-' else 'D' if q == '-' else '=' if q == t else 'X'
        if ops and ops[-1][1] == op:
            ops[-1] = (ops[-1][0]+1, op)
        else:
            ops.append((1, op))
    return ops


def parse(path, fmt, query_lengths, reference_lengths):
    def checked(a):
        if a.query not in query_lengths or a.reference not in reference_lengths:
            raise ValueError('Native names do not match public FASTA')
        if not (0 <= a.qstart < a.qend <= query_lengths[a.query] and
                0 <= a.tstart < a.tend <= reference_lengths[a.reference] and a.strand in ['+', '-']):
            raise ValueError('Native alignment bounds/strand invalid')
        # Force endpoint validation even for alignments below the paired threshold.
        for _ in a.segments():
            pass
        return a
    if fmt == 'maf':
        block = []
        def convert(rows):
            if len(rows) != 2:
                raise ValueError('Expected exactly two MAF sequence lines')
            t, q = rows
            tn, ts, tz, tstrand, tl, ta = t
            qn, qs, qz, qstrand, ql, qa = q
            ts, tz, tl, qs, qz, ql = map(int, [ts, tz, tl, qs, qz, ql])
            if tl != reference_lengths.get(tn) or ql != query_lengths.get(qn):
                raise ValueError('MAF record sizes differ from inputs')
            if len(ta.replace('-', '')) != tz or len(qa.replace('-', '')) != qz:
                raise ValueError('MAF aligned size mismatch')
            ts = ts if tstrand == '+' else tl-ts-tz
            qs = qs if qstrand == '+' else ql-qs-qz
            if tstrand == '-':
                qa, ta = rc(qa), rc(ta)
            return checked(Alignment(qn, tn, qs, qs+qz, ts, ts+tz,
                                     '+' if qstrand == tstrand else '-', string_ops(qa, ta)))
        with open(path) as f:
            for line in f:
                if line.startswith('a ') or not line.strip():
                    if block:
                        yield convert(block); block = []
                if line.startswith('s '):
                    block.append(line.split()[1:7])
            if block:
                yield convert(block)
        return
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#', '@')):
                continue
            v = line.rstrip().split('\t')
            if fmt == 'paf':
                if len(v) < 12:
                    raise ValueError('PAF columns missing')
                tags = {x[:2]: x[5:] for x in v[12:] if len(x) >= 5}
                if 'cg' not in tags:
                    raise ValueError('PAF lacks base-level CIGAR; bounding boxes are insufficient')
                # FastGA may retain the descriptive part of an original FASTA header.
                if v[0] not in query_lengths and v[0].split()[0] in query_lengths:
                    v[0] = v[0].split()[0]
                if v[5] not in reference_lengths and v[5].split()[0] in reference_lengths:
                    v[5] = v[5].split()[0]
                if int(v[1]) != query_lengths.get(v[0]) or int(v[6]) != reference_lengths.get(v[5]):
                    raise ValueError('PAF record lengths differ')
                yield checked(Alignment(v[0], v[5], int(v[2]), int(v[3]), int(v[7]), int(v[8]), v[4], cigar_ops(tags['cg'])))
            elif fmt == 'sam':
                if len(v) < 11:
                    raise ValueError('SAM columns missing')
                flag = int(v[1])
                if flag & 4:
                    continue
                ops = cigar_ops(v[5])
                leading = sum(n for n, op in ops[:2] if op in 'SH') if ops[0][1] in 'SH' else 0
                trailing = sum(n for n, op in ops[-2:] if op in 'SH') if ops[-1][1] in 'SH' else 0
                ops = [(n, op) for n, op in ops if op not in 'SHP']
                qlen = query_lengths[v[0]]
                reverse = bool(flag & 16)
                qs, qe = (trailing, qlen-leading) if reverse else (leading, qlen-trailing)
                ts = int(v[3])-1
                te = ts+sum(n for n, op in ops if op in 'MDN=X')
                yield checked(Alignment(v[0], v[2], qs, qe, ts, te, '-' if reverse else '+', ops))
            elif fmt == 'blast':
                if len(v) != 8:
                    raise ValueError('Expected BLAST qseqid sseqid qstart qend sstart send qseq sseq')
                qn, tn = v[:2]
                if tn.startswith('lcl|') and tn not in reference_lengths:
                    tn = tn[4:]
                q1, q2, t1, t2 = map(int, v[2:6])
                qa, ta = v[6:]
                if t1 > t2:
                    qa, ta = rc(qa), rc(ta)
                strand = '+' if (q1 <= q2) == (t1 <= t2) else '-'
                yield checked(Alignment(qn, tn, min(q1, q2)-1, max(q1, q2), min(t1, t2)-1,
                                        max(t1, t2), strand, string_ops(qa, ta)))
            else:
                raise ValueError('Unknown native format')
