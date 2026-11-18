"""Paired-position region90; minimum paired length applies to each native HSP."""
import json
from collections import defaultdict
from .alignment import parse
from .common import read_fasta, read_json, root, sha, write_json
from .dataset import batch_paths, verify_batch_file, verify_dataset
from .fasta import Fasta


def score(truth, alignments, minimum=100, coverage=.9):
    qcovered, tcovered, strict = set(), set(), set()
    generating = {}
    q, t = 0, truth['start']
    for qb, tb in zip(truth['qaln'], truth['taln']):
        if qb != '-' and tb != '-':
            qq = truth['qstart']+q if truth['strand'] == '+' else truth['qend']-1-q
            generating[qq] = t
        q += qb != '-'; t += tb != '-'
    for a in alignments:
        if (a.query != truth['query_id'] or a.reference != truth['record'] or
            a.strand != truth['strand'] or a.paired < minimum):
            continue
        for q0, t0, size, step in a.segments():
            # Intersect the segment with BOTH truth axes before taking position unions.
            lo, hi = max(0, truth['start']-t0), min(size, truth['end']-t0)
            if step == 1:
                lo, hi = max(lo, truth['qstart']-q0), min(hi, truth['qend']-q0)
            else:
                lo, hi = max(lo, q0-truth['qend']+1), min(hi, q0-truth['qstart']+1)
            for i in range(lo, hi):
                qq, tt = q0+step*i, t0+i
                qcovered.add(qq); tcovered.add(tt)
                if generating.get(qq) == tt:
                    strict.add(qq)
    qcov = len(qcovered)/(truth['qend']-truth['qstart'])
    tcov = len(tcovered)/(truth['end']-truth['start'])
    return dict(query_coverage=qcov, source_coverage=tcov,
                region90=int(qcov >= coverage and tcov >= coverage),
                strict_generating_pair_recall=len(strict)/len(generating) if generating else None)


def evaluate_job(c, tool, batch):
    dataset = verify_dataset(c)
    batch_id = batch['id']
    folder = root(c)/'runs'/tool/batch_id
    record = read_json(folder/'run.json')
    qpath, tpath = batch_paths(c, batch)
    verify_batch_file(c, tpath, dataset)
    if record['query_sha256'] != sha(qpath):
        raise RuntimeError('Query changed since mapping')
    queries = read_fasta(qpath)
    ref = next(r for r in c['references'] if r['id'] == batch['reference_id'])
    with Fasta(ref['fasta'], root(c)/f"cache/{ref['id']}.fai.json") as fasta:
        alignments = defaultdict(list)
        if record['status'] == 'SUCCESS':
            if sha(record['raw_path']) != record['raw_sha256']:
                raise RuntimeError('Native output changed')
            for a in parse(record['raw_path'], record['format'], {k: len(v) for k, v in queries.items()}, fasta.lengths()):
                alignments[a.query].append(a)
        rows = []
        for line in tpath.read_text().splitlines():
            t = json.loads(line)
            result = score(t, alignments[t['query_id']], c['metrics']['minimum_paired_bp'], c['metrics']['coverage'])
            if record['status'] in ['NOT_SUPPORTED', 'INDEX_FAILED', 'DEPENDENCY_FAILED']:
                result = {k: None for k in result}
            rows.append({k: t[k] for k in ['fragment_id', 'query_id', 'source_id', 'reference_id', 'length',
                         'model', 'identity', 'context', 'dependency_group']} | result | dict(
                tool=tool, batch_id=batch_id, native_status=record['status'],
                background_group=t.get('background_group'), control=t['identity'] == 100))
    write_json(folder/'evaluation.json', dict(rows=rows, native_status=record['status'],
               raw_sha256=record.get('raw_sha256'), metric=c['metrics']))
