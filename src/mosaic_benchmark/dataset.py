"""Distinct loci, one mutation per condition, paired native-background samples."""
import bisect
import json
import math
import os
from pathlib import Path
import platform
import sys

from .common import (atomic_text, digest, public_config, random_stream, rc,
                     read_manifest, root, sha, write_json)
from .fasta import Fasta
from .models import mutate


def source_code_hashes():
    return {p.name: sha(p) for p in sorted(Path(__file__).parent.glob('*.py'))}


def cases(c):
    """Static metadata: Snakemake can construct the complete DAG before generation."""
    for r in c['references']:
        for n in c['lengths']:
            for identity in c['identities']:
                for model in (['P'] if identity == 100 else c['models']):
                    yield dict(id=f"{r['id']}_L{n}_{model}_{identity}_short", reference_id=r['id'],
                               length=n, model=model, identity=identity, context='short',
                               records=c['loci_per_reference_length'])
            if n in c['background']['lengths']:
                for identity in c['background']['identities']:
                    for model in c['models']:
                        yield dict(id=f"{r['id']}_L{n}_{model}_{identity}_native", reference_id=r['id'],
                                   length=n, model=model, identity=identity, context='native',
                                   records=c['background']['selected_loci_per_reference_length']//c['background']['targets_per_sample'])


def batches(c):
    for case in cases(c):
        size = c['batch'][case['context']]
        for i in range(math.ceil(case['records']/size)):
            yield dict(case, batch=i, id=case['id']+f"_b{i:04d}")


def batch_paths(c, batch):
    base = root(c)/'dataset'/'batches'/batch['id']
    return base.with_suffix('.fa'), base.with_suffix('.jsonl')


def draw(fasta, length, rng, max_attempts, forbidden=(), used_starts=None, occupied=None):
    choices = [(k, max(0, v['length']-length+1)) for k, v in fasta.records.items()]
    ends, total = [], 0
    for _, count in choices:
        total += count
        ends.append(total)
    if total == 0:
        raise ValueError('No reference record can hold the requested length')
    for attempt in range(1, max_attempts+1):
        u = rng.randrange(total)
        j = bisect.bisect_right(ends, u)
        name = choices[j][0]
        start = u-(ends[j-1] if j else 0)
        end = start+length
        if used_starts is not None and (name, start) in used_starts:
            continue
        intervals = list(forbidden.get(name, [])) if isinstance(forbidden, dict) else []
        if occupied is not None:
            entries = occupied.get(name, [])
            at = bisect.bisect_left(entries, (start, end))
            intervals += entries[max(0, at-1):at+1]
        if any(start < b and a < end for a, b in intervals):
            continue
        raw = fasta.fetch(name, start, end)
        if set(raw.upper())-set('ACGT'):
            continue
        return name, start, raw, attempt
    raise RuntimeError(f'Input sampling exhausted {max_attempts} candidates; do not change seed or tool denominator')


def excluded_regions(ref):
    out = {}
    if ref.get('exclude_bed'):
        for line in Path(ref['exclude_bed']).read_text().splitlines():
            if line and not line.startswith('#'):
                name, a, b, *_ = line.split()
                if not 0 <= int(a) < int(b):
                    raise ValueError('Invalid 0-based exclusion BED')
                out.setdefault(name, []).append((int(a), int(b)))
    return out


def sample_sources(c, ref, fasta):
    sources, used, occupied = {}, set(), {}
    n_loci = c['loci_per_reference_length']
    excluded = excluded_regions(ref)
    for length in c['lengths']:
        orientation_rng, _ = random_stream(c['seed'], 'orientation', ref['id'], length)
        signs = ['+']*(n_loci//2) + ['-']*(n_loci-n_loci//2)
        orientation_rng.shuffle(signs)
        rows = []
        for i in range(n_loci):
            rng, h = random_stream(c['seed'], 'source', 'final', ref['id'], length, i+1)
            name, start, raw, attempts = draw(fasta, length, rng, c['sampling_attempt_limit'],
                excluded, used,
                occupied if c.get('nonoverlapping_loci', True) else None)
            used.add((name, start))
            bisect.insort(occupied.setdefault(name, []), (start, start+length))
            sid = f"{ref['id']}_L{length}_s{i+1:04d}"
            rows.append(dict(source_id=sid, reference_id=ref['id'], record=name, start=start,
                end=start+length, strand=signs[i], length=length, sequence=raw.upper(),
                masked_bases=sum(x.islower() for x in raw), source_seed_sha256=h,
                candidates_tried=attempts, dependency_group=sid, split='final'))
        sources[length] = rows
    assert len(used) == len(c['lengths'])*n_loci
    return sources


def truth_variant(c, source, model, identity):
    seq, path = mutate(source['sequence'], identity, 'U' if model == 'P' else model,
                       c['seed'], source['source_id'])
    qid = f"{source['source_id']}_{model}_{identity}"
    if source['strand'] == '-':
        seq = rc(seq)
    truth = {k: v for k, v in source.items() if k != 'sequence'}
    truth.update(path, query_id=qid, fragment_id=qid, model=model, identity=identity,
                 qstart=0, qend=len(seq), query_length=len(seq), context='short',
                 strict_path_orientation='qaln/taln run along forward reference; query reversed if strand=-')
    return seq, truth


def emit_batches(c, folder, case, records):
    size = c['batch'][case['context']]
    for i, begin in enumerate(range(0, len(records), size)):
        name = case['id']+f'_b{i:04d}'
        block = records[begin:begin+size]
        # Bowtie1's FASTA reader stops at the first sequence newline. All tools
        # receive the same unwrapped records, preserving the complete query.
        atomic_text(folder/'batches'/f'{name}.fa', ''.join('>'+qid+'\n'+seq+'\n' for qid, seq, _ in block))
        atomic_text(folder/'batches'/f'{name}.jsonl', ''.join(json.dumps(t, separators=(',', ':'))+'\n'
                    for _, _, truths in block for t in truths))


def input_fingerprints(c):
    result = {}
    for r in [*c['references'], c['background']['donor']]:
        p = Path(r['fasta']).resolve()
        h = sha(p)
        if r.get('sha256') and r['sha256'] != h:
            raise RuntimeError('FASTA checksum differs: '+str(p))
        st = p.stat()
        result[str(p)] = dict(sha256=h, size=st.st_size, mtime_ns=st.st_mtime_ns)
        if r.get('exclude_bed'):
            result[str(Path(r['exclude_bed']).resolve())] = dict(sha256=sha(r['exclude_bed']))
    return result


def prepare(c):
    out = root(c)/'dataset'
    if (out/'manifest.json').exists():
        verify_dataset(c, full=True)
        return
    fingerprints = input_fingerprints(c)
    folder = root(c)/f'.dataset.{os.getpid()}'
    folder.mkdir(parents=True, exist_ok=False)
    case_lookup = {x['id']: x for x in cases(c)}
    source_rows, backgrounds, exact = [], [], {}
    with Fasta(c['background']['donor']['fasta'], root(c)/'cache/donor.fai.json') as donor:
        for ref in c['references']:
            with Fasta(ref['fasta'], root(c)/f"cache/{ref['id']}.fai.json") as fasta:
                if ref.get('record') and set(fasta.records) != {ref['record']}:
                    raise ValueError('Chromosome-only reference contains unexpected records')
                sources = sample_sources(c, ref, fasta)
                for length, rows in sources.items():
                    source_rows.extend(rows)
                    selection, selection_hash = random_stream(c['seed'], 'background_selection', ref['id'], length)
                    chosen = list(range(len(rows)))
                    selection.shuffle(chosen)
                    chosen = chosen[:c['background']['selected_loci_per_reference_length']]
                    k = c['background']['targets_per_sample']
                    for identity in c['identities']:
                        for model in (['P'] if identity == 100 else c['models']):
                            variants = [truth_variant(c, s, model, identity) for s in rows]
                            case_id = f"{ref['id']}_L{length}_{model}_{identity}_short"
                            emit_batches(c, folder, case_lookup[case_id], [(t['query_id'], q, [t]) for q, t in variants])
                            for q, t in variants:
                                h = digest(q)
                                exact.setdefault(h, []).append(t['fragment_id'])
                            if (length not in c['background']['lengths'] or
                                identity not in c['background']['identities'] or model == 'P'):
                                continue
                            composites = []
                            for group_number, begin in enumerate(range(0, len(chosen), k)):
                                group = chosen[begin:begin+k]
                                # Background independent of identity/model: same sampled blocks/layout
                                # for a given four-source group, target variants change as specified.
                                rng, h = random_stream(c['seed'], 'background', ref['id'], length, group_number)
                                count, rem = divmod(c['background']['total_bp'], k+1)
                                chunks, positions, targets, offset = [], [], [], 0
                                qid = f"{ref['id']}_L{length}_{model}_{identity}_g{group_number:04d}"
                                for slot in range(k+1):
                                    width = count+(slot < rem)
                                    if width:
                                        rec, start, raw, tries = draw(donor, width, rng, c['sampling_attempt_limit'])
                                        chunks.append(raw.upper()); offset += width
                                        positions.append(dict(slot=slot, donor_record=rec, start=start,
                                                              end=start+width, candidates_tried=tries))
                                    if slot < k:
                                        q, t = variants[group[slot]]
                                        tt = dict(t, query_id=qid, qstart=offset, qend=offset+len(q),
                                                  context='native', sample_id=qid,
                                                  background_group=f"{ref['id']}_L{length}_g{group_number:04d}")
                                        targets.append(tt); chunks.append(q); offset += len(q)
                                sequence = ''.join(chunks)
                                for t in targets:
                                    t['query_length'] = len(sequence)
                                assert len(sequence) == c['background']['total_bp']+sum(t['qend']-t['qstart'] for t in targets)
                                composites.append((qid, sequence, targets))
                                backgrounds.append(dict(sample_id=qid, seed_sha256=h,
                                    selection_seed_sha256=selection_hash, slots=positions,
                                    fragment_ids=[t['fragment_id'] for t in targets]))
                            emit_batches(c, folder, case_lookup[case_id.replace('_short', '_native')], composites)
    atomic_text(folder/'sources.jsonl', ''.join(json.dumps(s, separators=(',', ':'))+'\n' for s in source_rows))
    atomic_text(folder/'backgrounds.jsonl', ''.join(json.dumps(x, separators=(',', ':'))+'\n' for x in backgrounds))
    write_json(folder/'sequence_aliases.json', {h: ids for h, ids in exact.items() if len(ids) > 1})
    # P100 aliases are semantic model aliases, not three extra query records.
    write_json(folder/'perfect_aliases.json', {'models': ['U', 'C', 'I'], 'canonical_model': 'P',
        'rule': 'One P100 record per distinct locus, all three perfect conditions alias it'})
    files = {str(p.relative_to(folder)): sha(p) for p in sorted(folder.rglob('*')) if p.is_file()}
    write_json(folder/'manifest.json', dict(config_sha256=digest(public_config(c)), config=public_config(c),
        code_sha256=source_code_hashes(), input_files=fingerprints, files_sha256=files,
        source_loci=len(source_rows), planned_batches=len(list(batches(c))),
        python=sys.version, platform=platform.platform(), master_seed=c['seed'],
        truth_role='Generating source/path, not natural orthology or certified equivalent loci',
        scientific_acceptance='NOT_ESTABLISHED'))
    if out.exists():
        # Snakemake creates output parent directories before running a rule.
        # Only remove that empty directory; never discard partial science files.
        if any(out.iterdir()):
            raise RuntimeError('Dataset output already exists without a valid manifest; preserve partial evidence')
        out.rmdir()
    folder.replace(out)


def verify_dataset(c, full=False):
    manifest = read_manifest(root(c)/'dataset/manifest.json')
    if manifest['config_sha256'] != digest(public_config(c)) or manifest['code_sha256'] != source_code_hashes():
        raise RuntimeError('Frozen configuration or source changed; use a new work_dir')
    if full:
        if manifest['input_files'] != input_fingerprints(c):
            raise RuntimeError('Frozen input changed')
        for rel, h in manifest['files_sha256'].items():
            if sha(root(c)/'dataset'/rel) != h:
                raise RuntimeError('Frozen dataset changed: '+rel)
    else:
        for path, record in manifest['input_files'].items():
            if 'size' in record:
                st = Path(path).stat()
                if (st.st_size, st.st_mtime_ns) != (record['size'], record['mtime_ns']):
                    raise RuntimeError('Input metadata changed after initial SHA verification')
    return manifest


def verify_batch_file(c, path, manifest):
    relative = str(path.relative_to(root(c)/'dataset'))
    if sha(path) != manifest['files_sha256'][relative]:
        raise RuntimeError('Frozen batch changed: '+relative)
