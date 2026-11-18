"""Bounded native execution, immutable attempts, indexes, and six tool policies."""
import fcntl
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

from .common import read_fasta, read_json, read_manifest, root, sha, write_json
from .dataset import batch_paths, batches, verify_batch_file, verify_dataset
from .fasta import Fasta

POLICIES = {
    'minimap2': ['-x', 'asm20', '-k11', '-w5', '-c', '--cs=long', '-A1', '-B1',
                 '--score-N=0', '-O4,8', '-E2,1', '-s20', '-m20', '-n2', '-p0.8'],
    'last': ['-j3', '-fMAF', '-m100'],
    'lastz': ['--format=maf'],
    'blastn': ['-task', 'blastn', '-word_size', '7', '-reward', '1', '-penalty', '-1',
               '-gapopen', '2', '-gapextend', '1', '-dust', 'yes', '-soft_masking', 'true',
               '-strand', 'both', '-evalue', '1e-6', '-max_target_seqs', '500', '-parse_deflines'],
    # Explicit new task policy: not the historical long-genome -s1000 policy.
    'fastga': ['-k', '-i0.7', '-f10', '-c85', '-s100', '-l100'],
    'bowtie': ['-f', '-S', '-v', '3', '-a', '--best', '--strata'],
}
EXECUTABLES = {
    'minimap2': ['minimap2'], 'last': ['lastdb', 'lastal'], 'lastz': ['lastz'],
    'blastn': ['makeblastdb', 'blastn'], 'bowtie': ['bowtie-build', 'bowtie'],
    'fastga': ['FAtoGDB', 'GIXmake', 'GDBtoFA', 'FastGA', 'ALNtoPAF'],
}
FORMATS = dict(minimap2='paf', last='maf', lastz='maf', blastn='blast', fastga='paf', bowtie='sam')


def executable(c, name):
    candidate = c.get('executables', {}).get(name, name)
    found = shutil.which(candidate)
    if not found:
        raise FileNotFoundError(f'Missing executable {name}: configure executables or PATH')
    return str(Path(found).resolve())


def hardware():
    cpu = Path('/proc/cpuinfo').read_text() if Path('/proc/cpuinfo').exists() else ''
    model = next((s.split(':', 1)[1].strip() for s in cpu.splitlines() if s.startswith('model name')), None)
    return dict(cpu_model=model, affinity=sorted(os.sched_getaffinity(0)), hostname=os.uname().nodename,
                captured_at=time.time(), load_average=os.getloadavg())


def doctor(c):
    if (root(c)/'toolchain.json').exists():
        for tool in c['enabled_tools']:
            verify_tool(c, tool)
        return read_json(root(c)/'toolchain.json')['tools']
    records = {}
    for tool in c['enabled_tools']:
        entries = {}
        for name in EXECUTABLES[tool]:
            path = executable(c, name)
            result = subprocess.run([path, '--version'], capture_output=True, text=True, timeout=15)
            version = (result.stdout+result.stderr)[:12000]
            if name == 'bowtie' and ('bowtie2' in version.lower() or 'bowtie2' in path.lower()):
                raise ValueError('Bowtie2 is not original Bowtie')
            entries[name] = dict(path=path, sha256=sha(path), version_output=version,
                                 version_returncode=result.returncode)
            # Bowtie's launcher executes these adjacent binaries; pin those too.
            if name == 'bowtie':
                for p in sorted(Path(path).parent.glob('bowtie*-?')):
                    if p.is_file():
                        entries[p.name] = dict(path=str(p), sha256=sha(p))
        records[tool] = entries
    write_json(root(c)/'toolchain.json', dict(tools=records, policies=POLICIES, hardware=hardware(),
        fastga_build_label=c.get('fastga_build_label', 'upstream'),
        qualification='Executable inventory only; perfect controls and parser tests provide separate evidence'))
    return records


def verify_tool(c, tool):
    inventory = read_json(root(c)/'toolchain.json')
    for entry in inventory['tools'][tool].values():
        if sha(entry['path']) != entry['sha256']:
            raise RuntimeError('Pinned executable/helper changed')
    if inventory['policies'] != POLICIES:
        raise RuntimeError('Frozen native policies changed')


def time_values(path):
    result = {'gnu_wall_seconds': None, 'cpu_seconds': None, 'max_rss_kib': None}
    if not path.exists():
        return result
    vals = {}
    for line in path.read_text().splitlines():
        if ': ' in line:
            key, value = line.strip().rsplit(': ', 1)
            vals[key] = value
    try:
        result['cpu_seconds'] = float(vals['User time (seconds)'])+float(vals['System time (seconds)'])
        result['max_rss_kib'] = int(vals['Maximum resident set size (kbytes)'])
        wall = vals['Elapsed (wall clock) time (h:mm:ss or m:ss)']
        total = 0.
        for part in wall.split(':'):
            total = total*60+float(part)
        result['gnu_wall_seconds'] = total
    except (KeyError, ValueError):
        pass
    return result


def measured(c, argv, folder, timeout, index=False):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    previous = sorted(folder.glob('attempt_*'))
    attempt = folder/f'attempt_{len(previous)+1:04d}'
    attempt.mkdir()
    write_json(attempt/'command.json', dict(argv=argv, timeout_seconds=timeout,
        hardware=hardware(), memory_limit_gib=c['resources']['memory_gib'],
        purpose='index_or_preparation' if index else 'mapping_or_conversion'))
    # Shared lock protects the total budget even outside Snakemake's native_slot.
    lock = root(c)/'native.lock'
    lock.parent.mkdir(parents=True, exist_ok=True)
    with lock.open('a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        bp = root(c)/'budget.json'
        budget = read_json(bp) if bp.exists() else dict(spent_seconds=0., completed_attempts=0)
        if budget.get('reserved_attempt'):
            raise RuntimeError('Unsettled prior attempt: reconcile evidence before resuming')
        if budget['spent_seconds']+timeout > c['resources']['total_native_budget_seconds']:
            raise RuntimeError('Total native budget cannot admit the next timeout reservation')
        budget['reserved_attempt'] = str(attempt)
        budget['reserved_seconds'] = timeout
        write_json(bp, budget)
        raw, err, tv = attempt/'stdout.partial', attempt/'stderr.txt', attempt/'time.txt'
        cpus = c['resources'].get('cpu_affinity')
        if cpus is None:
            cpus = sorted(os.sched_getaffinity(0))[:c['resources']['threads']]
        memory = int(c['resources']['memory_gib']*1024**3)
        size = c['resources']['index_output_limit_mib' if index else 'output_limit_mib']*1024**2
        worker = [sys.executable, '-m', 'mosaic_benchmark.worker', '--cpus', ','.join(map(str, cpus)),
                  '--memory-bytes', str(memory), '--output-bytes', str(size)]
        # Apply the named1MiB stack contract to GIXmake and all FastGA helpers,
        # too: GIXmake launches many joinable workers under the AS budget.
        is_fastga = Path(argv[0]).name in EXECUTABLES['fastga']
        if is_fastga:
            worker += ['--stack-bytes', str(1024**2)]
        command_path = attempt/'command.json'
        write_json(command_path, read_json(command_path) | dict(resource_wrapper_argv=worker,
            effective_cpu_affinity=cpus, fastga_stack_bytes=1024**2 if is_fastga else None))
        started = time.monotonic()
        status, returncode = 'FAILED', None
        env = dict(os.environ, LC_ALL='C')
        helper_dirs = {str(Path(executable(c, n)).parent) for n in EXECUTABLES.get('fastga', [])} if 'fastga' in c['enabled_tools'] else set()
        env['PATH'] = os.pathsep.join(sorted(helper_dirs)+[env.get('PATH', '')])
        with raw.open('wb') as stdout, err.open('wb') as stderr:
            process = subprocess.Popen([*worker, '--', '/usr/bin/time', '-v', '-o', str(tv), *argv], stdout=stdout,
                stderr=stderr, env=env, start_new_session=True)
            try:
                returncode = process.wait(timeout=timeout)
                status = 'SUCCESS' if returncode == 0 else 'FAILED'
            except subprocess.TimeoutExpired:
                status = 'TIMEOUT'
                os.killpg(process.pid, signal.SIGTERM)
                try:
                    process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    process.wait()
                returncode = process.returncode
        elapsed = time.monotonic()-started
        if status == 'SUCCESS':
            output = attempt/'stdout.txt'
            raw.replace(output); raw = output
        receipt = dict(status=status, returncode=returncode, wrapper_wall_seconds=elapsed,
            raw_path=str(raw), raw_sha256=sha(raw), command_path=str(attempt/'command.json'),
            stderr_path=str(err), time_path=str(tv), **time_values(tv))
        write_json(attempt/'receipt.json', receipt)
        budget['spent_seconds'] += elapsed
        budget['completed_attempts'] += 1
        budget['last_receipt'] = str(attempt/'receipt.json')
        budget.pop('reserved_attempt'); budget.pop('reserved_seconds')
        write_json(bp, budget)
        return receipt


def fastga_database(c, fasta, prefix, folder):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    result = []
    for name, args in [
        ('FAtoGDB', [str(fasta), str(prefix)]),
        ('GIXmake', [f"-T{c['resources']['threads']}", '-k40', '-P'+str(folder), str(prefix)+'.1gdb']),
        ('GDBtoFA', ['-U', str(prefix)+'.1gdb'])]:
        stat = measured(c, [executable(c, name), *args], Path(folder)/name,
                        c['resources']['index_timeout_seconds'], index=True)
        result.append(stat)
        if stat['status'] != 'SUCCESS':
            return result
    # Validate sequence names, lengths, and all bases, streaming for large genomes.
    with Fasta(fasta) as expected, Fasta(result[-1]['raw_path']) as observed:
        if expected.lengths() != observed.lengths():
            raise RuntimeError('FastGA GDB readback record/length mismatch')
        for name, n in expected.lengths().items():
            for start in range(0, n, 65536):
                if expected.fetch(name, start, min(n, start+65536)).upper() != observed.fetch(name, start, min(n, start+65536)).upper():
                    raise RuntimeError('FastGA GDB full sequence readback mismatch')
    return result


def build_index(c, tool, reference_id):
    dataset = verify_dataset(c); verify_tool(c, tool)
    ref = next(r for r in c['references'] if r['id'] == reference_id)
    folder = root(c)/'indexes'/reference_id/tool
    folder.mkdir(parents=True, exist_ok=True)
    prefix = folder/'db'
    threads = str(c['resources']['threads'])
    argv = {
        'minimap2': ['minimap2', '-x', 'asm20', '-k11', '-w5', '-t', threads, '-d', str(prefix)+'.mmi', ref['fasta']],
        'last': ['lastdb', '-P'+threads, '-uYASS', '-c', str(prefix), ref['fasta']],
        'blastn': ['makeblastdb', '-dbtype', 'nucl', '-blastdb_version', '4', '-parse_seqids', '-in', ref['fasta'], '-out', str(prefix)],
        'bowtie': ['bowtie-build', '-f', ref['fasta'], str(prefix)],
    }
    if tool == 'lastz':
        stats = []
    elif tool == 'fastga':
        stats = fastga_database(c, ref['fasta'], prefix, folder)
    else:
        command = argv[tool]; command[0] = executable(c, command[0])
        stats = [measured(c, command, folder/'build', c['resources']['index_timeout_seconds'], index=True)]
    files = sorted({p for pattern in ['db*', '.db*'] for p in folder.glob(pattern) if p.is_file()})
    write_json(folder/'manifest.json', dict(tool=tool, reference_id=reference_id, prefix=str(prefix),
        reference_sha256=dataset['input_files'][str(Path(ref['fasta']).resolve())]['sha256'],
        status='SUCCESS' if all(s['status'] == 'SUCCESS' for s in stats) else 'INDEX_FAILED', measurements=stats,
        files_sha256={str(p): sha(p) for p in files},
        files_stat={str(p): [p.stat().st_size, p.stat().st_mtime_ns] for p in files}))


def verify_index(index, full=False):
    for name, expected in index['files_sha256'].items():
        path = Path(name)
        stat = path.stat()
        changed = [stat.st_size, stat.st_mtime_ns] != index['files_stat'][name]
        if changed or (full and sha(path) != expected):
            raise RuntimeError('Frozen reference index changed: '+name)


def mapping_argv(c, tool, reference, query, prefix):
    threads = str(c['resources']['threads'])
    if tool == 'minimap2':
        return [executable(c, 'minimap2'), *POLICIES[tool], '-t', threads, str(prefix)+'.mmi', str(query)]
    if tool == 'last':
        return [executable(c, 'lastal'), *POLICIES[tool], '-P'+threads, str(prefix), str(query)]
    if tool == 'lastz':
        return [executable(c, 'lastz'), reference+'[multiple]', str(query), *POLICIES[tool]]
    if tool == 'blastn':
        return [executable(c, 'blastn'), *POLICIES[tool], '-num_threads', threads, '-outfmt',
                '6 qseqid sseqid qstart qend sstart send qseq sseq', '-query', str(query), '-db', str(prefix)]
    if tool == 'bowtie':
        return [executable(c, 'bowtie'), *POLICIES[tool], '-p', threads, str(prefix), str(query)]
    raise ValueError('FastGA has a separate indexed-query pipeline')


def map_job(c, tool, b):
    dataset = verify_dataset(c); verify_tool(c, tool)
    batch_id = b['id']
    query, _ = batch_paths(c, b)
    verify_batch_file(c, query, dataset)
    index = read_manifest(root(c)/'indexes'/b['reference_id']/tool/'manifest.json')
    folder = root(c)/'runs'/tool/batch_id
    folder.mkdir(parents=True, exist_ok=True)
    record = dict(tool=tool, batch=b, query_sha256=sha(query), format=FORMATS[tool],
                  index_manifest_sha256=sha(root(c)/'indexes'/b['reference_id']/tool/'manifest.json'))
    if tool == 'bowtie' and any(len(x) > 1024 for x in read_fasta(query).values()):
        record.update(status='NOT_SUPPORTED', reason='Original Bowtie maximum read length1024; full composite query retained')
    elif index['status'] != 'SUCCESS':
        record.update(status='INDEX_FAILED')
    else:
        verify_index(index)
        ref = next(r for r in c['references'] if r['id'] == b['reference_id'])
        if tool != 'fastga':
            stat = measured(c, mapping_argv(c, tool, ref['fasta'], query, index['prefix']),
                            folder/'mapping', c['resources']['mapping_timeout_seconds'])
        else:
            qp = folder/'query_db'/'db'
            prep = fastga_database(c, query, qp, folder/'query_preparation')
            record['query_preparation_measurements'] = prep
            if not all(s['status'] == 'SUCCESS' for s in prep):
                record.update(status='DEPENDENCY_FAILED')
                write_json(folder/'run.json', record)
                return
            aln = folder/'alignment.1aln'
            argv = [executable(c, 'FastGA'), '-T'+str(c['resources']['threads']), *POLICIES[tool],
                    '-P'+str(folder), '-1:'+str(aln), str(qp)+'.gix', index['prefix']+'.gix']
            stat = measured(c, argv, folder/'mapping', c['resources']['mapping_timeout_seconds'])
            if stat['status'] == 'SUCCESS':
                record['mapping_measurement'] = stat
                conversion = measured(c, [executable(c, 'ALNtoPAF'), '-x', '-T'+str(c['resources']['threads']), str(aln)],
                                      folder/'conversion', c['resources']['conversion_timeout_seconds'])
                record['conversion_measurement'] = conversion
                if conversion['status'] == 'SUCCESS':
                    stat = dict(stat, raw_path=conversion['raw_path'], raw_sha256=conversion['raw_sha256'])
                else:
                    stat = dict(conversion, status='DEPENDENCY_FAILED')
        record.update(stat)
    write_json(folder/'run.json', record)


def perfect_gate(c, tool, reference_id):
    checks = []
    for b in batches(c):
        if b['reference_id'] == reference_id and b['model'] == 'P':
            rows = read_json(root(c)/'runs'/tool/b['id']/'evaluation.json')['rows']
            checks.extend(r['native_status'] == 'SUCCESS' and r['region90'] == 1 for r in rows)
    passed = bool(checks) and all(checks)
    write_json(root(c)/'gates'/f'{tool}_{reference_id}.json',
               dict(status='PASS' if passed else 'FAIL', perfect_queries=len(checks)))
    if not passed:
        raise RuntimeError(f'Perfect source controls failed for {tool}/{reference_id}; raw attempts preserved')
