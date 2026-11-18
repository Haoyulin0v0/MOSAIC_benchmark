"""Plan, generate, run/resume a single DAG, and locate independent-service logs."""
import argparse
import fcntl
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from .common import digest, load_config, read_json, root, write_fasta, write_json
from .dataset import batches, cases, prepare, verify_dataset
from .materials import fetch
from .native import doctor, hardware, verify_index
from .report import plot


def plan(c):
    lengths, identities, refs, n = c['lengths'], c['identities'], c['references'], c['loci_per_reference_length']
    source_count = len(refs)*len(lengths)*n
    short_mutated = source_count*sum(i != 100 for i in identities)*len(c['models'])
    short_perfect = source_count if 100 in identities else 0
    long_count = sum(x['records'] for x in cases(c) if x['context'] == 'native')
    all_batches = list(batches(c))
    per_tool = {}
    for tool in c['enabled_tools']:
        applicable = [b for b in all_batches if not (tool == 'bowtie' and b['context'] == 'native')]
        mapping_calls = len(applicable)
        # FastGA: FAtoGDB + GIXmake + readback per query batch; separate PAF conversion.
        query_preparation_calls = 3*mapping_calls if tool == 'fastga' else 0
        conversions = mapping_calls if tool == 'fastga' else 0
        index_calls = len(refs)*(3 if tool == 'fastga' else 0 if tool == 'lastz' else 1)
        per_tool[tool] = dict(mapping_calls_upper=mapping_calls, index_calls_upper=index_calls,
            query_preparation_calls_upper=query_preparation_calls, conversion_calls_upper=conversions,
            not_supported_long_batches=sum(b['context'] == 'native' for b in all_batches) if tool == 'bowtie' else 0)
    res = c['resources']
    bound = sum(x['mapping_calls_upper']*res['mapping_timeout_seconds'] +
                (x['index_calls_upper']+x['query_preparation_calls_upper'])*res['index_timeout_seconds'] +
                x['conversion_calls_upper']*res['conversion_timeout_seconds'] for x in per_tool.values())
    return dict(source_loci=source_count, nonperfect_short_queries=short_mutated,
        canonical_perfect_queries=short_perfect, long_native_samples=long_count,
        long_target_evaluations=long_count*c['background']['targets_per_sample'],
        cases=len(list(cases(c))), batches_per_tool=len(all_batches), tools=per_tool,
        timeout_reservation_bound_seconds=bound, authorized_native_budget_seconds=res['total_native_budget_seconds'],
        timing_note='Timeout reservation bound, not predicted duration; actual deployment must be resource-qualified',
        scope='Design only; no download, generation, indexing or native call')


def snakefile():
    source = Path(__file__).resolve().parents[2]/'workflow/Snakefile'
    installed = Path(sys.prefix)/'share/mosaic-benchmark/workflow/Snakefile'
    path = source if source.exists() else installed
    if not path.exists():
        raise FileNotFoundError('Packaged workflow/Snakefile missing')
    return path


def foreground(c, dry_run=False):
    folder = root(c)
    folder.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, '-m', 'snakemake', '--snakefile', str(snakefile()), '--configfile',
        c['_path'], '--cores', str(c['resources']['threads']), '--resources', 'native_slot=1',
        '--rerun-incomplete', '--force-use-threads']
    if dry_run:
        return subprocess.call(command+['--dry-run'])
    lock = (folder/'supervisor.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        raise RuntimeError('A supervisor already owns this work_dir')
    # Check hashes before scheduler skip logic can reuse an old output.
    if (folder/'dataset/manifest.json').exists():
        verify_dataset(c, full=True)
    if (folder/'toolchain.json').exists():
        doctor(c)
    for path in (folder/'indexes').glob('*/*/manifest.json'):
        verify_index(read_json(path), full=True)
    write_json(folder/'execution_service.json', dict(status='RUNNING', config=c['_path'],
        command=command, hardware=hardware(), started_at=datetime.now(timezone.utc).isoformat()))
    start = time.monotonic()
    with (folder/'execution.log').open('a') as log:
        result = subprocess.call(command, stdout=log, stderr=subprocess.STDOUT)
    receipt = dict(status='COMPLETED' if result == 0 else 'FAILED', returncode=result,
                   elapsed_seconds=time.monotonic()-start, log=str(folder/'execution.log'),
                   finished_at=datetime.now(timezone.utc).isoformat())
    write_json(folder/'execution_exit.json', receipt)
    write_json(folder/'execution_service.json', receipt)
    with (folder/'execution_history.jsonl').open('a') as history:
        history.write(json.dumps(receipt)+'\n')
    lock.close()
    print(json.dumps(receipt))
    return result


def background(c):
    folder = root(c)
    folder.mkdir(parents=True, exist_ok=True)
    info = folder/'execution_submission.json'
    if info.exists():
        old = read_json(info)
        unit = old.get('unit')
        if unit:
            status = subprocess.run(['systemctl', '--user', 'show', unit, '--property=ActiveState', '--value'],
                                    capture_output=True, text=True)
            if status.returncode == 0 and status.stdout.strip() in ['active', 'activating']:
                print(json.dumps(dict(status='ALREADY_RUNNING', unit=unit, log=str(folder/'execution.log'))))
                return 0
    unit = 'mosaic-'+digest(str(folder))[:12]+'-'+str(int(time.time()))
    argv = ['systemd-run', '--user', '--unit', unit,
        '--property=StandardOutput=append:'+str(folder/'execution.log'),
        '--property=StandardError=append:'+str(folder/'execution.log'),
        '--working-directory='+str(Path.cwd())]
    for variable in ['PATH', 'PYTHONPATH']:
        if os.environ.get(variable):
            value = os.environ[variable]
            if variable == 'PYTHONPATH':
                value = os.pathsep.join(str(Path(x).resolve()) for x in value.split(os.pathsep) if x)
            argv.append('--setenv='+variable+'='+value)
    argv += [sys.executable, '-m', 'mosaic_benchmark',
        'run', '--config', c['_path'], '--foreground']
    subprocess.run(argv, check=True)
    state = subprocess.run(['systemctl', '--user', 'show', unit, '--property=ActiveState,SubState,MainPID'],
                           capture_output=True, text=True, check=True)
    receipt = dict(unit=unit, status='SERVICE_SUBMITTED', systemd_state=state.stdout,
        log=str(folder/'execution.log'), service=str(info), exit=str(folder/'execution_exit.json'))
    # The service writes RUNNING/exit state; retain submission separately to avoid races.
    write_json(folder/'execution_submission.json', receipt)
    print(json.dumps(receipt, indent=2))
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ['plan', 'prepare', 'doctor', 'verify', 'run', 'resume', 'logs', 'plot']:
        p = sub.add_parser(name)
        p.add_argument('--config', required=True)
        if name in ['run', 'resume']:
            p.add_argument('--foreground', action='store_true')
            p.add_argument('--dry-run', action='store_true')
    p = sub.add_parser('fetch-references'); p.add_argument('--manifest', default='config/references.json'); p.add_argument('--directory', default='data')
    args = parser.parse_args(argv)
    if args.command == 'fetch-references':
        fetch(args.manifest, args.directory)
        return 0
    c = load_config(args.config)
    if args.command == 'plan':
        print(json.dumps(plan(c), indent=2))
    elif args.command == 'prepare':
        prepare(c)
    elif args.command == 'doctor':
        doctor(c)
    elif args.command == 'verify':
        print(json.dumps(dict(status='PASS_FROZEN_DATASET', source_loci=verify_dataset(c, full=True)['source_loci'])))
    elif args.command in ['run', 'resume']:
        return foreground(c, args.dry_run) if args.foreground or args.dry_run else background(c)
    elif args.command == 'logs':
        print(json.dumps({name: str(root(c)/name) for name in ['execution.log', 'execution_service.json',
            'execution_submission.json', 'execution_exit.json', 'budget.json', 'summary/summary.json']}, indent=2))
    elif args.command == 'plot':
        plot(c)
    return 0
