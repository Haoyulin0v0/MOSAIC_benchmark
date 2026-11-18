"""Streaming per-target results, strata, paired changes, and separate costs."""
from collections import defaultdict
import csv
from pathlib import Path
import sqlite3

from .common import read_json, root, write_json
from .dataset import batches, verify_dataset


def tsv(path, rows, columns):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: '' if row.get(k) is None else row.get(k) for k in columns})


def summarize(c):
    manifest = verify_dataset(c, full=True)
    out = root(c)/'summary'
    out.mkdir(exist_ok=True)
    db = sqlite3.connect(':memory:')
    db.execute('create table scores(tool text, fragment text, source text, ref text, length integer, model text, identity integer, context text, recovered integer, qcov real, tcov real, status text, background_group text)')
    grouped, costs, index_costs = defaultdict(list), [], []
    columns = ['tool', 'fragment_id', 'query_id', 'source_id', 'reference_id', 'length', 'model', 'identity',
               'context', 'native_status', 'region90', 'query_coverage', 'source_coverage',
               'strict_generating_pair_recall', 'dependency_group', 'background_group', 'control', 'batch_id']
    def rows():
        for tool in c['enabled_tools']:
            for ref in c['references']:
                index = read_json(root(c)/'indexes'/ref['id']/tool/'manifest.json')
                for i, stat in enumerate(index['measurements']):
                    index_costs.append(dict(tool=tool, reference_id=ref['id'], stage=i, **stat))
            for b in batches(c):
                folder = root(c)/'runs'/tool/b['id']
                run = read_json(folder/'run.json')
                costs.append(dict(tool=tool, batch_id=b['id'], reference_id=b['reference_id'],
                    context=b['context'], length=b['length'], status=run['status'],
                    mapping_wall_seconds=run.get('gnu_wall_seconds'),
                    mapping_cpu_seconds=run.get('cpu_seconds'), mapping_max_rss_kib=run.get('max_rss_kib'),
                    query_preparation_wall_seconds=sum(s['wrapper_wall_seconds'] for s in run.get('query_preparation_measurements', [])),
                    conversion_wall_seconds=run.get('conversion_measurement', {}).get('wrapper_wall_seconds')))
                evaluation = read_json(folder/'evaluation.json')
                for r in evaluation['rows']:
                    if r['native_status'] in ['DEPENDENCY_FAILED', 'INDEX_FAILED']:
                        raise RuntimeError('Technical dependency failed; retain partial results, do not publish complete benchmark')
                    key = (tool, r['reference_id'], r['length'], r['model'], r['identity'], r['context'])
                    grouped[key].append(r['region90'])
                    db.execute('insert into scores values (?,?,?,?,?,?,?,?,?,?,?,?,?)',
                        (tool, r['fragment_id'], r['source_id'], r['reference_id'], r['length'], r['model'],
                         r['identity'], r['context'], r['region90'], r['query_coverage'], r['source_coverage'],
                         r['native_status'], r['background_group']))
                    yield r
    tsv(out/'per_query.tsv', rows(), columns)
    cells = []
    for key, values in sorted(grouped.items()):
        eligible = [v for v in values if v is not None]
        cells.append(dict(zip(['tool', 'reference_id', 'length', 'model', 'identity', 'context'], key)) |
                     dict(planned=len(values), evaluated=len(eligible), not_supported=len(values)-len(eligible),
                          recovered=sum(eligible), recovery_fraction=sum(eligible)/len(eligible) if eligible else None))
    tsv(out/'grouped.tsv', cells, ['tool', 'reference_id', 'length', 'model', 'identity', 'context',
                                 'planned', 'evaluated', 'not_supported', 'recovered', 'recovery_fraction'])
    db.execute('create index matching on scores(tool,fragment,context)')
    pairs = []
    query = '''select n.tool,n.ref,n.length,n.model,n.identity,count(*),sum(s.recovered),sum(n.recovered),
               sum(case when s.recovered=0 and n.recovered=1 then 1 else 0 end),
               sum(case when s.recovered=1 and n.recovered=0 then 1 else 0 end)
               from scores n join scores s on n.tool=s.tool and n.fragment=s.fragment and s.context='short'
               where n.context='native' and n.recovered is not null and s.recovered is not null
               group by n.tool,n.ref,n.length,n.model,n.identity'''
    for tool, ref, length, model, identity, n, short, native, gain, loss in db.execute(query):
        pairs.append(dict(tool=tool, reference_id=ref, length=length, model=model, identity=identity,
            paired_targets=n, short_recovered=short, native_recovered=native, gain=gain, loss=loss,
            native_minus_short=(native-short)/n))
    tsv(out/'paired_background.tsv', pairs, ['tool', 'reference_id', 'length', 'model', 'identity',
        'paired_targets', 'short_recovered', 'native_recovered', 'gain', 'loss', 'native_minus_short'])
    tsv(out/'mapping_costs.tsv', costs, ['tool', 'batch_id', 'reference_id', 'context', 'length', 'status',
        'mapping_wall_seconds', 'mapping_cpu_seconds', 'mapping_max_rss_kib',
        'query_preparation_wall_seconds', 'conversion_wall_seconds'])
    tsv(out/'index_costs.tsv', index_costs, ['tool', 'reference_id', 'stage', 'status',
        'gnu_wall_seconds', 'cpu_seconds', 'max_rss_kib', 'wrapper_wall_seconds'])
    macro = {}
    for tool in c['enabled_tools']:
        for context in c['contexts']:
            selected = [x for x in cells if x['tool'] == tool and x['context'] == context and
                        x['model'] != 'P' and x['recovery_fraction'] is not None]
            macro[tool+':'+context] = dict(eligible_cells=len(selected),
                macro_recovery=sum(x['recovery_fraction'] for x in selected)/len(selected) if selected else None)
    write_json(out/'summary.json', dict(status='COMPLETED_DECLARED_MATRIX', source_loci=manifest['source_loci'],
        macro=macro, contexts=c['contexts'], clustered=c['clustered'],
        budget=read_json(root(c)/'budget.json'), source_config_sha256=manifest['config_sha256'],
        metric=c['metrics'], scientific_acceptance='NOT_AUTOMATICALLY_ESTABLISHED',
        uncertainty='Descriptive source-stratified results; no query-independent binomial CI or natural-family independence claim',
        bowtie_limit='Original ungapped Bowtie -v3; queries>1024bp NOT_SUPPORTED, excluded from recovery denominator and explicitly counted'))
    db.close()


def plot(c):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out = root(c)/'summary'
    with (out/'grouped.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    colors = dict(zip(c['enabled_tools'], plt.get_cmap('tab10').colors))
    for model in c['models']:
        fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharey=True, constrained_layout=True)
        for ax, length in zip(axes.flat, c['lengths']):
            for tool in c['enabled_tools']:
                series = defaultdict(list)
                for r in rows:
                    if (r['tool'] == tool and r['model'] == model and r['context'] == 'short' and
                        int(r['length']) == length and r['recovery_fraction']):
                        series[int(r['identity'])].append(float(r['recovery_fraction']))
                ids = sorted(series)
                ax.plot(ids, [sum(series[i])/len(series[i]) for i in ids], 'o-', color=colors[tool], label=tool, markersize=3)
            ax.set(title=f'{length} bp', xlabel='Nominal identity (%)', ylabel='Region90 recovery')
            ax.set_ylim(-.02, 1.02); ax.grid(alpha=.2)
        axes.flat[0].legend(fontsize=8)
        fig.suptitle(f'{model}: isolated queries, equal reference weights')
        for ext in ['png', 'pdf']:
            fig.savefig(out/f'recovery_{model}.{ext}', dpi=180)
        plt.close(fig)
