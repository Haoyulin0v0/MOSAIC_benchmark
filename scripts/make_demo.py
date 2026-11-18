#!/usr/bin/env python3
"""Small synthetic interface/workflow fixture; never biological benchmark data."""
import argparse
import json
from pathlib import Path
import random
from mosaic_benchmark.common import write_fasta, write_json


def make_demo(directory, tools=None):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    config = json.loads((Path(__file__).resolve().parents[1]/'config/benchmark.json').read_text())
    rng = random.Random(271828)
    write_fasta(directory/'reference.fa', {'chr_one': ''.join(rng.choices('ACGT', k=120000)),
                                          'chr_two': ''.join(rng.choices('ACGT', k=80000))})
    write_fasta(directory/'donor.fa', {'unrelated_synthetic': ''.join(rng.choices('ACGT', k=50000))})
    config.update(profile='synthetic_interface_fixture', work_dir=str(directory/'run'),
        references=[dict(id='demo', fasta=str(directory/'reference.fa'))],
        identities=[70, 100], loci_per_reference_length=4)
    config['background'].update(donor=dict(id='synthetic', fasta=str(directory/'donor.fa')),
                               selected_loci_per_reference_length=4, identities=[70])
    config['resources'].update(threads=4, memory_gib=8, index_timeout_seconds=60,
        mapping_timeout_seconds=30, conversion_timeout_seconds=30, total_native_budget_seconds=7200)
    if tools is not None:
        config['enabled_tools'] = tools
    write_json(directory/'config.json', config)
    return directory/'config.json'


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', default='work/demo')
    p.add_argument('--tools', nargs='+')
    a = p.parse_args()
    print(make_demo(a.directory, a.tools))
