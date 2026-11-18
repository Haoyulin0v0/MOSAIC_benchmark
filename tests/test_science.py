import json
from pathlib import Path
import sys
import tempfile
import unittest

from mosaic_benchmark.alignment import Alignment, parse
from mosaic_benchmark.common import load_config, random_stream, rc, read_fasta, read_json, write_fasta
from mosaic_benchmark.dataset import prepare, batches, batch_paths, cases, verify_dataset
from mosaic_benchmark.evaluate import score
from mosaic_benchmark.fasta import Fasta
from mosaic_benchmark.models import mutate

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
from make_demo import make_demo


class ScienceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_seed_domains_do_not_collide(self):
        a, h = random_stream(123, 'source', 'ab', 'c')
        b, hh = random_stream(123, 'source', 'a', 'bc')
        self.assertNotEqual(h, hh)
        self.assertEqual(a.random(), random_stream(123, 'source', 'ab', 'c')[0].random())

    def test_all_edit_budgets_and_replay(self):
        rng, _ = random_stream(23, 'fixture')
        for length in [150, 250, 500, 1000]:
            source = ''.join(rng.choices('ACGT', k=length))
            for identity in list(range(50, 91, 5))+[92, 94, 96, 98, 100]:
                for model in ['U', 'C', 'I']:
                    with self.subTest(length=length, identity=identity, model=model):
                        query, path = mutate(source, identity, model, 123, 'fixture')
                        self.assertEqual(path['qaln'].replace('-', ''), query)
                        self.assertEqual(path['taln'].replace('-', ''), source)
                        self.assertEqual(path['M']+path['X']+path['D'], length)
                        self.assertEqual(path['M']+path['X']+path['I'], len(query))
                        self.assertEqual(path['M'], sum(q == t for q, t in zip(path['qaln'], path['taln'])))
                        if model == 'C' and identity != 100:
                            self.assertEqual(len(path['hot_positions']), round(length*.5))
                            self.assertEqual(len(set(path['changed_positions']) & set(path['hot_positions'])), round(.8*path['X']))
                        self.assertEqual((query, path), mutate(source, identity, model, 123, 'fixture'))

    def test_float_rounding_is_preserved(self):
        _, p = mutate('A'*150, 85, 'U', 123, 'rounding')
        self.assertEqual(p['X'], 22)

    def test_seekable_fasta_and_mask(self):
        write_fasta(self.base/'ref.fa', {'a': 'acGT'*51, 'b': 'N'*11+'ACGT'})
        with Fasta(self.base/'ref.fa', self.base/'cache.json') as f:
            self.assertEqual(f.fetch('a', 77, 123), ('acGT'*51)[77:123])
            self.assertEqual(f.fetch('b', 10, 15), 'NACGT')
        with Fasta(self.base/'ref.fa', self.base/'cache.json') as f:
            self.assertEqual(f.fetch('a', 203, 204), 'T')

    def test_distinct_loci_and_paired_background(self):
        config_path = make_demo(self.base/'demo', [])
        c = load_config(config_path)
        # The scheduler creates output parents before running the generation rule.
        (Path(c['work_dir'])/'dataset').mkdir(parents=True)
        prepare(c)
        source_rows = [json.loads(x) for x in (Path(c['work_dir'])/'dataset/sources.jsonl').read_text().splitlines()]
        self.assertEqual(len(source_rows), 16)
        self.assertEqual(len({(x['record'], x['start']) for x in source_rows}), 16)
        for a in source_rows:
            for b in source_rows:
                if a['source_id'] != b['source_id'] and a['record'] == b['record']:
                    self.assertFalse(a['start'] < b['end'] and b['start'] < a['end'])
        short = {}
        for b in batches(c):
            qp, tp = batch_paths(c, b)
            seqs = read_fasta(qp)
            # Original Bowtie reads only the first sequence line: never wrap queries.
            self.assertEqual(len(qp.read_text().splitlines()), 2*len(seqs))
            truths = [json.loads(x) for x in tp.read_text().splitlines()]
            if b['context'] == 'short':
                short.update({t['fragment_id']: seqs[t['query_id']] for t in truths})
        native_samples = 0
        for b in batches(c):
            if b['context'] == 'native':
                qp, tp = batch_paths(c, b); seqs = read_fasta(qp)
                truths = [json.loads(x) for x in tp.read_text().splitlines()]
                native_samples += len(seqs)
                self.assertEqual(len(truths), 4*len(seqs))
                for t in truths:
                    self.assertEqual(seqs[t['query_id']][t['qstart']:t['qend']], short[t['fragment_id']])
        self.assertEqual(native_samples, 6)
        verify_dataset(c, full=True)
        self.assertEqual(sum(b['records'] for b in batches(c) if b['model'] == 'P'), 16)

    def test_dataset_tampering_is_detected(self):
        c = load_config(make_demo(self.base/'demo', [])); prepare(c)
        qp, _ = batch_paths(c, next(iter(batches(c))))
        qp.write_text(qp.read_text()+'\n')
        with self.assertRaises(RuntimeError):
            verify_dataset(c, full=True)

    def test_production_background_covers_six_identities_and_full_locus_cohorts(self):
        config_path = Path(__file__).resolve().parents[1] / 'config/benchmark.json'
        config = load_config(config_path)
        native = [case for case in cases(config) if case['context'] == 'native']
        self.assertEqual({case['identity'] for case in native}, {50, 60, 70, 80, 90, 94})
        self.assertEqual(len(native), 216)
        self.assertTrue(all(case['records'] == 250 for case in native))
        self.assertEqual(sum(case['records'] for case in native), 54000)
        original = read_json(config_path)
        for identities, count in [([60, 80, 94], 1000), ([50, 60, 70, 80, 90, 94], 200)]:
            with self.subTest(identities=identities, count=count):
                original['background'].update(identities=identities,
                                              selected_loci_per_reference_length=count)
                path = self.base / 'invalid.json'
                path.write_text(json.dumps(original))
                with self.assertRaisesRegex(ValueError, 'at least six identity levels'):
                    load_config(path)

    def test_region_clipping_whole_hsp_and_wrong_strand(self):
        t = dict(query_id='q', record='r', strand='+', qstart=50, qend=150,
                 start=1000, end=1100, qaln='A'*100, taln='A'*100)
        # A long HSP can pass m150 even though its clipped source region is100bp.
        a = Alignment('q', 'r', 0, 200, 950, 1150, '+', [(200, 'M')])
        self.assertEqual(score(t, [a], minimum=150)['region90'], 1)
        short = Alignment('q', 'r', 50, 150, 1000, 1100, '+', [(100, 'M')])
        self.assertEqual(score(t, [short], minimum=150)['region90'], 0)
        wrong = Alignment('q', 'r', 50, 150, 1000, 1100, '-', [(100, 'M')])
        self.assertEqual(score(t, [wrong], minimum=50)['region90'], 0)

    def test_cannot_combine_small_hsps_to_pass_gate(self):
        t = dict(query_id='q', record='r', strand='+', qstart=0, qend=150,
                 start=0, end=150, qaln='A'*150, taln='A'*150)
        alignments = [Alignment('q', 'r', 0, 75, 0, 75, '+', [(75, 'M')]),
                      Alignment('q', 'r', 75, 150, 75, 150, '+', [(75, 'M')])]
        self.assertEqual(score(t, alignments, minimum=100)['region90'], 0)
        self.assertEqual(score(t, alignments, minimum=50)['region90'], 1)

    def test_reverse_indel_truth_oracle(self):
        t = dict(query_id='q', record='r', strand='-', qstart=10, qend=15,
                 start=100, end=105, qaln='AC-GTA', taln='ACCG-A')
        a = Alignment('q', 'r', 10, 15, 100, 105, '-', [(2, 'M'), (1, 'D'), (1, 'M'), (1, 'I'), (1, 'M')])
        result = score(t, [a], minimum=1, coverage=.8)
        self.assertEqual(result['strict_generating_pair_recall'], 1)
        self.assertEqual(result['query_coverage'], .8)
        self.assertEqual(result['source_coverage'], .8)

    def test_formats_agree_on_reverse_indel_coordinates(self):
        examples = {
            'paf': 'q\t5\t0\t5\t-\tr\t10\t2\t7\t4\t6\t60\tcg:Z:2M1D1M1I1M\n',
            'sam': 'q\t16\tr\t3\t60\t2M1D1M1I1M\t*\t0\t0\tACGTA\t*\n',
            'maf': 'a score=1\ns r 2 5 + 10 ACCG-A\ns q 0 5 - 5 AC-GTA\n\n',
            'blast': 'q\tr\t5\t1\t3\t7\tAC-GTA\tACCG-A\n',
        }
        normalized = []
        for fmt, text in examples.items():
            p = self.base/fmt; p.write_text(text)
            a = list(parse(p, fmt, {'q': 5}, {'r': 10}))[0]
            normalized.append(list(a.segments()))
        self.assertTrue(all(x == normalized[0] for x in normalized))

    def test_reject_bbox_only_paf_and_bad_cigar(self):
        p = self.base/'bad'; p.write_text('q\t5\t0\t5\t+\tr\t10\t2\t7\t5\t5\t60\n')
        with self.assertRaises(ValueError):
            list(parse(p, 'paf', {'q': 5}, {'r': 10}))
        with self.assertRaises(ValueError):
            list(Alignment('q', 'r', 0, 5, 2, 7, '+', [(4, 'M')]).segments())


if __name__ == '__main__':
    unittest.main()
