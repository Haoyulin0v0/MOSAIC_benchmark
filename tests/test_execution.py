from pathlib import Path
import sys
import tempfile
import unittest
from mosaic_benchmark.common import read_json
from mosaic_benchmark.native import measured


class ExecutionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.base = Path(self.temp.name)
        self.config = dict(work_dir=str(self.base), enabled_tools=[], resources=dict(
            memory_gib=1, threads=1, output_limit_mib=1, index_output_limit_mib=1,
            total_native_budget_seconds=60))

    def tearDown(self):
        self.temp.cleanup()

    def test_success_and_timeout_keep_evidence_and_charge_budget(self):
        a = measured(self.config, [sys.executable, '-c', 'print("ok")'], self.base/'ok', 5)
        self.assertEqual(a['status'], 'SUCCESS')
        self.assertEqual(Path(a['raw_path']).read_text().strip(), 'ok')
        b = measured(self.config, [sys.executable, '-c', 'import time; time.sleep(10)'], self.base/'timeout', .2)
        self.assertEqual(b['status'], 'TIMEOUT')
        self.assertTrue(Path(b['raw_path']).exists())
        budget = read_json(self.base/'budget.json')
        self.assertEqual(budget['completed_attempts'], 2)
        self.assertGreater(budget['spent_seconds'], .2)
        self.assertNotIn('reserved_attempt', budget)

    def test_budget_rejects_before_launch_and_keeps_prior_attempts(self):
        self.config['resources']['total_native_budget_seconds'] = .1
        with self.assertRaises(RuntimeError):
            measured(self.config, [sys.executable, '-c', 'raise Exception()'], self.base/'rejected', 5)
        self.assertFalse((self.base/'rejected/attempt_0001/receipt.json').exists())


if __name__ == '__main__':
    unittest.main()
