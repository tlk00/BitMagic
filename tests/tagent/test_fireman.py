"""Check fireman reports against deterministic CPU and memory samples."""
import unittest
from unittest.mock import patch
import fireman


class FiremanTests(unittest.TestCase):
    def report(self, baseline, current, memory):
        clock = fireman.Schedule()
        clock.baseline = baseline
        clock.deadline = 0
        with patch.object(clock, 'cpu', return_value=current), patch.object(fireman.time, 'monotonic', return_value=100), patch.object(fireman.random, 'uniform', return_value=900):
            message = clock.due(memory)
            self.assertIsNone(clock.due(memory))
            return message

    def test_high_cpu_and_memory_thresholds(self):
        high = self.report((100, 10), (200, 15), {'total': 100, 'available': 80})
        self.assertIn('pulling hard', high)
        self.assertNotIn('Coal', high)
        low = self.report((100, 10), (200, 100), {'total': 100, 'available': 50})
        self.assertIn('Easy steaming', low)
        self.assertIn('Coal reserves', low)
        critical = self.report((100, 10), (200, 15), {'total': 100, 'available': 15})
        self.assertIn('critically low', critical)

    def test_cgroup_and_unknown_measurements(self):
        message = self.report(None, None, {'total': 1000, 'available': 600, 'cgroup_limit': 100, 'cgroup_current': 90})
        self.assertIn('critically low', message)
        self.assertIn('fire is steady', message)
        message = self.report(None, None, {})
        self.assertNotIn('Coal', message)


if __name__ == '__main__':
    unittest.main()
