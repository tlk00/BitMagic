"""Focused memory-monitor threshold and snapshot checks."""
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch
import memory_monitor


class Agent:
    def __init__(self, root):
        self.home = Path(root)
        self.active_run = None
        self.children = set()
        self.perf_logs = {}
        self.mutex = threading.RLock()
        self.jobs = 2
        self.jemalloc = None
        self.messages = []

    def say(self, message):
        self.messages.append(message)


class MonitorTests(unittest.TestCase):
    def test_thresholds_escalation_and_oom(self):
        with tempfile.TemporaryDirectory() as root:
            monitor = memory_monitor.MemoryMonitor(Agent(root))
            gib = 1024 ** 3
            data = {'total': 16 * gib, 'available': 4 * gib, 'cgroup_events': {'oom_kill': '0'}}
            with patch.object(monitor, 'sample', side_effect=lambda: dict(data)), patch.object(monitor, 'snapshot') as snap:
                monitor.tick()
                self.assertEqual(snap.call_count, 0)
                data['available'] = 2 * gib
                monitor.tick()
                self.assertEqual(snap.call_args.args[0], 'warning')
                monitor.last_snapshot = memory_monitor.time.monotonic()
                monitor.tick()
                self.assertEqual(snap.call_count, 1)
                data['available'] = 500 * 1024**2
                monitor.tick()
                self.assertEqual(snap.call_args.args[0], 'critical')
                data['cgroup_events'] = {'oom_kill': '1'}
                monitor.tick()
                self.assertEqual(snap.call_args.args[0], 'oom-counter')
            self.assertEqual(len(Path(root, 'resources.jsonl').read_text().splitlines()), 5)

    def test_cgroup_headroom_and_snapshot(self):
        with tempfile.TemporaryDirectory() as root:
            agent = Agent(root)
            monitor = memory_monitor.MemoryMonitor(agent)
            gib = 1024**3
            data = {'total': 32*gib, 'available': 20*gib, 'cgroup_limit': 4*gib,
                    'cgroup_current': 3*gib, 'cgroup_events': {}}
            with patch.object(monitor, 'sample', return_value=data):
                monitor.tick()
            directory = next((Path(root) / 'diagnostics').iterdir())
            self.assertTrue((directory / 'processes.txt').exists())
            self.assertEqual(json.loads((directory / 'summary.json').read_text())['cgroup_limit'], 4*gib)
            self.assertTrue(any('WARNING' in x for x in agent.messages))

    def test_macos_unavailable_values_are_explicit(self):
        with tempfile.TemporaryDirectory() as root:
            monitor = memory_monitor.MemoryMonitor(Agent(root))
            with patch.object(memory_monitor.platform, 'system', return_value='Darwin'), \
                 patch.object(memory_monitor.subprocess, 'check_output', return_value='34359738368\n'):
                sample = monitor.sample()
            self.assertIsNone(sample['available'])
            self.assertIn('unavailable', sample['note'])


if __name__ == '__main__':
    unittest.main()
