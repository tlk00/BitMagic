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
    def test_perf_groups_are_attributed_to_owned_runner(self):
        """Only group processes in tracked runner process groups are reported."""
        with tempfile.TemporaryDirectory() as root:
            monitor = memory_monitor.MemoryMonitor(Agent(root))
            output = ('101 100 2048 /tmp/perf_release --group sparse-vector\n'
                      '102 100 4096 /tmp/perf_release --group bit-core --fast\n'
                      '103 200 8192 /tmp/perf_release --group unrelated\n'
                      '104 100 100 bash runner\n')
            with patch.object(memory_monitor.subprocess, 'check_output', return_value=output):
                groups = monitor.perf_group_memory({100: 'perf_release'})
            self.assertEqual([(g['variant'], g['group'], g['rss_bytes']) for g in groups],
                             [('perf_release', 'bit-core', 4096 * 1024),
                              ('perf_release', 'sparse-vector', 2048 * 1024)])

    def test_perf_samples_include_group_observations(self):
        """Samples attach group memory for currently tracked perf runners."""
        with tempfile.TemporaryDirectory() as root:
            agent = Agent(root)
            agent.memory_phase = 'perf'
            runner = unittest.mock.Mock()
            runner.pid = 100
            agent.perf_logs[runner] = Path(root) / 'perf_release' / 'perf.log'
            monitor = memory_monitor.MemoryMonitor(agent)
            group = {'pid': 101, 'variant': 'perf_release', 'group': 'bit-core', 'rss_bytes': 2048}
            with patch.object(memory_monitor.platform, 'system', return_value='Darwin'), \
                 patch.object(memory_monitor.subprocess, 'check_output', return_value='34359738368\n'), \
                 patch.object(monitor, 'perf_group_memory', return_value=[group]) as read:
                sample = monitor.sample()
            read.assert_called_once_with({100: 'perf_release'})
            self.assertEqual(sample['perf_groups'], [group])

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
