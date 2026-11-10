"""Bounded, best-effort memory diagnostics for TAgent (no external packages)."""
import datetime
import json
import os
from pathlib import Path
import platform
import subprocess
import threading
import time


class MemoryMonitor:
    """Sample memory independently and preserve evidence before exhaustion."""

    def __init__(self, agent, interval=5):
        self.agent = agent
        self.interval = interval
        self.stop_event = threading.Event()
        self.lock = threading.Lock()
        self.level = 0
        self.last_snapshot = 0
        self.last_oom = None
        self.thread = None

    def sample(self):
        """Read host-visible memory and, where available, cgroup v2 limits."""
        result = {'time': datetime.datetime.now().astimezone().isoformat(),
                  'platform': platform.system(), 'cycle': str(self.agent.active_run or ''),
                  'workers': len(self.agent.children), 'phase': getattr(self.agent, 'memory_phase', 'unknown')}
        if platform.system() == 'Linux':
            values = {}
            for line in Path('/proc/meminfo').read_text().splitlines():
                parts = line.split()
                if len(parts) >= 2 and parts[1].isdigit():
                    values[parts[0].rstrip(':')] = int(parts[1]) * 1024
            result.update(total=values.get('MemTotal'), available=values.get('MemAvailable'),
                          swap_total=values.get('SwapTotal'), swap_free=values.get('SwapFree'))
            try:
                result['pressure'] = Path('/proc/pressure/memory').read_text().strip()
            except OSError:
                pass
            try:
                rel = next(line.split(':', 2)[2] for line in Path('/proc/self/cgroup').read_text().splitlines()
                           if line.startswith('0::'))
                group = Path('/sys/fs/cgroup') / rel.lstrip('/')
                limit = (group / 'memory.max').read_text().strip()
                current = int((group / 'memory.current').read_text().strip())
                events = dict(line.split() for line in (group / 'memory.events').read_text().splitlines())
                result.update(cgroup_current=current, cgroup_limit=None if limit == 'max' else int(limit),
                              cgroup_events=events)
            except (OSError, ValueError, StopIteration):
                pass
        elif platform.system() == 'Darwin':
            total = subprocess.check_output(['sysctl', '-n', 'hw.memsize'], text=True, timeout=5)
            result.update(total=int(total.strip()), available=None, swap_total=None, swap_free=None,
                          note='macOS available-memory and cgroup metrics unavailable')
        else:
            result['note'] = 'memory metrics unavailable on this platform'
        return result

    def history(self, sample):
        """Append flushed samples, rotating at 2 MiB with one previous file."""
        path = self.agent.home / 'resources.jsonl'
        if path.exists() and path.stat().st_size >= 2 * 1024 * 1024:
            path.replace(path.with_suffix('.jsonl.1'))
        with path.open('a') as stream:
            stream.write(json.dumps(sample) + '\n')
            stream.flush()

    def snapshot(self, reason, sample=None):
        """Write a bounded process snapshot; errors never disrupt testing."""
        try:
            with self.lock:
                sample = sample or self.sample()
                base = self.agent.active_run or self.agent.home
                directory = base / 'diagnostics' / ('memory-' + datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f') + '-' + reason)
                directory.mkdir(parents=True, exist_ok=True)
                (directory / 'reason.txt').write_text('Trigger: ' + reason + '\n')
                (directory / 'summary.json').write_text(json.dumps(sample, indent=2) + '\n')
                with self.agent.mutex:
                    children = [{'pid': p.pid, 'args': p.args} for p in self.agent.children]
                    logs = {str(p.pid): str(log) for p, log in self.agent.perf_logs.items()}
                (directory / 'owned-jobs.json').write_text(json.dumps({'children': children, 'perf_logs': logs}, indent=2, default=str))
                try:
                    ps = subprocess.run(['ps', '-eo', 'pid,ppid,pgid,rss,args'], capture_output=True, text=True, timeout=8)
                    lines = ps.stdout.splitlines() or ['Process listing unavailable: ' + ps.stderr.strip()]
                except (OSError, subprocess.SubprocessError) as error:
                    lines = ['Process listing unavailable: ' + str(error)]
                rows = sorted(lines[1:], key=lambda x: int(x.split(None, 4)[3]) if len(x.split(None, 4)) > 3 and x.split(None, 4)[3].isdigit() else 0, reverse=True)
                (directory / 'processes.txt').write_text('\n'.join(lines[:1] + rows[:100]) + '\n')
                if platform.system() == 'Linux':
                    for name in ('meminfo', 'swaps'):
                        try:
                            (directory / (name + '.txt')).write_text(Path('/proc', name).read_text())
                        except OSError:
                            pass
                    # VmSwap is not supplied by ps; obtain it for the ranked processes.
                    swaps = []
                    for row in rows[:100]:
                        try:
                            pid = int(row.split(None, 1)[0])
                            status = Path('/proc', str(pid), 'status').read_text()
                            value = next((x for x in status.splitlines() if x.startswith('VmSwap:')), 'VmSwap: unavailable')
                            swaps.append('{} {}'.format(pid, value))
                        except (OSError, ValueError, IndexError):
                            continue
                    (directory / 'process-swap.txt').write_text('\n'.join(swaps) + '\n')
                for name in ('resources.jsonl',):
                    path = self.agent.home / name
                    if path.exists():
                        (directory / name).write_bytes(path.read_bytes())
                available = sample.get('available')
                free_gib = '{:.2f} GiB'.format(available / 1024**3) if available is not None else 'unknown'
                swap_total, swap_free = sample.get('swap_total'), sample.get('swap_free')
                swap = '{:.2f}/{:.2f} GiB'.format((swap_total - swap_free) / 1024**3, swap_total / 1024**3) if swap_total is not None and swap_free is not None else 'unknown'
                self.agent.say('MEMORY {} | available {} | swap used {} | evidence: {}'.format(reason.upper(), free_gib, swap, directory))
                self.last_snapshot = time.monotonic()
                snapshots = sorted((base / 'diagnostics').glob('memory-*'))
                for old in snapshots[:-20]:
                    if old.is_dir() and not old.is_symlink():
                        import shutil
                        shutil.rmtree(old)
        except Exception as error:
            try:
                self.agent.say('Memory snapshot unavailable: ' + str(error))
            except Exception:
                pass

    def tick(self):
        sample = self.sample()
        self.history(sample)
        total, available = sample.get('total'), sample.get('available')
        limit, current = sample.get('cgroup_limit'), sample.get('cgroup_current')
        if limit is not None and current is not None and (total is None or limit < total):
            total, available = limit, max(0, limit - current)
        level = 0
        if total and available is not None:
            if available < max(total * .05, 768 * 1024**2):
                level = 2
            elif available < max(total * .15, 2 * 1024**3):
                level = 1
            # Hysteresis prevents repeated crossings at a threshold.
            elif self.level and available < max(total * .18, 2.5 * 1024**3):
                level = self.level
        oom = int(sample.get('cgroup_events', {}).get('oom_kill', 0))
        increased = self.last_oom is not None and oom > self.last_oom
        self.last_oom = oom
        if increased:
            self.snapshot('oom-counter', sample)
        elif level > self.level or (level and time.monotonic() - self.last_snapshot >= 300):
            self.snapshot('critical' if level == 2 else 'warning', sample)
        self.level = level
        return sample

    def start(self):
        try:
            sample = self.tick()
            self.agent.say('MEMORY startup | total={} available={} swap_total={} swap_free={} cgroup_limit={} | workers={} allocator={}'.format(
                sample.get('total'), sample.get('available'), sample.get('swap_total'), sample.get('swap_free'),
                sample.get('cgroup_limit'), self.agent.jobs, self.agent.jemalloc or 'system'))
        except Exception as error:
            self.agent.say('Memory monitoring unavailable: ' + str(error))
        self.thread = threading.Thread(target=self.run, name='memory-monitor', daemon=True)
        self.thread.start()

    def run(self):
        while not self.stop_event.wait(self.interval):
            try:
                self.tick()
            except Exception as error:
                try:
                    self.agent.say('Memory monitoring unavailable: ' + str(error))
                except Exception:
                    pass
                # Back off persistent operational errors.
                self.stop_event.wait(60)

    def close(self):
        self.stop_event.set()
        if self.thread:
            self.thread.join(timeout=2)
