"""Fireman reports based on sampled host CPU and available memory."""
import random
import threading
import time
from pathlib import Path


class Schedule:
    """Share one 15–30 minute clock and CPU baseline across test workers."""

    def __init__(self):
        self.lock = threading.Lock()
        self.deadline = float('inf')
        self.baseline = None

    @staticmethod
    def cpu():
        """Return Linux aggregate CPU counters; unavailable on other platforms."""
        try:
            fields = Path('/proc/stat').read_text().splitlines()[0].split()
            values = [int(x) for x in fields[1:9]]
            if fields[0] != 'cpu' or len(values) < 4:
                return None
            return sum(values), values[3] + (values[4] if len(values) > 4 else 0)
        except (OSError, ValueError, IndexError):
            return None

    def start(self):
        """Start a cycle and establish CPU counter baseline."""
        with self.lock:
            self.baseline = self.cpu()
            self.deadline = time.monotonic() + random.uniform(900, 1800)

    def due(self, memory):
        """Return a factual, figurative report when due; reserve one slot."""
        with self.lock:
            now = time.monotonic()
            if now < self.deadline:
                return None
            current, prior = self.cpu(), self.baseline
            self.baseline = current
            self.deadline = now + random.uniform(900, 1800)
            usage = None
            if current and prior and current[0] > prior[0] and current[1] >= prior[1]:
                usage = max(0.0, min(1.0, 1 - (current[1] - prior[1]) / (current[0] - prior[0])))
            total, free = memory.get('total'), memory.get('available')
            limit, used = memory.get('cgroup_limit'), memory.get('cgroup_current')
            if limit and used is not None and (not total or limit < total):
                total, free = limit, min(free, max(0, limit - used)) if free is not None else max(0, limit - used)
            ratio = free / total if total and free is not None else None
            if usage is None:
                line = 'The fire is steady. Keeping watch on the boiler.'
            elif usage >= .8:
                line = 'We are pulling hard up this grade! I am feeding the fire, but this engine wants every pound of steam we can give her!'
            elif usage >= .4:
                line = 'The engine is working steadily. The fire and steam are holding.'
            else:
                line = 'Easy steaming on this stretch. The fire is steady.'
            if ratio is not None and ratio <= .15:
                line += ' Coal is critically low! Requesting a fuel stop urgently!'
            elif ratio is not None and ratio <= .5:
                line += ' Coal reserves are getting low. Keeping a close watch on our supply.'
            return 'Fireman Sam to Engineer Joe. ' + line
