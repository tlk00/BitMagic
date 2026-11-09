"""Original engineer chatter and a thread-safe, nonrepeating test-phase clock."""
from pathlib import Path
import random
import threading
import time


class Schedule:
    """Shuffle the complete joke deck before reusing any entry."""

    def __init__(self):
        """Keep the deck across test phases and build cycles."""
        self.lines = Path(__file__).with_suffix('.txt').read_text().splitlines()
        self.deck = []
        self.deadline = float('inf')
        self.count = 0
        self.lock = threading.Lock()

    def start(self):
        """Give each new test phase ten to twenty minutes before chatter."""
        with self.lock:
            self.deadline = time.monotonic() + random.uniform(600, 1200)

    def due(self):
        """Return one due line, reserving the next slot across worker threads."""
        with self.lock:
            if time.monotonic() < self.deadline:
                return None
            self.deadline = time.monotonic() + random.uniform(600, 1200)
            if not self.deck:
                self.deck = list(self.lines)
                random.shuffle(self.deck)
            prefix = 'Train engineer Joe here. ' if self.count % 4 == 0 else ''
            self.count += 1
            return prefix + self.deck.pop()
