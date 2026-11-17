"""In-memory, thread-safe voyage ledger; never reconstructs events from logs."""
import datetime
import threading
import time
from pathlib import Path


class TravelJournal:
    """Own timed voyage entries and atomically publish an English operator report."""

    def __init__(self, path, console):
        """Set the output path and console callback without touching an older report."""
        self.path = Path(path)
        self.console = console
        self.lock = threading.RLock()
        self.rows = []
        self.enabled = False
        self.commit = ''
        self.diagnostic = ''
        self.next_file = float('inf')
        self.next_console = float('inf')

    def start(self, name, parent=None, log=None):
        """Record a running entry and return its identity; parent controls indentation."""
        with self.lock:
            row = dict(name=name, parent=parent, log=str(log or ''), start=time.monotonic(),
                       end=None, finished=None, status='RUNNING')
            self.rows.append(row)
            return row

    def finish(self, row, status='PASS'):
        """Complete an entry exactly once using monotonic elapsed time."""
        with self.lock:
            if row['end'] is None:
                row.update(end=time.monotonic(), finished=datetime.datetime.now().astimezone(), status=status)

    def render(self):
        """Render a consistent hierarchical snapshot including active elapsed times."""
        with self.lock:
            now = time.monotonic()
            active = [r for r in self.rows if r['end'] is None]
            current = '; '.join('{} ({:.1f} min)'.format(r['name'], (now-r['start'])/60) for r in active)
            lines = ['BitMagic TAgent - Travel Journal',
                     'Generated: ' + datetime.datetime.now().astimezone().strftime('%m-%d %H:%M:%S %Z'),
                     'Commit: ' + self.commit,
                     ('DERAILED: ' + self.diagnostic) if self.diagnostic else 'Current: ' + (current or 'Waiting / completed'),
                     '', '{:<48} {:<12} {:<16} {:>10}'.format('Phase / test', 'Status', 'Finished', 'Minutes')]
            def append(parent, depth):
                """Render children beneath their parent in recorded order."""
                for row in self.rows:
                    if row['parent'] is not parent:
                        continue
                    label = '  '*depth + row['name']
                    finished = row['finished'].strftime('%m-%d %H:%M:%S') if row['finished'] else '-'
                    elapsed = ((row['end'] if row['end'] is not None else now)-row['start'])/60
                    lines.append('{:<48} {:<12} {:<16} {:>10.1f}'.format(label, row['status'], finished, elapsed))
                    append(row, depth+1)
            append(None, 0)
            return '\n'.join(lines)+'\n'

    def publish(self, first=False, emergency=None, force=False):
        """Replace the previous complete report; print only on the 30-minute schedule."""
        with self.lock:
            now = time.monotonic()
            if emergency is not None and not self.diagnostic:
                active = [r for r in self.rows if r['end'] is None]
                self.diagnostic = str(emergency) + ' | phase: ' + ' / '.join(r['name'] for r in active)
                logs = [r['log'] for r in active if r['log']]
                if logs:
                    self.diagnostic += ' | logs: ' + ', '.join(logs)
                for row in self.rows:
                    if row['end'] is None:
                        self.finish(row, 'FAILED')
            if not (self.enabled or first or emergency is not None):
                return
            if first:
                self.enabled = True
                self.next_console = now+1800
            if not (first or force or emergency is not None or now >= self.next_file):
                return
            report = self.render()
            temporary = self.path.with_suffix('.txt.tmp')
            try:
                temporary.write_text(report)
                temporary.replace(self.path)
            except OSError as error:
                self.console('TRAVEL JOURNAL ERROR | {} | {}'.format(self.path, error))
                self.next_file = now+60
                return
            self.next_file = now+600
            if now >= self.next_console:
                self.console(report)
                self.next_console = now+1800
