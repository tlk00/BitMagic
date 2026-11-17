"""Check journal scheduling, parallel timings, hierarchy and emergency replacement."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch
from travel_journal import TravelJournal


class TravelJournalTests(unittest.TestCase):
    """Use a controlled monotonic clock without waiting for publication intervals."""

    def test_schedule_and_emergency(self):
        """First publication waits for preparation; failures replace reports immediately."""
        with tempfile.TemporaryDirectory() as directory, patch('travel_journal.time.monotonic') as clock:
            clock.return_value = 0
            path = Path(directory) / 'travel-journal.txt'
            path.write_text('previous voyage')
            console = Mock()
            journal = TravelJournal(path, console)
            cache = journal.start('Reference cache')
            journal.publish()
            self.assertEqual(path.read_text(), 'previous voyage')
            clock.return_value = 60
            journal.finish(cache)
            journal.publish(first=True)
            parent = journal.start('Compilation')
            first = journal.start('perf', parent, '/logs/build-perf.log')
            second = journal.start('stress', parent, '/logs/build-stress.log')
            clock.return_value = 120
            journal.finish(first)
            clock.return_value = 180
            journal.finish(second)
            journal.finish(parent)
            clock.return_value = 659
            journal.publish()
            self.assertNotIn('Compilation', path.read_text())
            clock.return_value = 660
            journal.publish()
            report = path.read_text()
            self.assertIn('  perf', report)
            self.assertIn('  stress', report)
            self.assertIn('2.0', report)
            console.assert_not_called()
            clock.return_value = 1860
            journal.publish()
            console.assert_called_once()
            stage = journal.start('Tests (Norm)')
            journal.start('stress_release', stage, '/logs/stress/runner.log')
            clock.return_value = 1861
            journal.publish(emergency='SIGSEGV')
            report = path.read_text()
            self.assertIn('DERAILED: SIGSEGV', report)
            self.assertIn('/logs/stress/runner.log', report)
            self.assertIn('FAILED', report)
            self.assertFalse(path.with_suffix('.txt.tmp').exists())


if __name__ == '__main__':
    unittest.main()
