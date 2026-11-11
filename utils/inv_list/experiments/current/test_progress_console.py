#!/usr/bin/env python3
"""Check immediate captured progress and single-row interactive rendering."""
import os
from pathlib import Path
import pty
import selectors
import subprocess
import unittest

SCRIPT = Path(__file__).resolve().parent / 'progress_console.py'


class ProgressConsoleTest(unittest.TestCase):
    """Exercise the renderer as a streaming subprocess, without dataset I/O."""

    def test_captured_updates_are_immediate(self):
        """Each update must arrive before the next input line or end of input."""
        with subprocess.Popen(['python3', str(SCRIPT)], stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE) as process:
            with selectors.DefaultSelector() as selector:
                selector.register(process.stdout, selectors.EVENT_READ)
                for count in (1, 2):
                    line = (f'Progress compress vectors={count} prepared=10 '
                            f'completed=4 elapsed_s=30 output_bytes=1024\n')
                    process.stdin.write(line.encode())
                    process.stdin.flush()
                    self.assertTrue(selector.select(3), 'progress was retained')
                    text = process.stdout.readline().decode()
                    self.assertIn(f'{count}v', text)
                    self.assertIn('done4 prep10', text)
                    self.assertNotIn('\x1b', text)
                process.stdin.write(b'Warning: test message\n')
                process.stdin.flush()
                self.assertTrue(selector.select(3))
                self.assertEqual(process.stdout.readline(), b'Warning: test message\n')
            process.stdin.close()
            process.wait(timeout=3)
            self.assertEqual(process.returncode, 0)
            self.assertEqual(process.stdout.read(), b'')

    def test_silent_keeps_results_and_errors(self):
        """Silent presentation discards progress while retaining ordinary diagnostic lines."""
        result = subprocess.run(['python3', str(SCRIPT), '--silent'],
                                input='Progress verify vectors=1 elapsed_s=30\nError: visible\nVerification complete.\n',
                                capture_output=True, text=True, check=True)
        self.assertEqual(result.stdout, 'Error: visible\nVerification complete.\n')
        self.assertNotIn('\x1b', result.stdout)

    def test_terminal_updates_reuse_one_row(self):
        """TTY progress clears/replaces its row and leaves ordinary text visible."""
        master, slave = pty.openpty()
        try:
            with subprocess.Popen(['python3', str(SCRIPT)], stdin=subprocess.PIPE,
                                  stdout=slave) as process:
                os.close(slave)
                slave = None
                with selectors.DefaultSelector() as selector:
                    selector.register(master, selectors.EVENT_READ)
                    process.stdin.write(b'Progress compress vectors=1 elapsed_s=30\n')
                    process.stdin.flush()
                    self.assertTrue(selector.select(3))
                    output = os.read(master, 4096)
                    self.assertIn(b'\r\x1b[2Kcompress 1v', output)
                    self.assertNotIn(b'\n', output)
                    process.stdin.write(b'Warning: visible\n')
                    process.stdin.flush()
                    self.assertTrue(selector.select(3))
                    self.assertIn(b'Warning: visible', os.read(master, 4096))
                process.stdin.close()
                process.wait(timeout=3)
                self.assertEqual(process.returncode, 0)
        finally:
            os.close(master)
            if slave is not None:
                os.close(slave)


if __name__ == '__main__':
    unittest.main()
