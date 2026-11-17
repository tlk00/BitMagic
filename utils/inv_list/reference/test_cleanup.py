#!/usr/bin/env python3
"""Check metadata cleanup boundaries and generator failure/lock handling."""
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

from cleanup_reference import cleanup_output, is_appledouble
import gov2_reference


def metadata(path):
    """Write a minimal valid AppleDouble container to path for cleanup fixtures."""
    path.write_bytes(struct.pack('>II16sHIII', 0x00051607, 0x00020000,
                                 b'Mac OS X', 1, 9, 38, 1) + b'x')


class CleanupTests(unittest.TestCase):
    """Preserve payloads, ordinary dotfiles, tools, failed work and symlink targets."""

    def test_cleanup_boundaries(self):
        """Delete only valid metadata, including orphans, and empty staging dirs."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / 'reference'
            root.mkdir()
            outside = Path(temporary) / 'outside'
            outside.mkdir()
            metadata(outside / '._keep')
            (root / 'linked-directory').symlink_to(outside, target_is_directory=True)
            (root / '._linked-file').symlink_to(outside / '._keep')
            (root / '.tools').mkdir()
            (root / '.tools' / 'generator').write_bytes(b'keep executable')
            (root / '.partial').mkdir()
            (root / '.user-settings').write_bytes(b'keep settings')
            (root / '._notes').write_bytes(b'keep ordinary file')
            metadata(root / '._orphan')
            metadata(root / '.tools' / '._generator')
            metadata(root / '.partial' / '._leftover')
            malformed = root / '._truncated'
            metadata(malformed)
            malformed.write_bytes(malformed.read_bytes()[:-1])
            self.assertFalse(is_appledouble(malformed))
            self.assertEqual(cleanup_output(root), (3, 117, 1))
            self.assertEqual((root / '.tools' / 'generator').read_bytes(), b'keep executable')
            self.assertTrue((root / '.user-settings').exists())
            self.assertTrue((root / '._notes').exists())
            self.assertTrue(malformed.exists())
            self.assertTrue((outside / '._keep').exists())
            self.assertTrue((root / '._linked-file').is_symlink())
            self.assertEqual(cleanup_output(root), (0, 0, 0))

    def test_failed_generation(self):
        """Clean lock metadata after failure while retaining incomplete work evidence."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / 'reference'

            def fail(args, sets, levels):
                """Create lock/staging sidecars, then emulate a failed reader."""
                metadata(args.out / '.reference-lock' / '._owner.json')
                metadata(args.out / '._.reference-lock')
                work = args.out / '.partial' / 'failed-case'
                work.mkdir(parents=True)
                (work / 'FAILURE.txt').write_text('reader failed')
                metadata(work / '._FAILURE.txt')
                raise RuntimeError('reader failed')

            with patch.object(gov2_reference, 'execute', side_effect=fail):
                with self.assertRaisesRegex(RuntimeError, 'reader failed'):
                    gov2_reference.main(['--allow-local-output', '--out', str(root)])
            self.assertFalse((root / '.reference-lock').exists())
            self.assertFalse(list(root.rglob('._*')))
            self.assertEqual((root / '.partial' / 'failed-case' / 'FAILURE.txt').read_text(), 'reader failed')


if __name__ == '__main__':
    unittest.main()
