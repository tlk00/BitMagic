"""Small-fixture checks for optional inv_list reference caching."""
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import inv_cache


class ReferenceCacheTests(unittest.TestCase):
    """Exercise publication and reuse without real Gov2 data."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'reference'
        self.source.mkdir()
        (self.source / 'Gov2Flat').mkdir()
        (self.source / 'references.tsv').write_text('reference data')
        (self.source / 'Gov2Flat' / 'gov2.sorted').write_bytes(b'input data')
        self.cache = self.root / 'ssd'
        self.cache.mkdir()
        self.messages = []

    def copy(self):
        return inv_cache.copy_reference(self.source, self.cache, 0, self.messages.append)

    def test_verified_copy_and_reuse(self):
        first = self.copy()
        self.assertEqual((first / 'Gov2Flat/gov2.sorted').read_bytes(), b'input data')
        self.assertTrue(any('verify ' in line and 'bytes' in line for line in self.messages))
        with patch.object(inv_cache.shutil, 'disk_usage', side_effect=AssertionError('space check on reuse')):
            self.assertEqual(self.copy(), first)
        self.assertTrue(any('reuse ' in line for line in self.messages))

    def test_changed_source_creates_new_copy(self):
        first = self.copy()
        (self.source / 'references.tsv').write_text('changed reference')
        second = self.copy()
        self.assertNotEqual(first, second)
        self.assertEqual((second / 'references.tsv').read_text(), 'changed reference')

    def test_corrupted_copy_is_quarantined_and_rebuilt(self):
        first = self.copy()
        (first / 'references.tsv').write_bytes(b'x')
        rebuilt = self.copy()
        self.assertEqual(first, rebuilt)
        self.assertEqual((rebuilt / 'references.tsv').read_text(), 'reference data')
        self.assertTrue(any('altered cache preserved' in line for line in self.messages))

    def test_space_failure_does_not_publish(self):
        with patch.object(inv_cache.shutil, 'disk_usage', return_value=type('Usage', (), {'free': 0})()):
            with self.assertRaisesRegex(ValueError, 'insufficient cache space'):
                self.copy()
        self.assertEqual(list((self.cache / 'inv-reference-v1').iterdir()), [])

    def test_copy_failure_leaves_no_partial_cache(self):
        original = Path.open
        def bad_open(path, *args, **kwargs):
            if path.name == 'gov2.sorted' and 'r' in args[0]:
                raise OSError('simulated read error')
            return original(path, *args, **kwargs)
        with patch.object(Path, 'open', bad_open):
            with self.assertRaisesRegex(OSError, 'simulated read error'):
                self.copy()
        self.assertEqual(list((self.cache / 'inv-reference-v1').iterdir()), [])

    def test_digest_mismatch_does_not_publish(self):
        original = Path.open
        def corrupt_before_verification(path, *args, **kwargs):
            if path.name == 'gov2.sorted' and args and args[0] == 'rb' and self.cache.resolve() in path.resolve().parents:
                path.write_bytes(b'other data')
            return original(path, *args, **kwargs)
        with patch.object(Path, 'open', corrupt_before_verification):
            with self.assertRaisesRegex(ValueError, 'verification'):
                self.copy()
        self.assertEqual(list((self.cache / 'inv-reference-v1').iterdir()), [])

    def test_reuse_does_not_read_reference_contents(self):
        first = self.copy()
        original = Path.open
        def guard_data(path, *args, **kwargs):
            if path.name in ('gov2.sorted', 'references.tsv'):
                raise AssertionError('reference content was reread')
            return original(path, *args, **kwargs)
        with patch.object(Path, 'open', guard_data):
            self.assertEqual(self.copy(), first)

    def test_reject_symlink_and_overlapping_locations(self):
        os.symlink('references.tsv', self.source / 'link')
        with self.assertRaisesRegex(ValueError, 'symlink'):
            self.copy()
        with self.assertRaisesRegex(ValueError, 'must not contain'):
            inv_cache.copy_reference(self.source, self.root, 0, self.messages.append)
