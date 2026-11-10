#!/usr/bin/env python3
import json
import struct
import tempfile
import unittest
from pathlib import Path
from sample_gov2 import sample, records

class SamplingTest(unittest.TestCase):
    def test_paired_midpoints_preserve_records(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=root/'input';source.mkdir()
            expected={}
            for ds in ('sorted','unsorted'):
                rows=[]
                for i in range(33):
                    ids=[i*100+j+(ds=='unsorted') for j in range(i+1)]
                    rows.append(struct.pack('<I',len(ids))+struct.pack('<'+'I'*len(ids),*ids))
                expected[ds]=rows
                (source/('gov2.'+ds)).write_bytes(b''.join(rows))
            meta=sample(source,root/'out',5)
            indices=[3,9,16,23,29]
            self.assertEqual(meta['indices_zero_based'],indices)
            for ds in expected:
                self.assertEqual((root/'out'/('gov2.'+ds)).read_bytes(),b''.join(expected[ds][i] for i in indices))
            self.assertEqual(json.loads((root/'out/sample_manifest.json').read_text()),meta)
            with self.assertRaises(ValueError):sample(source,root/'too_many',34)
            meta2=sample(source,root/'all',33)
            self.assertEqual(meta2['indices_zero_based'],list(range(33)))
            for ds in expected:self.assertEqual((root/'all'/('gov2.'+ds)).read_bytes(),(source/('gov2.'+ds)).read_bytes())
    def test_truncated_and_empty_records(self):
        import io
        for data in (b'\x01',struct.pack('<I',2)+struct.pack('<I',3),struct.pack('<I',0)):
            with self.assertRaises(ValueError):list(records(io.BytesIO(data),len(data)))

if __name__=='__main__': unittest.main()
