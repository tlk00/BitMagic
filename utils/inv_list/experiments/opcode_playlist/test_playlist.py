#!/usr/bin/env python3
import functools,subprocess,tempfile,unittest
from pathlib import Path
from playlist import rank_compare,byte_rank_key
TOOL=Path(__file__).resolve().parents[2]/'gov2_playlist'
class PlaylistTest(unittest.TestCase):
    def test_byte_contribution_ranking(self):
        rows=[dict(bytes=100,ids=1,vector_index=0),dict(bytes=200,ids=100,vector_index=7),dict(bytes=200,ids=200,vector_index=3)]
        self.assertEqual([r['vector_index'] for r in sorted(rows,key=byte_rank_key)],[3,7,0])
    def test_exact_ranking_and_ties(self):
        a=dict(bytes=2**60+1,ids=2**60,vector_index=5)
        b=dict(bytes=2**60,ids=2**60,vector_index=0)
        self.assertLess(rank_compare(a,b),0)
        rows=[dict(bytes=100,ids=10,vector_index=8),dict(bytes=200,ids=20,vector_index=7),dict(bytes=200,ids=20,vector_index=3)]
        self.assertEqual([r['vector_index'] for r in sorted(rows,key=functools.cmp_to_key(rank_compare))],[3,7,8])
    def test_membership_and_reject_bad_streams(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);txt=p/'indices';bv=p/'list.bv';decoded=p/'decoded'
            txt.write_text('1107204\n0\n65536\n7\n')
            subprocess.run([str(TOOL),'pack',str(txt),str(bv)],check=True)
            subprocess.run([str(TOOL),'unpack',str(bv),str(decoded)],check=True)
            self.assertEqual(decoded.read_text(),'0\n7\n65536\n1107204\n')
            good=bv.read_bytes()
            for data in (good[:-1],good+b'junk'):
                bv.write_bytes(data)
                self.assertNotEqual(subprocess.run([str(TOOL),'unpack',str(bv),str(decoded)]).returncode,0)
            txt.write_text('1\n1\n')
            self.assertNotEqual(subprocess.run([str(TOOL),'pack',str(txt),str(bv)]).returncode,0)
if __name__=='__main__':unittest.main()
