"""Python access to compressed inv_list sidecars through the BM C++ codec.

Builds a temporary scalar inspection tool on first use unless INV_INDEX_TOOL
names an existing executable. Decimal boundaries cross subprocess pipes only;
sidecar files always use the compressed V2 format.
"""
from functools import lru_cache
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
_tools = []


@lru_cache(maxsize=1)
def index_tool():
    """Return the configured tool or build one, retaining its temporary directory."""
    override = os.environ.get('INV_INDEX_TOOL')
    if override:
        return Path(override)
    temporary = tempfile.TemporaryDirectory(prefix='bm-index-tool-')
    _tools.append(temporary)
    tool = Path(temporary.name) / 'inv_index'
    subprocess.run([*shlex.split(os.environ.get('CXX', 'c++')), '-std=c++17', '-O2',
                    '-I' + str(HERE.parents[2] / 'src'), str(HERE / 'inv_index.cpp'),
                    '-o', str(tool)], check=True)
    return tool


def read_offsets(path, file_bytes=None):
    """Return all decoded boundaries, checking the associated BLOB's exact size.

    path names a compressed sidecar; file_bytes overrides the size of the .bv
    sibling for synthetic tests. The returned list owns its decoded values.
    """
    path = Path(path)
    if file_bytes is None:
        file_bytes = path.with_suffix('.bv').stat().st_size
    output = subprocess.check_output([str(index_tool()), 'dump', str(path), str(file_bytes)], text=True)
    return list(map(int, output.splitlines()))


def write_offsets(path, offsets, level=6):
    """Replace path with a compressed sidecar for borrowed offsets at the given level."""
    subprocess.run([str(index_tool()), 'write', str(path), str(level)],
                   input=''.join(f'{value}\n' for value in offsets), text=True,
                   check=True, stdout=subprocess.PIPE)
