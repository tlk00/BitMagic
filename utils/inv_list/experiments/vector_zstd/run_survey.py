#!/usr/bin/env python3
"""Find individually zstd-compressible L8 BLOBs in a paired distributed Gov2 sample."""
import argparse
import csv
import ctypes
import ctypes.util
from datetime import datetime
import hashlib
import json
from pathlib import Path
import time


class Zstd:
    """Own reusable native zstd buffers; each call produces an independent frame."""

    def __init__(self, level):
        """Load libzstd and retain the requested compression level."""
        library = ctypes.util.find_library("zstd") or "/opt/homebrew/lib/libzstd.dylib"
        self.lib = ctypes.CDLL(library)
        for name, args, result in (
            ("ZSTD_compressBound", [ctypes.c_size_t], ctypes.c_size_t),
            ("ZSTD_compress", [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p,
                               ctypes.c_size_t, ctypes.c_int], ctypes.c_size_t),
            ("ZSTD_decompress", [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p,
                                 ctypes.c_size_t], ctypes.c_size_t),
            ("ZSTD_isError", [ctypes.c_size_t], ctypes.c_uint),
            ("ZSTD_getErrorName", [ctypes.c_size_t], ctypes.c_char_p),
            ("ZSTD_versionString", [], ctypes.c_char_p),
        ):
            function = getattr(self.lib, name)
            function.argtypes = args
            function.restype = result
        self.level = level

    def checked(self, value):
        """Return a native result or raise with the library's diagnostic."""
        if self.lib.ZSTD_isError(value):
            raise RuntimeError(self.lib.ZSTD_getErrorName(value).decode())
        return value

    def compress(self, source):
        """Return an independent frame and verify exact decompression equality."""
        target = ctypes.create_string_buffer(self.lib.ZSTD_compressBound(len(source)))
        begin = time.perf_counter()
        size = self.checked(self.lib.ZSTD_compress(target, len(target), source,
                                                  len(source), self.level))
        encode_s = time.perf_counter() - begin
        frame = target.raw[:size]
        restored = ctypes.create_string_buffer(len(source))
        begin = time.perf_counter()
        restored_size = self.checked(self.lib.ZSTD_decompress(restored, len(source), frame, size))
        decode_s = time.perf_counter() - begin
        if restored_size != len(source) or restored.raw != source:
            raise ValueError("zstd round-trip mismatch")
        return frame, encode_s, decode_s


def index(path):
    """Validate the frozen collection index and return its vector boundaries."""
    with path.open() as stream:
        signature, size, count = stream.readline().split()
        positions = list(map(int, stream))
    if signature != "BM_INV_BLOB_INDEX_V1" or len(positions) != int(count)+1:
        raise ValueError("Invalid reference index")
    if positions[0] or positions[-1] != int(size):
        raise ValueError("Invalid index endpoints")
    if any(b <= a for a, b in zip(positions, positions[1:])):
        raise ValueError("Invalid index order")
    return positions


def dataset(args, codec, name, ordinals):
    """Survey one reference and retain both absolute and percentage winners."""
    prefix = args.reference / f"gov2.{name}.level-{args.bm_level:02}"
    offsets = index(Path(str(prefix)+".idx"))
    source = Path(str(prefix)+".bv")
    if source.stat().st_size != offsets[-1] or ordinals[-1] >= len(offsets)-1:
        raise ValueError("Reference size or ordinal mismatch")
    out = args.out / name
    out.mkdir()
    rows = []
    total_encode = total_decode = 0.
    with source.open("rb") as stream:
        for i, ordinal in enumerate(ordinals):
            size = offsets[ordinal+1]-offsets[ordinal]
            stream.seek(offsets[ordinal])
            blob = stream.read(size)
            if len(blob) != size:
                raise ValueError("Truncated reference")
            frame, encode_s, decode_s = codec.compress(blob)
            saved = size-len(frame)
            rows.append(dict(ordinal=ordinal, bm_bytes=size, zstd_bytes=len(frame),
                             saved_bytes=saved, saved_pct=100*saved/size,
                             encode_s=encode_s, decode_s=decode_s,
                             sha256=hashlib.sha256(blob).hexdigest()))
            total_encode += encode_s
            total_decode += decode_s
            if (i+1) % 2000 == 0:
                print(f"{name}: completed {i+1}/{len(ordinals)}", flush=True)
    with (out/"vectors.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    winning = [r for r in rows if r["saved_bytes"] > 0]
    absolute = sorted(winning, key=lambda r: (-r["saved_bytes"], r["ordinal"]))[:args.keep]
    percentage = sorted((r for r in winning if r["bm_bytes"] >= args.ratio_min_bytes),
                        key=lambda r: (-r["saved_pct"], r["ordinal"]))[:args.keep]
    chosen = {r["ordinal"]: r for r in absolute+percentage}
    (out/"absolute.ordinals.txt").write_text("".join(f'{r["ordinal"]}\n' for r in absolute))
    (out/"percentage.ordinals.txt").write_text("".join(f'{r["ordinal"]}\n' for r in percentage))
    (out/"selected.ordinals.txt").write_text("".join(f"{i}\n" for i in sorted(chosen)))
    examples = out/"examples"
    examples.mkdir()
    with source.open("rb") as stream:
        for ordinal in sorted(chosen):
            row = chosen[ordinal]
            stream.seek(offsets[ordinal])
            blob = stream.read(row["bm_bytes"])
            frame, _, _ = codec.compress(blob)
            if len(frame) != row["zstd_bytes"]:
                raise ValueError("Non-reproducible frame size")
            stem = examples/f"gov2.{name}.ordinal-{ordinal:07}.level-{args.bm_level:02}.bv"
            stem.write_bytes(blob)
            Path(str(stem)+".zst").write_bytes(frame)
    baseline = sum(r["bm_bytes"] for r in rows)
    saved = sum(r["saved_bytes"] for r in winning)
    summary = dict(dataset=name, vectors=len(rows), baseline_bytes=baseline,
                   all_frames_bytes=sum(r["zstd_bytes"] for r in rows),
                   winning_vectors=len(winning), selected_saved_bytes=saved,
                   selective_saved_pct=100*saved/baseline, encode_s=total_encode,
                   decode_s=total_decode, retained_vectors=len(chosen),
                   top_absolute_saved_bytes=sum(r["saved_bytes"] for r in absolute))
    return summary


def main():
    """Run both datasets with saved ordinals and write provenance and weighted totals."""
    usb = Path("/Volumes/USB_1TB")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=usb/"gov2-reference-Nov-2026.v3")
    parser.add_argument("--ordinals", type=Path, default=usb/"BitMagic-experiments/gov2-op68-split-20261113-094058/ordinals.txt")
    parser.add_argument("--out", type=Path, default=usb/"BitMagic-experiments"/
                        ("gov2-vector-zstd-"+datetime.now().strftime("%Y%m%d-%H%M%S")))
    parser.add_argument("--bm-level", type=int, choices=range(1, 11), default=8)
    parser.add_argument("--zstd-level", type=int, default=3)
    parser.add_argument("--keep", type=int, default=100)
    parser.add_argument("--ratio-min-bytes", type=int, default=1024)
    args = parser.parse_args()
    ordinals = list(map(int, args.ordinals.read_text().split()))
    if not ordinals or ordinals != sorted(set(ordinals)) or args.keep <= 0:
        parser.error("Invalid sample or retention count")
    if usb in args.out.resolve().parents and not usb.is_mount():
        parser.error("USB_1TB unavailable")
    args.out.mkdir(parents=True, exist_ok=False)
    codec = Zstd(args.zstd_level)
    provenance = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    provenance["zstd_version"] = codec.lib.ZSTD_versionString().decode()
    (args.out/"settings.json").write_text(json.dumps(provenance, indent=2)+"\n")
    (args.out/"ordinals.txt").write_text(args.ordinals.read_text())
    print(f"Reports: {args.out}", flush=True)
    summaries = [dataset(args, codec, name, ordinals) for name in ("sorted", "unsorted")]
    with (args.out/"report.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summaries[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summaries)
    for row in summaries:
        print(row, flush=True)


if __name__ == "__main__":
    main()
