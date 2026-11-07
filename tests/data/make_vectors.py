"""Writes vectors.bin, the small fixed inputs of tests/vectors.rs (which also runs under
Miri, where libdeflate's C code cannot): streams compressed by zlib, each with its data.

usage: python3 -I make_vectors.py [real BGZF file]...

Record: format (0 raw DEFLATE, 1 zlib, 2 gzip), compressed length, data length (u32 LE
each), compressed bytes, data bytes. Each real file contributes a block from its middle, as
raw DEFLATE. The committed file used public data: GIAB HG002 PacBio HiFi (chr20) and the
1000 Genomes high-coverage chr22 VCF.
"""
import os
import random
import struct
import sys
import zlib

rng = random.Random(7)


def inputs():
    for n in (0, 1, 2, 100, 1000, 4096):
        yield bytes(rng.choice(b"ACGT") for _ in range(n))  # pairs: compresses about 4:1
        yield bytes(rng.randrange(256) for _ in range(n))  # stored or nearly
        yield bytes(i % 7 for i in range(n))  # runs and short offsets
        words = [b"chr1\t", b"10583\t", b"rs58108140\t", b"G\tA\t", b"PASS\n"]
        text = b"".join(rng.choice(words) for _ in range(n))[:n]
        yield text  # compresses well: no pairs


def records():
    for data in inputs():
        for level, strategy in ((0, zlib.Z_DEFAULT_STRATEGY), (1, zlib.Z_DEFAULT_STRATEGY),
                                (6, zlib.Z_DEFAULT_STRATEGY), (9, zlib.Z_FIXED),
                                (9, zlib.Z_HUFFMAN_ONLY), (9, zlib.Z_RLE)):
            for fmt, wbits in ((0, -15), (1, 15), (2, 31)):
                c = zlib.compressobj(level, zlib.DEFLATED, wbits, 9, strategy)
                yield fmt, c.compress(data) + c.flush(), data
    for path in sys.argv[1:]:
        # One whole block from the middle of each file.
        with open(path, "rb") as f:
            f.seek(0, 2)
            f.seek(f.tell() // 2)
            window = f.read(1 << 18)
        pos = 0
        while True:
            pos = window.index(b"\x1f\x8b\x08\x04", pos)
            bsize = struct.unpack_from("<H", window, pos + 16)[0] + 1
            if window[pos + bsize:pos + bsize + 4] == b"\x1f\x8b\x08\x04":
                break
            pos += 1
        xlen = struct.unpack_from("<H", window, pos + 10)[0]
        deflate = window[pos + 12 + xlen:pos + bsize - 8]
        yield 0, deflate, zlib.decompress(deflate, -15)


def main():
    out = bytearray()
    n = 0
    for fmt, z, data in records():
        out += struct.pack("<BII", fmt, len(z), len(data)) + z + data
        n += 1
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "vectors.bin")
    with open(path, "wb") as f:
        f.write(out)
    print(f"{n} vectors, {len(out)} bytes, to {path}")


main()
