#!/usr/bin/env bash
# Checks zflate on a 64-bit Arm machine (Apple silicon, Graviton, Ampere): the test suite,
# which compares every CRC32 implementation and both codecs with libdeflate, then speed
# against libdeflate. Writes everything to arm-check-<time>.txt to send back.
#
# usage: benches/arm-check.sh [file ...]   (BGZF files: BAM, BCF or bgzipped VCF)
# needs: Rust (rustup) and a C compiler for libdeflate (macOS: xcode-select --install)
set -euo pipefail
cd "$(dirname "$0")/.."
out="arm-check-$(date +%Y%m%d-%H%M%S).txt"
{
    echo "== machine"
    uname -a
    sysctl -n machdep.cpu.brand_string 2>/dev/null || grep -m1 -i "model name" /proc/cpuinfo || true
    rustc --version

    echo "== tests"
    cargo test --release 2>&1 | grep -E "test result|panicked|FAILED|^error" || true
    cargo test --release --test differential every_crc32 -- --nocapture 2>&1 \
        | grep "crc32 implementations" || true

    echo "== crc32 speed"
    cargo bench --bench crc32 -- --warm-up-time 1 --measurement-time 2 2>&1 \
        | grep -E "^crc32/|thrpt:" || true

    for file in "$@"; do
        echo "== $file: compression, levels 1, 5 and 6"
        ZFLATE_WINDOW_MB=4 cargo run -q --release --example compress -- compare "$file" 1,5,6 3
        echo "== $file: decompression"
        cargo run -q --release --example validate -- bgzf-time "$file"
    done
} 2>&1 | tee "$out"
echo "Results in $out"
