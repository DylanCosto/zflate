//! Adler-32 with SSE2, AVX2, AVX-VNNI and AVX-512 VNNI, ported from libdeflate's
//! `x86/adler32_template.h` and `x86/adler32_impl.h` (MIT).
//!
//! Over a chunk of `n` bytes `d_0 .. d_n-1`, Adler-32 adds `Σ d_i` to `s1` and
//! `n·s1 + Σ (n - i)·d_i` to `s2`. The vector code takes `VL` bytes a step. Per-lane sums of the
//! bytes give `s1`. For `s2`, each step's bytes are weighted `VL, VL - 1, .., 1`, and the `s1`
//! lanes so far are added into `s1_sums` before each step, which, multiplied by `VL` at the end,
//! weights each earlier byte by the `VL` positions every later step adds.
//!
//! The VNNI versions get both sums from `vpdpbusd`, which multiplies unsigned bytes by signed
//! bytes and adds the products in groups of four into 32-bit lanes: by the weights for `s2`,
//! and by ones for `s1`. The SSE2 and AVX2 versions sum bytes with `psadbw`, and keep each byte
//! position's total in a 16-bit lane, multiplied by its weight with `pmaddwd` once per chunk.

use std::arch::x86_64::*;

use super::{DIVISOR, Impl, MAX_CHUNK_LEN, scalar};

/// The implementations this CPU supports, fastest first, as libdeflate chooses them.
pub(super) fn available() -> Vec<(&'static str, Impl)> {
    let mut impls = Vec::new();
    if is_x86_feature_detected!("avx512bw")
        && is_x86_feature_detected!("avx512vnni")
        && wide_vectors_pay()
    {
        // SAFETY: returned only once the CPU features `vl512_vnni` needs are detected.
        let f: Impl = |adler, data| unsafe { vl512_vnni(adler, data) };
        impls.push(("avx512-vnni-512", f));
    }
    if is_x86_feature_detected!("avx512bw")
        && is_x86_feature_detected!("avx512vl")
        && is_x86_feature_detected!("avx512vnni")
    {
        // SAFETY: returned only once the CPU features `vl256_vnni` needs are detected.
        let f: Impl = |adler, data| unsafe { vl256_vnni(adler, data) };
        impls.push(("avx512-vnni-256", f));
    }
    if is_x86_feature_detected!("avx2") && is_x86_feature_detected!("avxvnni") {
        // SAFETY: returned only once the CPU features `avx2_vnni` needs are detected.
        let f: Impl = |adler, data| unsafe { avx2_vnni(adler, data) };
        impls.push(("avx-vnni", f));
    }
    if is_x86_feature_detected!("avx2") {
        // SAFETY: returned only once AVX2 is detected.
        let f: Impl = |adler, data| unsafe { avx2(adler, data) };
        impls.push(("avx2", f));
    }
    if is_x86_feature_detected!("sse2") {
        // SAFETY: returned only once SSE2 is detected.
        let f: Impl = |adler, data| unsafe { sse2(adler, data) };
        impls.push(("sse2", f));
    }
    impls
}

/// Whether 512-bit vectors pay: not on the Intel cores that lower their clock for a while after
/// using them, from Skylake-SP and Cascade Lake to Ice Lake and Tiger Lake, which would slow
/// the code around a short checksum (libdeflate's `allow_512bit_vectors`).
#[allow(
    unused_unsafe,
    reason = "`__cpuid` is unsafe in Rust 1.91, the oldest zflate builds with, and safe later"
)]
fn wide_vectors_pay() -> bool {
    // SAFETY: CPUID is on every x86-64 CPU.
    let (vendor, info) = unsafe { (__cpuid(0), __cpuid(1)) };
    let intel = [vendor.ebx, vendor.edx, vendor.ecx] == [0x756e_6547, 0x4965_6e69, 0x6c65_746e];
    let family = (info.eax >> 8) & 0xf;
    let model = (info.eax >> 4) & 0xf | (info.eax >> 12) & 0xf0;
    // Skylake-SP and its successors to Cooper Lake; Ice Lake server and client; Tiger Lake.
    !(intel && family == 6 && matches!(model, 85 | 106 | 108 | 126 | 140 | 141))
}

/// The weights of a step's bytes for the VNNI versions: `[VL, VL - 1, .., 1]`.
const fn descending<const VL: usize>() -> [u8; VL] {
    let mut w = [0; VL];
    let mut i = 0;
    while i < VL {
        w[i] = (VL - i) as u8;
        i += 1;
    }
    w
}

/// Aligns long inputs to `VL` bytes, adding the bytes before the boundary to the sums: an
/// unaligned load that crosses a cache line costs more. Short inputs take that cost instead.
#[inline(always)]
fn align<const VL: usize>(s1: &mut u32, s2: &mut u32, data: &[u8]) -> usize {
    let skip = data.as_ptr().align_offset(VL);
    if data.len() <= 65_536 || skip == 0 || skip >= VL {
        return 0;
    }
    for &b in &data[..skip] {
        *s1 += u32::from(b);
        *s2 += *s1;
    }
    *s1 %= DIVISOR;
    *s2 %= DIVISOR;
    skip
}

/// Sums the 32-bit lanes of a `s1` and a `s2` vector.
#[inline]
#[target_feature(enable = "sse2")]
fn reduce128(v_s1: __m128i, v_s2: __m128i) -> (u32, u32) {
    let v_s1 = _mm_add_epi32(v_s1, _mm_shuffle_epi32::<0x31>(v_s1));
    let v_s2 = _mm_add_epi32(v_s2, _mm_shuffle_epi32::<0x31>(v_s2));
    let v_s1 = _mm_add_epi32(v_s1, _mm_shuffle_epi32::<0x02>(v_s1));
    let v_s2 = _mm_add_epi32(v_s2, _mm_shuffle_epi32::<0x02>(v_s2));
    (
        _mm_cvtsi128_si32(v_s1) as u32,
        _mm_cvtsi128_si32(v_s2) as u32,
    )
}

#[inline]
#[target_feature(enable = "avx2")]
fn reduce256(v_s1: __m256i, v_s2: __m256i) -> (u32, u32) {
    reduce128(
        _mm_add_epi32(
            _mm256_castsi256_si128(v_s1),
            _mm256_extracti128_si256::<1>(v_s1),
        ),
        _mm_add_epi32(
            _mm256_castsi256_si128(v_s2),
            _mm256_extracti128_si256::<1>(v_s2),
        ),
    )
}

#[inline]
#[target_feature(enable = "avx512f")]
fn reduce512(v_s1: __m512i, v_s2: __m512i) -> (u32, u32) {
    reduce256(
        _mm256_add_epi32(
            _mm512_castsi512_si256(v_s1),
            _mm512_extracti64x4_epi64::<1>(v_s1),
        ),
        _mm256_add_epi32(
            _mm512_castsi512_si256(v_s2),
            _mm512_extracti64x4_epi64::<1>(v_s2),
        ),
    )
}

/// The last 1 to 64 bytes of a chunk, the rest of the vector zero.
#[inline]
#[target_feature(enable = "avx512bw")]
fn tail512(c: &[u8]) -> __m512i {
    debug_assert!((1..=64).contains(&c.len()));
    // SAFETY: the mask selects the first `c.len()` bytes, which `c` holds; a masked load does
    // not touch the others.
    unsafe { _mm512_maskz_loadu_epi8(u64::MAX >> (64 - c.len()), c.as_ptr().cast()) }
}

/// The last 1 to 32 bytes of a chunk, the rest of the vector zero.
#[inline]
#[target_feature(enable = "avx512bw,avx512vl")]
fn tail256_masked(c: &[u8]) -> __m256i {
    debug_assert!((1..=32).contains(&c.len()));
    // SAFETY: as in `tail512`.
    unsafe { _mm256_maskz_loadu_epi8(u32::MAX >> (32 - c.len()), c.as_ptr().cast()) }
}

/// The last 1 to 32 bytes of a chunk, the rest of the vector zero, without AVX-512's masks.
#[inline]
#[target_feature(enable = "avx2")]
fn tail256_copied(c: &[u8]) -> __m256i {
    let mut buf = [0; 32];
    buf[..c.len()].copy_from_slice(c);
    // SAFETY: `buf` is 32 bytes.
    unsafe { _mm256_loadu_si256(buf.as_ptr().cast()) }
}

/// Defines a VNNI Adler-32 over `VL`-byte vectors (libdeflate's template with `USE_VNNI`).
macro_rules! vnni {
    (
        $(#[$doc:meta])*
        fn $name:ident = $features:literal, $vl:literal, $vec:ty {
            add8: $add8:ident, add32: $add32:ident, dpbusd: $dpbusd:ident,
            loadu: $loadu:ident, mullo32: $mullo32:ident, set1_8: $set1_8:ident,
            set1_32: $set1_32:ident, zero: $zero:ident, sll32: $sll32:ident,
            tail: $tail:ident, reduce: $reduce:ident $(,)?
        }
    ) => {
        $(#[$doc])*
        #[target_feature(enable = $features)]
        fn $name(adler: u32, data: &[u8]) -> u32 {
            const VL: usize = $vl;
            const WEIGHTS: [u8; VL] = descending::<VL>();
            let load = |p: &[u8]| -> $vec {
                debug_assert!(p.len() >= VL);
                // SAFETY: `p` holds at least `VL` bytes.
                unsafe { $loadu(p.as_ptr().cast()) }
            };
            let shl = |v: $vec, by: usize| $sll32(v, _mm_cvtsi32_si128(by.trailing_zeros() as i32));
            let ones = $set1_8(1);
            let zero = $zero();
            let (mut s1, mut s2) = (adler & 0xffff, adler >> 16);
            let mut p = &data[align::<VL>(&mut s1, &mut s2, data)..];

            while !p.is_empty() {
                // Chunks short enough that the sums cannot overflow.
                let (chunk, rest) = p.split_at(p.len().min(MAX_CHUNK_LEN & !(4 * VL - 1)));
                p = rest;
                let mut weights = load(&WEIGHTS);
                let mut v_s1 = zero;
                let mut v_s2 = zero;
                // The chunk's bytes all come after the `s1` it starts with.
                s2 += s1 * chunk.len() as u32;

                let (steps, mut c) = chunk.as_chunks::<{ 4 * VL }>();
                if !steps.is_empty() {
                    // Four independent sets of sums, for the latency of `vpdpbusd`.
                    let [mut v_s1_b, mut v_s1_c, mut v_s1_d] = [zero; 3];
                    let [mut v_s2_b, mut v_s2_c, mut v_s2_d] = [zero; 3];
                    let [mut sums_a, mut sums_b, mut sums_c, mut sums_d] = [zero; 4];
                    for step in steps {
                        let (a, b, cc, d) = (
                            load(&step[..]),
                            load(&step[VL..]),
                            load(&step[2 * VL..]),
                            load(&step[3 * VL..]),
                        );
                        v_s2 = $dpbusd(v_s2, a, weights);
                        v_s2_b = $dpbusd(v_s2_b, b, weights);
                        v_s2_c = $dpbusd(v_s2_c, cc, weights);
                        v_s2_d = $dpbusd(v_s2_d, d, weights);
                        sums_a = $add32(sums_a, v_s1);
                        sums_b = $add32(sums_b, v_s1_b);
                        sums_c = $add32(sums_c, v_s1_c);
                        sums_d = $add32(sums_d, v_s1_d);
                        v_s1 = $dpbusd(v_s1, a, ones);
                        v_s1_b = $dpbusd(v_s1_b, b, ones);
                        v_s1_c = $dpbusd(v_s1_c, cc, ones);
                        v_s1_d = $dpbusd(v_s1_d, d, ones);
                    }
                    // Within a step, the second vector's bytes come `VL` later than the
                    // first's, and so on: s2 += 4VL·sums + 3VL·s1 + 2VL·s1_b + VL·s1_c.
                    let ab = $add32(v_s1, v_s1_b);
                    let ac = $add32(v_s1, v_s1_c);
                    let sums = $add32($add32(sums_a, sums_b), $add32(sums_c, sums_d));
                    v_s1 = $add32($add32(ab, v_s1_c), v_s1_d);
                    v_s2 = $add32(
                        $add32(
                            $add32(shl(sums, 4 * VL), shl(ab, 2 * VL)),
                            $add32(shl(ac, VL), v_s2),
                        ),
                        $add32($add32(v_s2_b, v_s2_c), v_s2_d),
                    );
                }

                // The last 0 to 4VL - 1 bytes of the chunk.
                if c.len() >= 2 * VL {
                    let (a, b) = (load(c), load(&c[VL..]));
                    v_s2 = $add32(v_s2, shl(v_s1, 2 * VL));
                    v_s1 = $dpbusd(v_s1, a, ones);
                    v_s1 = $dpbusd(v_s1, b, ones);
                    v_s2 = $dpbusd(v_s2, a, $set1_8(VL as i8));
                    v_s2 = $dpbusd(v_s2, a, weights);
                    v_s2 = $dpbusd(v_s2, b, weights);
                    c = &c[2 * VL..];
                }
                if !c.is_empty() {
                    // The `s1` lanes so far count once for each byte left, and the bytes
                    // left weigh `n, n - 1, ..`: zero or less for the padding.
                    let n = c.len();
                    v_s2 = $add32(v_s2, $mullo32(v_s1, $set1_32(n as i32)));
                    weights = $add8(weights, $set1_8((n as i32 - VL as i32) as i8));
                    if n > VL {
                        let data = load(c);
                        v_s1 = $dpbusd(v_s1, data, ones);
                        v_s2 = $dpbusd(v_s2, data, weights);
                        c = &c[VL..];
                        weights = $add8(weights, $set1_8(-(VL as i8)));
                    }
                    let data = $tail(c);
                    v_s1 = $dpbusd(v_s1, data, ones);
                    v_s2 = $dpbusd(v_s2, data, weights);
                }

                let (a, b) = $reduce(v_s1, v_s2);
                s1 = (s1 + a) % DIVISOR;
                s2 = (s2 + b) % DIVISOR;
            }
            (s2 << 16) | s1
        }
    };
}

vnni! {
    /// AVX-512 VNNI on 512-bit vectors, for CPUs whose clock 512-bit vectors do not lower.
    fn vl512_vnni = "avx512bw,avx512vnni", 64, __m512i {
        add8: _mm512_add_epi8, add32: _mm512_add_epi32, dpbusd: _mm512_dpbusd_epi32,
        loadu: _mm512_loadu_si512, mullo32: _mm512_mullo_epi32, set1_8: _mm512_set1_epi8,
        set1_32: _mm512_set1_epi32, zero: _mm512_setzero_si512, sll32: _mm512_sll_epi32,
        tail: tail512, reduce: reduce512,
    }
}

vnni! {
    /// AVX-512 VNNI on 256-bit vectors, using its masked loads: for Ice Lake and Tiger Lake,
    /// which lower their clock for 512-bit vectors.
    fn vl256_vnni = "avx512bw,avx512vl,avx512vnni", 32, __m256i {
        add8: _mm256_add_epi8, add32: _mm256_add_epi32, dpbusd: _mm256_dpbusd_epi32,
        loadu: _mm256_loadu_si256, mullo32: _mm256_mullo_epi32, set1_8: _mm256_set1_epi8,
        set1_32: _mm256_set1_epi32, zero: _mm256_setzero_si256, sll32: _mm256_sll_epi32,
        tail: tail256_masked, reduce: reduce256,
    }
}

vnni! {
    /// AVX-VNNI on 256-bit vectors, for CPUs without AVX-512 (Alder Lake and later).
    fn avx2_vnni = "avx2,avxvnni", 32, __m256i {
        add8: _mm256_add_epi8, add32: _mm256_add_epi32, dpbusd: _mm256_dpbusd_avx_epi32,
        loadu: _mm256_loadu_si256, mullo32: _mm256_mullo_epi32, set1_8: _mm256_set1_epi8,
        set1_32: _mm256_set1_epi32, zero: _mm256_setzero_si256, sll32: _mm256_sll_epi32,
        tail: tail256_copied, reduce: reduce256,
    }
}

/// Defines an SSE2 or AVX2 Adler-32 over `VL`-byte vectors (libdeflate's template without
/// `USE_VNNI`).
macro_rules! sad {
    (
        $(#[$doc:meta])*
        fn $name:ident = $features:literal, $vl:literal, $vec:ty, $weights:expr, {
            add16: $add16:ident, add32: $add32:ident, loadu: $loadu:ident, madd16: $madd16:ident,
            sad8: $sad8:ident, zero: $zero:ident, sll32: $sll32:ident,
            unpacklo8: $unpacklo8:ident, unpackhi8: $unpackhi8:ident, reduce: $reduce:ident $(,)?
        }
    ) => {
        $(#[$doc])*
        #[target_feature(enable = $features)]
        fn $name(adler: u32, data: &[u8]) -> u32 {
            const VL: usize = $vl;
            /// The weights `2VL .. 1` of a step's bytes, in the order `punpcklbw` and
            /// `punpckhbw` leave them, lane by 128-bit lane.
            const WEIGHTS: [[u16; VL / 2]; 4] = $weights;
            let load = |p: *const [u16; VL / 2]| -> $vec {
                // SAFETY: `p` points to `VL` bytes.
                unsafe { $loadu(p.cast()) }
            };
            let weights = WEIGHTS.each_ref().map(|w| load(w));
            let zero = $zero();
            let (mut s1, mut s2) = (adler & 0xffff, adler >> 16);
            let mut p = &data[align::<VL>(&mut s1, &mut s2, data)..];

            while !p.is_empty() {
                // Chunks short enough that neither the sums nor the 16-bit byte totals, which
                // `pmaddwd` takes as signed, can overflow.
                let max = (2 * VL * (i16::MAX as usize / 255)).min(MAX_CHUNK_LEN) & !(2 * VL - 1);
                let (chunk, rest) = p.split_at(p.len().min(max));
                p = rest;

                let (steps, c) = chunk.as_chunks::<{ 2 * VL }>();
                if !steps.is_empty() {
                    let mut v_s1 = zero;
                    let mut v_s1_sums = zero;
                    let mut totals = [zero; 4];
                    s2 += s1 * (2 * VL * steps.len()) as u32;
                    for step in steps {
                        // SAFETY: a step is `2VL` bytes.
                        let (a, b): ($vec, $vec) = unsafe {
                            ($loadu(step.as_ptr().cast()), $loadu(step[VL..].as_ptr().cast()))
                        };
                        v_s1_sums = $add32(v_s1_sums, v_s1);
                        totals[0] = $add16(totals[0], $unpacklo8(a, zero));
                        totals[1] = $add16(totals[1], $unpackhi8(a, zero));
                        totals[2] = $add16(totals[2], $unpacklo8(b, zero));
                        totals[3] = $add16(totals[3], $unpackhi8(b, zero));
                        v_s1 = $add32(v_s1, $add32($sad8(a, zero), $sad8(b, zero)));
                    }
                    // s2 += 2VL·s1_sums + Σ weight·total.
                    let shift = _mm_cvtsi32_si128((2 * VL).trailing_zeros() as i32);
                    let v_s2 = $add32(
                        $add32($sll32(v_s1_sums, shift), $madd16(totals[0], weights[0])),
                        $add32(
                            $add32($madd16(totals[1], weights[1]), $madd16(totals[2], weights[2])),
                            $madd16(totals[3], weights[3]),
                        ),
                    );
                    let (a, b) = $reduce(v_s1, v_s2);
                    s1 += a;
                    s2 += b;
                }
                // The last 0 to 2VL - 1 bytes, which also reduces the sums.
                scalar::chunk(&mut s1, &mut s2, c);
            }
            (s2 << 16) | s1
        }
    };
}

sad! {
    /// AVX2, for CPUs without VNNI.
    fn avx2 = "avx2", 32, __m256i, [
        [64, 63, 62, 61, 60, 59, 58, 57, 48, 47, 46, 45, 44, 43, 42, 41],
        [56, 55, 54, 53, 52, 51, 50, 49, 40, 39, 38, 37, 36, 35, 34, 33],
        [32, 31, 30, 29, 28, 27, 26, 25, 16, 15, 14, 13, 12, 11, 10, 9],
        [24, 23, 22, 21, 20, 19, 18, 17, 8, 7, 6, 5, 4, 3, 2, 1],
    ], {
        add16: _mm256_add_epi16, add32: _mm256_add_epi32, loadu: _mm256_loadu_si256,
        madd16: _mm256_madd_epi16, sad8: _mm256_sad_epu8, zero: _mm256_setzero_si256,
        sll32: _mm256_sll_epi32, unpacklo8: _mm256_unpacklo_epi8,
        unpackhi8: _mm256_unpackhi_epi8, reduce: reduce256,
    }
}

sad! {
    /// SSE2, which every x86-64 CPU has.
    fn sse2 = "sse2", 16, __m128i, [
        [32, 31, 30, 29, 28, 27, 26, 25],
        [24, 23, 22, 21, 20, 19, 18, 17],
        [16, 15, 14, 13, 12, 11, 10, 9],
        [8, 7, 6, 5, 4, 3, 2, 1],
    ], {
        add16: _mm_add_epi16, add32: _mm_add_epi32, loadu: _mm_loadu_si128,
        madd16: _mm_madd_epi16, sad8: _mm_sad_epu8, zero: _mm_setzero_si128,
        sll32: _mm_sll_epi32, unpacklo8: _mm_unpacklo_epi8, unpackhi8: _mm_unpackhi_epi8,
        reduce: reduce128,
    }
}
