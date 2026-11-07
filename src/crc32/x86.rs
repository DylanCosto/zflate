//! CRC32 by folding with carry-less multiplication: PCLMULQDQ on 128-bit vectors, and
//! VPCLMULQDQ with AVX-512 on 512-bit ones. The algorithm and its structure follow
//! libdeflate's `x86/crc32_pclmul_template.h` (MIT); its `scripts/gen-crc32-consts.py`
//! derives the constants.
//!
//! The message is a polynomial over GF(2). Multiplying a 128-bit piece by `x^d mod G` moves
//! it `d` bits later without changing the remainder, so it can be XORed into the piece
//! there ("folding"). Eight accumulators fold in parallel to hide the multiplier's latency,
//! then fold into each other, and a Barrett reduction turns the last 128 bits into the CRC.

use std::arch::x86_64::*;

use super::{Impl, scalar, x_pow_mod_g};

/// The multipliers that fold a 128-bit lane `d` bits forward, for its low and high 64 bits.
///
/// Only for the `FOLD_*` constants: called at run time, it would loop `d` times per call.
const fn mults(d: u32) -> (u32, u32) {
    (x_pow_mod_g(d + 31), x_pow_mod_g(d - 33))
}

const FOLD_128: (u32, u32) = mults(128);
const FOLD_256: (u32, u32) = mults(256);
const FOLD_512: (u32, u32) = mults(512);
const FOLD_1024: (u32, u32) = mults(1024);
const FOLD_2048: (u32, u32) = mults(2048);
const FOLD_4096: (u32, u32) = mults(4096);

/// `floor(x^95 / G(x))` and `G(x)`, for the final Barrett reduction (libdeflate's
/// `CRC32_BARRETT_CONSTANT_1` and `_2`).
const BARRETT: (u64, u64) = (0xb4e5_b025_f701_1641, 0x0000_0001_db71_0641);

/// `pshufb` with `SHIFT_TAB[n..n + 16]` shifts left by `16 - n` bytes, and with
/// `SHIFT_TAB[n + 16..n + 32]` shifts right by `n` bytes, zero-filling (`0xff` selects zero).
static SHIFT_TAB: [u8; 48] = {
    let mut t = [0xff; 48];
    let mut i = 0;
    while i < 16 {
        t[16 + i] = i as u8;
        i += 1;
    }
    t
};

/// The implementations this CPU supports, fastest first.
pub(super) fn available() -> Vec<(&'static str, Impl)> {
    let mut impls = Vec::new();
    // Every feature `vl512` is built with, though AVX-512 CPUs have SSE4.1 and AVX2: a
    // virtual machine may mask features oddly.
    if is_x86_feature_detected!("vpclmulqdq")
        && is_x86_feature_detected!("pclmulqdq")
        && is_x86_feature_detected!("sse4.1")
        && is_x86_feature_detected!("avx2")
        && is_x86_feature_detected!("avx512f")
        && is_x86_feature_detected!("avx512bw")
        && is_x86_feature_detected!("avx512vl")
    {
        // SAFETY: returned only once the CPU features `vl512` needs are detected.
        let f: Impl = |crc, data| unsafe { vl512(crc, data) };
        impls.push(("vpclmulqdq-avx512", f));
    }
    if is_x86_feature_detected!("pclmulqdq") && is_x86_feature_detected!("sse4.1") {
        // SAFETY: returned only once the CPU features `vl16` needs are detected.
        let f: Impl = |crc, data| unsafe { vl16(crc, data) };
        impls.push(("pclmulqdq", f));
    }
    impls
}

#[inline]
#[target_feature(enable = "sse2")]
fn mults128((lo, hi): (u32, u32)) -> __m128i {
    _mm_set_epi64x(i64::from(hi), i64::from(lo))
}

#[inline]
#[target_feature(enable = "sse2")]
fn load128(b: &[u8; 16]) -> __m128i {
    // SAFETY: `b` is 16 readable bytes, and the load is unaligned.
    unsafe { _mm_loadu_si128(b.as_ptr().cast()) }
}

/// Loads `N` consecutive 128-bit lanes from the start of `b`, which holds at least `16 * N`
/// bytes.
#[inline]
#[target_feature(enable = "sse2")]
fn load128_lanes<const N: usize>(b: &[u8]) -> [__m128i; N] {
    let mut v = [_mm_setzero_si128(); N];
    for (lane, chunk) in v.iter_mut().zip(b.as_chunks::<16>().0) {
        *lane = load128(chunk);
    }
    v
}

/// Folds `src` forward onto `dst`.
#[inline]
#[target_feature(enable = "pclmulqdq,sse2")]
fn fold128(src: __m128i, dst: __m128i, mults: __m128i) -> __m128i {
    let lo = _mm_clmulepi64_si128::<0x00>(src, mults);
    let hi = _mm_clmulepi64_si128::<0x11>(src, mults);
    _mm_xor_si128(_mm_xor_si128(dst, lo), hi)
}

/// Folds the 128-bit `x` over the last `1 <= n <= 15` bytes of `data`: it splits `x` and the
/// bytes into an `n`-byte piece and a 16-byte one, and folds the first onto the second.
#[inline]
#[target_feature(enable = "pclmulqdq,sse4.1")]
fn fold_tail(x: __m128i, data: &[u8], n: usize) -> __m128i {
    let last = data.last_chunk::<16>().expect("folding needs 16 bytes");
    let lshift = load128(SHIFT_TAB[n..].first_chunk().expect("n < 16"));
    let rshift = load128(SHIFT_TAB[n + 16..].first_chunk().expect("n < 16"));
    let x0 = _mm_shuffle_epi8(x, lshift);
    // The high bit of each `rshift` byte picks the data byte over the shifted `x` byte.
    let x1 = _mm_blendv_epi8(_mm_shuffle_epi8(x, rshift), load128(last), rshift);
    fold128(x0, x1, mults128(FOLD_128))
}

/// Reduces the 128-bit `x`, multiplied by `x^32`, modulo `G`: the (inverted) CRC.
#[inline]
#[target_feature(enable = "pclmulqdq,sse4.1")]
fn reduce(x: __m128i) -> u32 {
    let barrett = _mm_set_epi64x(BARRETT.1.cast_signed(), BARRETT.0.cast_signed());
    let x = _mm_xor_si128(
        _mm_clmulepi64_si128::<0x10>(x, mults128(FOLD_128)),
        _mm_bsrli_si128::<8>(x),
    );
    let q = _mm_clmulepi64_si128::<0x00>(x, barrett);
    let q = _mm_clmulepi64_si128::<0x10>(q, barrett);
    _mm_extract_epi32::<2>(_mm_xor_si128(x, q)).cast_unsigned()
}

/// Folds what remains of `data` after the 128-bit accumulator `x`, in 16-byte lanes then a
/// shorter tail, and reduces. `data` is the whole input, at least 16 bytes, and `rest` its
/// suffix, of fewer than `MAX_LANES + 1` lanes; the bound lets the lane loop compile to a
/// few tests instead of a counted loop.
#[inline]
#[target_feature(enable = "pclmulqdq,sse4.1")]
fn finish<const MAX_LANES: usize>(mut x: __m128i, data: &[u8], mut rest: &[u8]) -> u32 {
    debug_assert!(
        rest.len() < 16 * (MAX_LANES + 1),
        "{} bytes left",
        rest.len()
    );
    let m = mults128(FOLD_128);
    for _ in 0..MAX_LANES {
        let Some((lane, r)) = rest.split_first_chunk::<16>() else {
            break;
        };
        x = fold128(x, load128(lane), m);
        rest = r;
    }
    let tail = rest;
    if !tail.is_empty() {
        x = fold_tail(x, data, tail.len());
    }
    reduce(x)
}

/// PCLMULQDQ on 128-bit vectors: eight lanes, then four, two, one.
#[target_feature(enable = "pclmulqdq,sse4.1")]
fn vl16(crc: u32, data: &[u8]) -> u32 {
    if data.len() < 16 {
        return scalar::slice1(crc, data);
    }
    let crc = _mm_cvtsi32_si128(crc.cast_signed());
    let (m8, m4, m2, m1) = (
        mults128(FOLD_1024),
        mults128(FOLD_512),
        mults128(FOLD_256),
        mults128(FOLD_128),
    );
    let lanes = match data.len() {
        128.. => 8,
        64.. => 4,
        32.. => 2,
        _ => 1,
    };
    let mut rest = &data[16 * lanes..];

    let mut v4: [__m128i; 4];
    let mut v2: [__m128i; 2];
    let x: __m128i;
    if lanes == 8 {
        let mut v8: [__m128i; 8] = load128_lanes(data);
        v8[0] = _mm_xor_si128(v8[0], crc);
        // Fold eight lanes while eight more remain; then the eight onto four.
        while rest.len() >= 128 {
            let (group, r) = rest.split_at(128);
            let next: [__m128i; 8] = load128_lanes(group);
            for (acc, n) in v8.iter_mut().zip(next) {
                *acc = fold128(*acc, n, m8);
            }
            rest = r;
        }
        v4 = [v8[0], v8[1], v8[2], v8[3]];
        for (acc, n) in v4.iter_mut().zip(&v8[4..]) {
            *acc = fold128(*acc, *n, m4);
        }
    } else if lanes == 4 {
        v4 = load128_lanes(data);
        v4[0] = _mm_xor_si128(v4[0], crc);
    } else {
        v4 = [_mm_setzero_si128(); 4];
    }

    if lanes >= 4 {
        rest = fold_group(&mut v4, rest, m4);
        v2 = [fold128(v4[0], v4[2], m2), fold128(v4[1], v4[3], m2)];
    } else if lanes == 2 {
        v2 = load128_lanes(data);
        v2[0] = _mm_xor_si128(v2[0], crc);
    } else {
        v2 = [_mm_setzero_si128(); 2];
    }

    if lanes >= 2 {
        rest = fold_group(&mut v2, rest, m2);
        x = fold128(v2[0], v2[1], m1);
    } else {
        x = _mm_xor_si128(load128(data.first_chunk().expect("len >= 16")), crc);
    }
    finish::<1>(x, data, rest)
}

/// Folds one group of `N` lanes from `rest` into `v` if `rest` has that many.
#[inline]
#[target_feature(enable = "pclmulqdq,sse2")]
fn fold_group<'d, const N: usize>(
    v: &mut [__m128i; N],
    rest: &'d [u8],
    mults: __m128i,
) -> &'d [u8] {
    if rest.len() < 16 * N {
        return rest;
    }
    let (group, r) = rest.split_at(16 * N);
    let next: [__m128i; N] = load128_lanes(group);
    for (acc, n) in v.iter_mut().zip(next) {
        *acc = fold128(*acc, n, mults);
    }
    r
}

#[inline]
#[target_feature(enable = "avx512f")]
fn mults512(m: (u32, u32)) -> __m512i {
    _mm512_broadcast_i32x4(mults128(m))
}

#[inline]
#[target_feature(enable = "avx512f")]
fn load512(b: &[u8; 64]) -> __m512i {
    // SAFETY: `b` is 64 readable bytes, and the load is unaligned.
    unsafe { _mm512_loadu_si512(b.as_ptr().cast()) }
}

/// Loads a group of `N` 512-bit vectors.
#[inline]
#[target_feature(enable = "avx512f")]
fn load512_group<const N: usize>(group: &[[u8; 64]; N]) -> [__m512i; N] {
    let mut v = [_mm512_setzero_si512(); N];
    for (lane, chunk) in v.iter_mut().zip(group) {
        *lane = load512(chunk);
    }
    v
}

/// Folds the accumulators `v` forward onto the next group of `N` vectors.
#[inline]
#[target_feature(enable = "avx512f,vpclmulqdq")]
fn fold512_onto<const N: usize>(v: &mut [__m512i; N], group: &[[u8; 64]; N], mults: __m512i) {
    for (acc, chunk) in v.iter_mut().zip(group) {
        *acc = fold512(*acc, load512(chunk), mults);
    }
}

/// Folds four 128-bit lanes of `src` forward onto `dst` at once; `vpternlogd` with `0x96`
/// is a three-way XOR.
#[inline]
#[target_feature(enable = "avx512f,vpclmulqdq")]
fn fold512(src: __m512i, dst: __m512i, mults: __m512i) -> __m512i {
    _mm512_ternarylogic_epi32::<0x96>(
        _mm512_clmulepi64_epi128::<0x00>(src, mults),
        _mm512_clmulepi64_epi128::<0x11>(src, mults),
        dst,
    )
}

#[inline]
#[target_feature(enable = "avx512f,avx512vl,vpclmulqdq")]
fn fold256(src: __m256i, dst: __m256i, mults: __m256i) -> __m256i {
    _mm256_ternarylogic_epi32::<0x96>(
        _mm256_clmulepi64_epi128::<0x00>(src, mults),
        _mm256_clmulepi64_epi128::<0x11>(src, mults),
        dst,
    )
}

/// Loads the first group of `N` vectors, with the CRC XORed into the first: the starting
/// accumulators.
#[inline]
#[target_feature(enable = "avx512f")]
fn start512<const N: usize>(group: &[[u8; 64]; N], crc: __m512i) -> [__m512i; N] {
    let mut v = load512_group(group);
    v[0] = _mm512_xor_si512(v[0], crc);
    v
}

/// One stage of the cascade: the `N` accumulators fold onto the next group of `N` vectors if
/// there is one, then their first half folds onto their second. Returns the `H = N / 2`
/// accumulators and the vectors left.
///
/// Each entry point in [`vl512`] calls the stages in a straight line, as libdeflate's
/// `goto`s do: merging accumulators that came from different stages sends them through the
/// stack.
#[inline]
#[target_feature(enable = "avx512f,vpclmulqdq")]
fn stage512<const N: usize, const H: usize>(
    mut v: [__m512i; N],
    mut vecs: &[[u8; 64]],
    mults_n: __m512i,
    mults_h: __m512i,
) -> ([__m512i; H], &[[u8; 64]]) {
    const { assert!(N == 2 * H, "a stage halves the accumulators") };
    if let Some((group, rest)) = vecs.split_first_chunk::<N>() {
        fold512_onto(&mut v, group, mults_n);
        vecs = rest;
    }
    let mut half = [_mm512_setzero_si512(); H];
    for (i, acc) in half.iter_mut().enumerate() {
        *acc = fold512(v[i], v[i + H], mults_h);
    }
    (half, vecs)
}

/// The end of [`vl512`]: folds the last accumulator onto a last vector if one remains, then
/// from 512 bits to 256, folding in a 32-byte lane if one remains, to 128, and finishes.
#[inline]
#[target_feature(enable = "pclmulqdq,sse4.1,avx2,avx512f,avx512vl,vpclmulqdq")]
fn end512(mut z: __m512i, vecs: &[[u8; 64]], mut rest: &[u8], data: &[u8]) -> u32 {
    if let [last] = vecs {
        z = fold512(z, load512(last), mults512(FOLD_512));
    }
    let m256 = _mm256_broadcastsi128_si256(mults128(FOLD_256));
    let mut y = fold256(
        _mm512_extracti64x4_epi64::<0>(z),
        _mm512_extracti64x4_epi64::<1>(z),
        m256,
    );
    if let Some((lane, r)) = rest.split_first_chunk::<32>() {
        // SAFETY: `lane` is 32 readable bytes, and the load is unaligned.
        let next = unsafe { _mm256_loadu_si256(lane.as_ptr().cast()) };
        y = fold256(y, next, m256);
        rest = r;
    }
    let x = fold128(
        _mm256_extracti128_si256::<0>(y),
        _mm256_extracti128_si256::<1>(y),
        mults128(FOLD_128),
    );
    finish::<1>(x, data, rest)
}

/// VPCLMULQDQ with AVX-512 on 512-bit vectors (four 128-bit lanes each): eight vectors,
/// then four, two, one, then down to 128 bits. Inputs under 64 bytes stay in 128-bit
/// vectors, so they are no slower than with [`vl16`].
#[target_feature(enable = "pclmulqdq,sse4.1,avx2,avx512f,avx512bw,avx512vl,vpclmulqdq")]
fn vl512(crc: u32, data: &[u8]) -> u32 {
    let len = data.len();
    let crc128 = _mm_cvtsi32_si128(crc.cast_signed());
    if len < 16 {
        if len < 4 {
            return scalar::slice1(crc, data);
        }
        // Load just the input bytes, XOR in the CRC, and shift the result to the top of
        // the register, where the reduction expects a 128-bit polynomial's low terms.
        let mask = (1u16 << len) - 1;
        // SAFETY: the mask selects only the `len` bytes of `data`; AVX-512 does not access
        // the bytes it masks off.
        let x = unsafe { _mm_maskz_loadu_epi8(mask, data.as_ptr().cast()) };
        let shift = load128(SHIFT_TAB[len..].first_chunk().expect("len < 16"));
        return reduce(_mm_shuffle_epi8(_mm_xor_si128(x, crc128), shift));
    }
    if len < 64 {
        let (first, rest) = data.split_first_chunk::<16>().expect("len >= 16");
        return finish::<3>(_mm_xor_si128(load128(first), crc128), data, rest);
    }

    let crc512 = _mm512_zextsi128_si512(crc128);
    let (m8, m4, m2, m1) = (
        mults512(FOLD_4096),
        mults512(FOLD_2048),
        mults512(FOLD_1024),
        mults512(FOLD_512),
    );
    let (vecs, rest) = data.as_chunks::<64>();
    // Start with as many accumulators as there are vectors for, up to eight.
    if let Some((group, vecs)) = vecs.split_first_chunk::<8>() {
        let mut v8 = start512(group, crc512);
        let (groups, vecs) = vecs.as_chunks::<8>();
        for group in groups {
            fold512_onto(&mut v8, group, m8);
        }
        let (v4, vecs) = stage512::<8, 4>(v8, vecs, m8, m4);
        let (v2, vecs) = stage512::<4, 2>(v4, vecs, m4, m2);
        let ([z], vecs) = stage512::<2, 1>(v2, vecs, m2, m1);
        end512(z, vecs, rest, data)
    } else if let Some((group, vecs)) = vecs.split_first_chunk::<4>() {
        let (v2, vecs) = stage512::<4, 2>(start512(group, crc512), vecs, m4, m2);
        let ([z], vecs) = stage512::<2, 1>(v2, vecs, m2, m1);
        end512(z, vecs, rest, data)
    } else if let Some((group, vecs)) = vecs.split_first_chunk::<2>() {
        let ([z], vecs) = stage512::<2, 1>(start512(group, crc512), vecs, m2, m1);
        end512(z, vecs, rest, data)
    } else {
        let (group, vecs) = vecs.split_first_chunk::<1>().expect("len >= 64");
        let [z] = start512(group, crc512);
        end512(z, vecs, rest, data)
    }
}
