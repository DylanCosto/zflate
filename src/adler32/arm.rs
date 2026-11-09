//! Adler-32 with NEON, and with the dot-product instructions, ported from libdeflate's
//! `arm/adler32_impl.h` (MIT). The sums are as in [`super::x86`], 64 bytes a step.
//!
//! The NEON version adds the bytes pairwise into `s1`, and keeps each of the 64 byte positions'
//! totals in a 16-bit lane, multiplied by its weight once per chunk. The dot-product version
//! gets both sums from `udot`, which multiplies bytes and adds them in groups of four into
//! 32-bit lanes: by ones for `s1`, and by the weights `64 .. 1` for `s2`.

use std::arch::aarch64::*;

use super::{DIVISOR, Impl, MAX_CHUNK_LEN, scalar};

/// The implementations this CPU supports, fastest first.
pub(super) fn available() -> Vec<(&'static str, Impl)> {
    let mut impls = Vec::new();
    if std::arch::is_aarch64_feature_detected!("neon")
        && std::arch::is_aarch64_feature_detected!("dotprod")
    {
        // SAFETY: returned only once NEON and the dot-product instructions are detected.
        let f: Impl = |adler, data| unsafe { neon_dotprod(adler, data) };
        impls.push(("neon-dotprod", f));
    }
    if std::arch::is_aarch64_feature_detected!("neon") {
        // SAFETY: returned only once NEON is detected.
        let f: Impl = |adler, data| unsafe { neon(adler, data) };
        impls.push(("neon", f));
    }
    impls
}

/// The weights of a step's 64 bytes.
const WEIGHTS: [u8; 64] = {
    let mut w = [0; 64];
    let mut i = 0;
    while i < 64 {
        w[i] = 64 - i as u8;
        i += 1;
    }
    w
};

/// The largest chunk: a multiple of a step's 64 bytes.
const CHUNK: usize = MAX_CHUNK_LEN & !63;

/// Aligns long inputs to 16 bytes, adding the bytes before the boundary to the sums. Short
/// inputs take the cost of unaligned loads instead.
#[inline(always)]
fn align(s1: &mut u32, s2: &mut u32, data: &[u8]) -> usize {
    let skip = data.as_ptr().align_offset(16);
    if data.len() <= 32_768 || skip == 0 || skip >= 16 {
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

/// Loads four 16-byte vectors.
#[inline]
#[target_feature(enable = "neon")]
fn load4(p: &[u8; 64]) -> [uint8x16_t; 4] {
    // SAFETY: each load reads 16 of the 64 bytes `p` holds.
    unsafe { [0, 16, 32, 48].map(|i| vld1q_u8(p[i..].as_ptr())) }
}

/// NEON.
#[target_feature(enable = "neon")]
fn neon(adler: u32, data: &[u8]) -> u32 {
    let w = load4(&WEIGHTS).map(|v| [vmovl_u8(vget_low_u8(v)), vmovl_high_u8(v)]);
    let (mut s1, mut s2) = (adler & 0xffff, adler >> 16);
    let mut p = &data[align(&mut s1, &mut s2, data)..];

    while !p.is_empty() {
        let (chunk, rest) = p.split_at(p.len().min(CHUNK));
        p = rest;

        let (steps, c) = chunk.as_chunks::<64>();
        if !steps.is_empty() {
            let mut v_s1 = vdupq_n_u32(0);
            let mut v_s2 = vdupq_n_u32(0);
            // The total of each of the 64 byte positions across the steps.
            let mut totals = [vdupq_n_u16(0); 8];
            s2 += s1 * (64 * steps.len()) as u32;

            for step in steps {
                let d = load4(step);
                // The `s1` lanes so far, multiplied by 64 at the end.
                v_s2 = vaddq_u32(v_s2, v_s1);
                let mut pairs = vpaddlq_u8(d[0]);
                for (i, v) in d.iter().enumerate() {
                    if i > 0 {
                        pairs = vpadalq_u8(pairs, *v);
                    }
                    totals[2 * i] = vaddw_u8(totals[2 * i], vget_low_u8(*v));
                    totals[2 * i + 1] = vaddw_high_u8(totals[2 * i + 1], *v);
                }
                v_s1 = vpadalq_u16(v_s1, pairs);
            }

            // s2 += 64·s2_lanes + Σ weight·total.
            v_s2 = vqshlq_n_u32::<6>(v_s2);
            for (total, weights) in totals.iter().zip(w.as_flattened()) {
                v_s2 = vmlal_u16(v_s2, vget_low_u16(*total), vget_low_u16(*weights));
                v_s2 = vmlal_high_u16(v_s2, *total, *weights);
            }
            s1 += vaddvq_u32(v_s1);
            s2 += vaddvq_u32(v_s2);
        }
        // The last 0 to 63 bytes, which also reduces the sums.
        scalar::chunk(&mut s1, &mut s2, c);
    }
    (s2 << 16) | s1
}

/// `acc + a·b`, four bytes into each 32-bit lane: `vdotq_u32`, whose intrinsic is not stable
/// in Rust 1.91, the oldest zflate builds with.
#[inline]
#[target_feature(enable = "neon,dotprod")]
fn udot(mut acc: uint32x4_t, a: uint8x16_t, b: uint8x16_t) -> uint32x4_t {
    // SAFETY: `udot` reads `a` and `b` and updates `acc`, nothing else; the function has the
    // dot-product feature, so the assembler accepts it.
    unsafe {
        std::arch::asm!(
            "udot {acc:v}.4s, {a:v}.16b, {b:v}.16b",
            acc = inout(vreg) acc,
            a = in(vreg) a,
            b = in(vreg) b,
            options(pure, nomem, nostack, preserves_flags),
        );
    }
    acc
}

/// NEON with the dot-product instructions.
#[target_feature(enable = "neon,dotprod")]
fn neon_dotprod(adler: u32, data: &[u8]) -> u32 {
    let w = load4(&WEIGHTS);
    let ones = vdupq_n_u8(1);
    let (mut s1, mut s2) = (adler & 0xffff, adler >> 16);
    let mut p = &data[align(&mut s1, &mut s2, data)..];

    while !p.is_empty() {
        let (chunk, rest) = p.split_at(p.len().min(CHUNK));
        p = rest;

        let (steps, c) = chunk.as_chunks::<64>();
        if !steps.is_empty() {
            // Four independent sets of sums, one per 16 bytes, for the latency of `udot`.
            let mut v_s1 = [vdupq_n_u32(0); 4];
            let mut v_s2 = [vdupq_n_u32(0); 4];
            let mut sums = [vdupq_n_u32(0); 4];
            s2 += s1 * (64 * steps.len()) as u32;

            for step in steps {
                let d = load4(step);
                for i in 0..4 {
                    sums[i] = vaddq_u32(sums[i], v_s1[i]);
                    v_s1[i] = udot(v_s1[i], d[i], ones);
                    v_s2[i] = udot(v_s2[i], d[i], w[i]);
                }
            }

            let add4 = |v: [uint32x4_t; 4]| vaddq_u32(vaddq_u32(v[0], v[1]), vaddq_u32(v[2], v[3]));
            // s2 += 64·sums + the weighted bytes.
            let v_s2 = vaddq_u32(add4(v_s2), vqshlq_n_u32::<6>(add4(sums)));
            s1 += vaddvq_u32(add4(v_s1));
            s2 += vaddvq_u32(v_s2);
        }
        // The last 0 to 63 bytes, which also reduces the sums.
        scalar::chunk(&mut s1, &mut s2, c);
    }
    (s2 << 16) | s1
}
