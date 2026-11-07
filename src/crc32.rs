//! The gzip CRC32 (ISO 3309, bit-reflected polynomial `0xedb88320`).
//!
//! The fast paths fold the data with carry-less multiplication, as libdeflate's
//! `crc32_pclmul_template.h` does; [`scalar`] is the fallback and handles the shortest inputs.
//! The best implementation for the CPU is chosen once, at the first call.

#[cfg(target_arch = "aarch64")]
mod arm;
mod scalar;
#[cfg(target_arch = "x86_64")]
mod x86;

use std::sync::OnceLock;

/// A CRC32 implementation, working on the inverted CRC (the register before and after the
/// final XOR), as the folding code wants it.
pub type Impl = fn(u32, &[u8]) -> u32;

/// The bit-reflected generator polynomial without its `x^32` term.
const POLY: u32 = 0xedb8_8320;

/// Returns `x^n mod G(x)`, bit-reflected: bit 31 is the coefficient of `x^0`.
///
/// These are the folding multipliers; computing them here, rather than copying them, keeps
/// them checkable.
#[cfg_attr(
    not(any(target_arch = "x86_64", target_arch = "aarch64")),
    allow(dead_code)
)]
const fn x_pow_mod_g(n: u32) -> u32 {
    let mut r = 1 << 31;
    let mut i = 0;
    while i < n {
        r = (r >> 1) ^ if r & 1 != 0 { POLY } else { 0 };
        i += 1;
    }
    r
}

/// Updates a CRC32 with `data`, starting from 0, as libdeflate's `libdeflate_crc32` does.
///
/// `crc32(crc32(0, a), b)` is the CRC of `a` followed by `b`.
#[must_use]
pub fn crc32(crc: u32, data: &[u8]) -> u32 {
    !best()(!crc, data)
}

fn best() -> Impl {
    static BEST: OnceLock<Impl> = OnceLock::new();
    *BEST.get_or_init(|| available().first().map_or(scalar::slice8, |&(_, f)| f))
}

/// The implementations this CPU supports, fastest first.
fn available() -> Vec<(&'static str, Impl)> {
    let mut impls = Vec::new();
    #[cfg(target_arch = "x86_64")]
    if !cfg!(feature = "force-portable") {
        impls.extend(x86::available());
    }
    // Not under Miri, which emulates neither PMULL nor EOR3: Apple's targets enable them
    // at compile time, so detection would pick them. (On x86-64 Miri detects no features.)
    #[cfg(target_arch = "aarch64")]
    if !cfg!(feature = "force-portable") && !cfg!(miri) {
        impls.extend(arm::available());
    }
    impls.push(("slice8", scalar::slice8 as Impl));
    impls
}

/// Each CRC32 implementation this CPU supports, by name, fastest first, for tests and
/// benchmarks. Not part of the stable API.
#[doc(hidden)]
pub fn implementations() -> Vec<(&'static str, Impl)> {
    available()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The multipliers libdeflate's `crc32_multipliers.h` lists, by exponent.
    #[test]
    fn multipliers_match_libdeflate() {
        for (n, want) in [
            (95, 0xccaa_009e),
            (159, 0xae68_9191),
            (223, 0x8125_6527),
            (287, 0xf1da_05aa),
            (479, 0x1d95_13d7),
            (543, 0x8f35_2d95),
            (991, 0x910e_eec1),
            (1055, 0x33ff_f533),
            (2015, 0xe95c_1271),
            (2079, 0xce33_71cb),
            (4063, 0x0c30_f51d),
            (4127, 0x1072_db28),
        ] {
            assert_eq!(x_pow_mod_g(n), want, "x^{n} mod G");
        }
    }

    #[test]
    fn check_value() {
        // The standard check value for CRC-32/ISO-HDLC.
        for (name, f) in available() {
            assert_eq!(!f(!0, b"123456789"), 0xcbf4_3926, "{name}");
        }
    }
}
