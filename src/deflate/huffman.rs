//! Length-limited canonical Huffman codes, built as libdeflate builds them (after 7-Zip's
//! `HuffEnc.c`): sort the symbols by frequency, build only the tree's inner nodes, in place,
//! from the sorted leaves, read off how many codewords each length gets, clamping depths to
//! the limit, and hand out canonical codewords in symbol order.

use super::tables::NUM_LITLEN_SYMS;

/// Each working entry holds a symbol in its low bits and a frequency, or later a parent index
/// or a depth, above them.
const NUM_SYMBOL_BITS: u32 = 10;
const SYMBOL_MASK: u32 = (1 << NUM_SYMBOL_BITS) - 1;
const FREQ_MASK: u32 = !SYMBOL_MASK;
/// The frequencies of a code may add up to this much at most.
pub(super) const MAX_FREQ: u32 = (1 << (32 - NUM_SYMBOL_BITS)) - 1;

/// Builds a canonical Huffman code for symbols with the frequencies `freqs`, whose codewords
/// are at most `max_len` bits, into `lens` and `codewords` (bit-reversed, as DEFLATE writes
/// them). Unused symbols get length 0. A code always has at least two codewords, as other
/// compressors' do and some decoders require: with fewer symbols in use, symbols 0 and 1 (or
/// the used one) get 1-bit codewords.
pub(super) fn make_code(max_len: u32, freqs: &[u32], lens: &mut [u8], codewords: &mut [u32]) {
    let num_syms = freqs.len();
    assert!(
        (2..=NUM_LITLEN_SYMS).contains(&num_syms)
            && lens.len() == num_syms
            && codewords.len() == num_syms
            && (1..=15).contains(&max_len),
        "{num_syms} symbols, {max_len}-bit codewords"
    );
    debug_assert!(freqs.iter().map(|&f| u64::from(f)).sum::<u64>() <= u64::from(MAX_FREQ));
    // The working array shares the codewords' storage, as in libdeflate.
    let a = codewords;
    let num_used = sort_symbols(freqs, lens, a);
    if num_used < 2 {
        let sym = if num_used == 1 { a[0] & SYMBOL_MASK } else { 0 };
        let other = if sym == 0 { 1 } else { sym as usize };
        a[0] = 0;
        lens[0] = 1;
        a[other] = 1;
        lens[other] = 1;
        return;
    }
    build_tree(&mut a[..num_used]);
    let mut len_counts = [0u32; 16];
    compute_length_counts(&mut a[..num_used - 1], &mut len_counts, max_len);
    gen_codewords(a, lens, &len_counts, max_len, num_used);
}

/// Sorts the used symbols by frequency, then by symbol, into `a` as `sym | freq << 10`;
/// sets the lengths of unused symbols to 0, and returns how many symbols are used. A counting
/// sort does most of the work, since most frequencies are small; the symbols in its last
/// bucket, the most frequent, are then sorted fully.
fn sort_symbols(freqs: &[u32], lens: &mut [u8], a: &mut [u32]) -> usize {
    let num_syms = freqs.len();
    let mut counters = [0u32; NUM_LITLEN_SYMS];
    let counters = &mut counters[..num_syms];
    let bucket = |freq: u32| (freq as usize).min(num_syms - 1);
    for &freq in freqs {
        counters[bucket(freq)] += 1;
    }
    // Make the counters the buckets' starts, skipping bucket 0, the unused symbols.
    let mut num_used = 0;
    for counter in &mut counters[1..] {
        let count = *counter;
        *counter = num_used;
        num_used += count;
    }
    for (sym, &freq) in freqs.iter().enumerate() {
        if freq == 0 {
            lens[sym] = 0;
        } else {
            debug_assert!(freq <= MAX_FREQ, "frequency {freq}");
            let slot = &mut counters[bucket(freq)];
            a[*slot as usize] = sym as u32 | freq << NUM_SYMBOL_BITS;
            *slot += 1;
        }
    }
    // The counters are now the buckets' ends: the last bucket starts where the one before
    // it ends. Its entries are distinct, so any sort gives libdeflate's order.
    a[counters[num_syms - 2] as usize..counters[num_syms - 1] as usize].sort_unstable();
    num_used as usize
}

/// Builds the inner nodes of the Huffman tree over `a`'s sorted leaves, in `a` itself:
/// `a[e]` for `e < a.len() - 1` becomes an inner node, the last of them the root, each
/// holding its parent's index (or, while unparented, its frequency) above the symbol bits,
/// which stay as they were. Leaves and unparented inner nodes are each kept in increasing
/// order of frequency, so the two lowest are always among the next two of each.
fn build_tree(a: &mut [u32]) {
    let last = a.len() - 1;
    // The next leaf without a parent, the next inner node without one, and the next inner
    // node to make.
    let (mut i, mut b, mut e) = (0, 0, 0);
    loop {
        let freq = |x: u32| x & FREQ_MASK;
        let new_freq;
        if i < last && (b == e || freq(a[i + 1]) <= freq(a[b])) {
            // Two leaves.
            new_freq = freq(a[i]) + freq(a[i + 1]);
            i += 2;
        } else if b + 2 <= e && (i > last || freq(a[b + 1]) < freq(a[i])) {
            // Two inner nodes.
            new_freq = freq(a[b]) + freq(a[b + 1]);
            a[b] = (e as u32) << NUM_SYMBOL_BITS | (a[b] & SYMBOL_MASK);
            a[b + 1] = (e as u32) << NUM_SYMBOL_BITS | (a[b + 1] & SYMBOL_MASK);
            b += 2;
        } else {
            // A leaf and an inner node.
            new_freq = freq(a[i]) + freq(a[b]);
            a[b] = (e as u32) << NUM_SYMBOL_BITS | (a[b] & SYMBOL_MASK);
            i += 1;
            b += 1;
        }
        a[e] = new_freq | (a[e] & SYMBOL_MASK);
        e += 1;
        // A tree with n leaves has n - 1 inner nodes.
        if e >= last {
            break;
        }
    }
}

/// Counts the codewords of each length from the inner nodes `nodes` that [`build_tree`]
/// made (the root last), overwriting the parent indices with depths. It starts from two
/// codewords of length 1, the root's children, and for each inner node, visited parents
/// first, turns one codeword of its depth into two one deeper. A node deeper than the limit
/// is moved up to the deepest length that still has a codeword to split, which keeps the
/// code complete; it is not optimal, but it is close and rarely needed.
fn compute_length_counts(nodes: &mut [u32], len_counts: &mut [u32; 16], max_len: u32) {
    let root = nodes.len() - 1;
    len_counts.fill(0);
    len_counts[1] = 2;
    nodes[root] &= SYMBOL_MASK;
    for node in (0..root).rev() {
        let parent = (nodes[node] >> NUM_SYMBOL_BITS) as usize;
        let parent_depth = nodes[parent] >> NUM_SYMBOL_BITS;
        let mut depth = parent_depth + 1;
        nodes[node] = (nodes[node] & SYMBOL_MASK) | depth << NUM_SYMBOL_BITS;
        if depth >= max_len {
            depth = max_len;
            loop {
                depth -= 1;
                if len_counts[depth as usize] != 0 {
                    break;
                }
            }
        }
        len_counts[depth as usize] -= 1;
        len_counts[depth as usize + 1] += 2;
    }
}

/// Gives the lengths out in decreasing order to the symbols in increasing order of
/// frequency (the first `num_used` entries of `a`), then assigns canonical codewords in
/// symbol order, overwriting `a`.
fn gen_codewords(
    a: &mut [u32],
    lens: &mut [u8],
    len_counts: &[u32; 16],
    max_len: u32,
    num_used: usize,
) {
    let mut i = 0;
    for len in (1..=max_len).rev() {
        for _ in 0..len_counts[len as usize] {
            lens[(a[i] & SYMBOL_MASK) as usize] = len as u8;
            i += 1;
        }
    }
    debug_assert_eq!(i, num_used);
    assign_codewords(lens, len_counts, a);
}

/// Assigns the canonical code with lengths `lens` (0 for unused symbols) to `codewords`,
/// bit-reversed.
pub(super) fn codewords_from_lens(lens: &[u8], codewords: &mut [u32]) {
    let mut len_counts = [0u32; 16];
    for &len in lens {
        len_counts[usize::from(len)] += 1;
    }
    assign_codewords(lens, &len_counts, codewords);
}

/// [`codewords_from_lens`] with the number of codewords of each length known.
fn assign_codewords(lens: &[u8], len_counts: &[u32; 16], codewords: &mut [u32]) {
    let mut next = [0u32; 16];
    for len in 2..16 {
        next[len] = (next[len - 1] + len_counts[len - 1]) << 1;
    }
    // No branch on unused symbols, which come and go unpredictably: they take `next[0]`,
    // which nothing else uses, and get 0.
    for (codeword, &len) in codewords.iter_mut().zip(lens) {
        let len = usize::from(len) & 15;
        *codeword = reverse_codeword(next[len], len as u32);
        next[len] += 1;
    }
}

/// `codeword`'s low `len` bits (at most 16), reversed; 0 for an unused symbol. x86 has no
/// bit-reverse instruction, so `reverse_bits` takes a dozen; a byte table takes two loads,
/// as in libdeflate.
#[inline]
pub(super) fn reverse_codeword(codeword: u32, len: u32) -> u32 {
    debug_assert!(len <= 16 && (len == 0 || codeword >> len == 0));
    let reversed = u32::from(REVERSED_BYTES[(codeword & 0xff) as usize]) << 8
        | u32::from(REVERSED_BYTES[((codeword >> 8) & 0xff) as usize]);
    reversed.checked_shr(16 - len).unwrap_or(0) & ((1 << len) - 1)
}

/// Each byte with its bits reversed.
static REVERSED_BYTES: [u8; 256] = {
    let mut table = [0; 256];
    let mut i = 0;
    while i < 256 {
        table[i] = (i as u8).reverse_bits();
        i += 1;
    }
    table
};

#[cfg(test)]
mod tests {
    use super::*;

    /// Checks that `lens` is a complete prefix code within `max_len` and that the codewords
    /// are canonical.
    fn check(freqs: &[u32], max_len: u32, lens: &[u8], codewords: &[u32]) {
        let used = freqs.iter().filter(|&&f| f != 0).count();
        let mut kraft = 0u64;
        for (sym, (&len, &freq)) in lens.iter().zip(freqs).enumerate() {
            assert!(u32::from(len) <= max_len, "symbol {sym}: length {len}");
            if used >= 2 {
                assert_eq!(len == 0, freq == 0, "symbol {sym}");
            }
            if len != 0 {
                kraft += 1 << (15 - len);
            }
        }
        assert_eq!(kraft, 1 << 15, "the code is not complete: {lens:?}");
        // Canonical: the codewords, unreversed, increase with (length, symbol).
        let mut syms: Vec<usize> = (0..lens.len()).filter(|&s| lens[s] != 0).collect();
        syms.sort_by_key(|&s| (lens[s], s));
        let natural = |s: usize| reverse_codeword(codewords[s], u32::from(lens[s]));
        for w in syms.windows(2) {
            let (a, b) = (w[0], w[1]);
            let (la, lb) = (u32::from(lens[a]), u32::from(lens[b]));
            assert!(natural(a) << (lb - la) < natural(b) || (la == lb && natural(a) < natural(b)));
        }
    }

    fn build(freqs: &[u32], max_len: u32) -> (Vec<u8>, Vec<u32>) {
        let mut lens = vec![0; freqs.len()];
        let mut codewords = vec![0; freqs.len()];
        make_code(max_len, freqs, &mut lens, &mut codewords);
        check(freqs, max_len, &lens, &codewords);
        (lens, codewords)
    }

    #[test]
    fn few_symbols_get_two_codewords() {
        assert_eq!(build(&[0; 30], 15).0[..3], [1, 1, 0]);
        assert_eq!(build(&[0, 0, 7, 0], 15).0, [1, 0, 1, 0]);
        assert_eq!(build(&[9, 0, 0, 0], 15).0, [1, 1, 0, 0]);
    }

    #[test]
    fn matches_textbook_huffman() {
        // Frequencies 1, 1, 2, 4, 8: lengths 4, 4, 3, 2, 1.
        assert_eq!(build(&[8, 4, 2, 1, 1], 15).0, [1, 2, 3, 4, 4]);
    }

    #[test]
    fn limits_lengths() {
        // Fibonacci frequencies make the deepest unlimited tree.
        let mut freqs = vec![1u32, 1];
        while freqs.len() < 30 {
            freqs.push(freqs[freqs.len() - 1] + freqs[freqs.len() - 2]);
        }
        for max_len in [7, 14, 15] {
            build(&freqs, max_len);
        }
        let (lens, _) = build(&freqs[..19], 7);
        assert_eq!(lens.iter().max(), Some(&7));
    }

    #[test]
    fn random_frequencies() {
        let mut x = 0x2545_f491_4f6c_dd1du64;
        for round in 0..2000 {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            let num_syms = [19, 30, 32, 288][round % 4];
            let max_len = [7, 15, 15, 14][round % 4];
            let freqs: Vec<u32> = (0..num_syms)
                .map(|i| {
                    let r = x.rotate_left(i as u32 * 7);
                    if r & 3 == 0 {
                        0
                    } else {
                        (r >> 8) as u32 % [2, 50, 10_000][round % 3]
                    }
                })
                .collect();
            build(&freqs, max_len);
        }
    }
}
