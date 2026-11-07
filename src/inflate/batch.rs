//! Many independent DEFLATE streams (BGZF blocks) decoded two at a time.
//!
//! One stream's fastloop is a chain of dependent lookups, each waiting on the last; a
//! second stream's chain is independent, so a core with room in its reorder window can run
//! both. Each round of [`fastloop_pair`] takes one step of each stream: the step is the
//! single fastloop's body, inlined twice, so each stream's branches have their own sites in
//! the predictor.
//!
//! [`inflate_batch`] keeps two lanes busy: when a lane's stream ends, the next stream in
//! the batch starts in it, so the lanes run alone only at the end of the batch. Headers,
//! table building and the generic loop near the ends of the buffers run one stream at a
//! time, as [`Stream::advance`]. A stream that meets an error is left for the caller to
//! decode again alone, so errors are exactly the single decoder's.

use super::*;

/// One stream's place between rounds. Its tables are its lane's, passed in.
struct Stream<'s> {
    bits: Bits<'s>,
    out: &'s mut [u8],
    pos: usize,
    state: State,
}

#[derive(Clone, Copy)]
enum State {
    /// The next block's header is next.
    Header,
    /// Inside a Huffman block.
    Block {
        is_final: bool,
        litlen_mask: u64,
    },
    Done,
}

impl<'s> Stream<'s> {
    fn new(t: &mut Tables, data: &'s [u8], out: &'s mut [u8]) -> Self {
        // As `inflate` decides, once per stream.
        let pairs = pairs_pay(data.len(), out.len());
        if pairs != t.pairs {
            t.static_loaded = false;
            t.pairs = pairs;
        }
        Self {
            bits: Bits::new(data),
            out,
            pos: 0,
            state: State::Header,
        }
    }

    /// Decodes until the fastloop can run, returning true, or until the stream ends,
    /// returning false.
    fn advance(&mut self, t: &mut Tables) -> Result<bool, DecompressionError> {
        loop {
            match self.state {
                State::Done => return Ok(false),
                State::Header => {
                    let block = start_block(t, &mut self.bits, self.out, &mut self.pos)?;
                    self.state = match block.huffman {
                        Some(litlen_mask) => State::Block {
                            is_final: block.is_final,
                            litlen_mask,
                        },
                        None if block.is_final => State::Done,
                        None => State::Header,
                    };
                }
                State::Block {
                    is_final,
                    litlen_mask,
                } => {
                    if fast_room(&self.bits, self.out.len(), self.pos) {
                        return Ok(true);
                    }
                    generic_loop(t, &mut self.bits, self.out, &mut self.pos, litlen_mask)?;
                    self.end_block(is_final);
                }
            }
        }
    }

    /// Decodes the rest alone.
    fn finish(&mut self, t: &mut Tables) -> Result<(), DecompressionError> {
        while self.advance(t)? {
            let litlen_mask = self.litlen_mask();
            if fastloop(t, &mut self.bits, self.out, &mut self.pos, litlen_mask)? {
                self.end_block(self.is_final());
            }
        }
        Ok(())
    }

    fn end_block(&mut self, is_final: bool) {
        self.state = if is_final { State::Done } else { State::Header };
    }

    fn litlen_mask(&self) -> u64 {
        match self.state {
            State::Block { litlen_mask, .. } => litlen_mask,
            _ => unreachable!("not in a Huffman block"),
        }
    }

    fn is_final(&self) -> bool {
        matches!(self.state, State::Block { is_final: true, .. })
    }

    fn inflated(&self) -> Result<Inflated, DecompressionError> {
        Ok(Inflated {
            consumed: self.bits.consumed()?,
            written: self.pos,
        })
    }
}

/// A lane: its tables, and the stream in it with that stream's index in the batch.
struct Lane<'t, 's> {
    t: &'t mut Tables,
    stream: Option<(usize, Stream<'s>)>,
}

impl<'s> Lane<'_, 's> {
    /// Makes the lane's stream ready for the fastloop, recording streams that end or fail
    /// and starting the next from `queue`, until one is ready or the queue is empty.
    fn fill(
        &mut self,
        queue: &mut impl Iterator<Item = (usize, (&'s [u8], &'s mut [u8]))>,
        done: &mut [Option<Inflated>],
    ) {
        loop {
            if let Some((i, s)) = &mut self.stream {
                match s.advance(self.t) {
                    Ok(true) => return,
                    Ok(false) => done[*i] = s.inflated().ok(),
                    Err(_) => {}
                }
                self.stream = None;
            }
            let Some((i, (data, out))) = queue.next() else {
                return;
            };
            self.stream = Some((i, Stream::new(self.t, data, out)));
        }
    }
}

/// Inflates each of `data` into the same-numbered `out`, two at a time, as calls of
/// [`inflate`] would, and records each stream's result in `done`, which must start as
/// `None`s, one per stream. A stream left at `None` met an error: its output then holds
/// anything, and it must be decoded again alone.
pub(crate) fn inflate_batch<'s, 'd: 's>(
    tables: [&mut Tables; 2],
    data: &[&'d [u8]],
    out: &'s mut [&mut [u8]],
    done: &mut [Option<Inflated>],
) {
    let mut queue = data
        .iter()
        .zip(out.iter_mut())
        .map(|(d, o)| (*d, &mut **o))
        .enumerate();
    let [ta, tb] = tables;
    let mut a = Lane {
        t: ta,
        stream: None,
    };
    let mut b = Lane {
        t: tb,
        stream: None,
    };
    loop {
        a.fill(&mut queue, done);
        b.fill(&mut queue, done);
        match (&mut a.stream, &mut b.stream) {
            (Some((_, sa)), Some((_, sb))) => {
                // SAFETY: a lane's stream is ready for the fastloop whenever `fill` leaves
                // one in it, and each lane's stream uses that lane's tables.
                if let Err(failed) = unsafe { fastloop_pair(sa, a.t, sb, b.t) } {
                    // The other stream's place is saved: it goes on in the next round.
                    match failed {
                        Failed::A => a.stream = None,
                        Failed::B => b.stream = None,
                    }
                }
            }
            (Some(_), None) | (None, Some(_)) => {
                // The queue is empty: the last stream finishes alone.
                let lane = if a.stream.is_some() { &mut a } else { &mut b };
                let (i, mut s) = lane.stream.take().expect("checked");
                if s.finish(lane.t).is_ok() {
                    done[i] = s.inflated().ok();
                }
            }
            (None, None) => return,
        }
    }
}

/// Which of [`fastloop_pair`]'s streams met an error.
enum Failed {
    A,
    B,
}

/// Both streams' fastloops, a step of each per round, until either reaches the end of its
/// block or of its fast region, or meets an error. Saves both streams' places, except that
/// of one that failed, which is to be given up.
///
/// # Safety
///
/// Both streams must be ready for the fastloop (their last [`Stream::advance`] returned
/// true), and `ta` and `tb` must be their tables.
#[inline(never)]
unsafe fn fastloop_pair(
    sa: &mut Stream<'_>,
    ta: &Tables,
    sb: &mut Stream<'_>,
    tb: &Tables,
) -> Result<(), Failed> {
    // SAFETY: as the caller promises: ready means inside a Huffman block with `fast_room`.
    let (mut fa, mut fb) = unsafe {
        (
            Fast::new(ta, &sa.bits, sa.out, sa.pos, sa.litlen_mask()),
            Fast::new(tb, &sb.bits, sb.out, sb.pos, sb.litlen_mask()),
        )
    };
    let (mut end_a, mut end_b) = (false, false);
    while fa.room() && fb.room() {
        // SAFETY: `room` holds for both.
        unsafe {
            match fa.step() {
                Ok(Step::Go) => {}
                Ok(Step::EndOfBlock) => {
                    end_a = true;
                    break;
                }
                Err(_) => {
                    fb.save(&mut sb.bits, &mut sb.pos);
                    return Err(Failed::A);
                }
            }
            match fb.step() {
                Ok(Step::Go) => {}
                Ok(Step::EndOfBlock) => {
                    end_b = true;
                    break;
                }
                Err(_) => {
                    fa.save(&mut sa.bits, &mut sa.pos);
                    return Err(Failed::B);
                }
            }
        }
    }
    fa.save(&mut sa.bits, &mut sa.pos);
    fb.save(&mut sb.bits, &mut sb.pos);
    if end_a {
        sa.end_block(sa.is_final());
    }
    if end_b {
        sb.end_block(sb.is_final());
    }
    Ok(())
}
