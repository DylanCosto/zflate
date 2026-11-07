//! Allocations that report failure as `None` instead of aborting, so that the constructors
//! can return an error under an address-space limit (`ulimit -v`), as libdeflate's do.

use std::{alloc::Layout, ptr::NonNull};

/// `N` copies of `value` in a box, built in place rather than on the stack, or `None` if
/// it cannot be allocated.
pub(crate) fn try_boxed_array<T: Clone, const N: usize>(value: T) -> Option<Box<[T; N]>> {
    let layout = Layout::new::<[T; N]>();
    if layout.size() == 0 {
        return Some(Box::new(std::array::from_fn(|_| value.clone())));
    }
    // SAFETY: the layout's size is not zero.
    let ptr = NonNull::new(unsafe { std::alloc::alloc(layout) }.cast::<T>())?;
    for i in 0..N {
        // SAFETY: the allocation holds `N` `T`s, and `i < N`.
        unsafe { ptr.as_ptr().add(i).write(value.clone()) };
    }
    // SAFETY: the allocation, from the global allocator with `[T; N]`'s layout, holds `N`
    // initialized `T`s, so `Box` may own it.
    Some(unsafe { Box::from_raw(ptr.as_ptr().cast::<[T; N]>()) })
}

/// `Box::new(value)`, or `None` if the box cannot be allocated.
pub(crate) fn try_box<T>(value: T) -> Option<Box<T>> {
    let layout = Layout::new::<T>();
    if layout.size() == 0 {
        return Some(Box::new(value));
    }
    // SAFETY: the layout's size is not zero.
    let ptr = NonNull::new(unsafe { std::alloc::alloc(layout) }.cast::<T>())?;
    // SAFETY: `ptr` is a new allocation from the global allocator with `T`'s layout, so
    // writing `value` initializes it, and `Box` may own it.
    unsafe {
        ptr.as_ptr().write(value);
        Some(Box::from_raw(ptr.as_ptr()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn allocations_hold_their_values() {
        assert_eq!(
            try_boxed_array::<_, 4>([1i16, 2]).as_deref(),
            Some(&[[1, 2]; 4])
        );
        assert_eq!(try_box([5u32; 1000]).as_deref(), Some(&[5; 1000]));
        assert_eq!(try_box(()), Some(Box::new(())));
        assert_eq!(try_boxed_array::<u8, 0>(1).as_deref(), Some(&[]));
        assert_eq!(try_boxed_array::<(), 3>(()).as_deref(), Some(&[(); 3]));
    }
}
