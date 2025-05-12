// MIT License

// Copyright 2024 Evan Schwartz

// Permission is hereby granted, free of charge, to any person obtaining a copy of this
// software and associated documentation files (the "Software"), to deal in the Software
// without restriction, including without limitation the rights to use, copy, modify,
// merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
// INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
// PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

use std::os::unix::raw::uid_t;

use crate::distances;

#[inline(always)]
pub fn hamming_bitwise_fast(x: &[u8], y: &[u8]) -> u32 {
    assert_eq!(x.len(), y.len());

    // Process 8 bytes at a time using u64
    let mut distance = x
        .chunks_exact(8)
        .zip(y.chunks_exact(8))
        .map(|(x_chunk, y_chunk)| {
            // This is safe because we know the chunks are exactly 8 bytes.
            // Also, we don't care whether the platform uses little-endian or big-endian
            // byte order. Since we're only XORing values, we just care that the
            // endianness is the same for both.
            let x_val = u64::from_ne_bytes(x_chunk.try_into().unwrap());
            let y_val = u64::from_ne_bytes(y_chunk.try_into().unwrap());
            (x_val ^ y_val).count_ones()
        })
        .sum::<u32>();

    if x.len() % 8 != 0 {
        distance += x
            .chunks_exact(8)
            .remainder()
            .iter()
            .zip(y.chunks_exact(8).remainder())
            .map(|(x_byte, y_byte)| (x_byte ^ y_byte).count_ones())
            .sum::<u32>();
    }

    distance
}


pub fn jaccard_distance_fast(x: &[u8], y: &[u8]) -> (u32, u32) {
    assert_eq!(x.len(), y.len());

    let mut intersection = 0;
    let mut union = 0;

    // Process 8 bytes at a time using u64
    for (x_chunk, y_chunk) in x.chunks_exact(8).zip(y.chunks_exact(8)) {
        let x_val = u64::from_ne_bytes(x_chunk.try_into().unwrap());
        let y_val = u64::from_ne_bytes(y_chunk.try_into().unwrap());

        intersection += (x_val & y_val).count_ones();
        union += (x_val | y_val).count_ones();
    }

    // Handle remaining bytes
    for (x_byte, y_byte) in x.chunks_exact(8).remainder().iter().zip(y.chunks_exact(8).remainder()) {
        intersection += (x_byte & y_byte).count_ones();
        union += (x_byte | y_byte).count_ones();
    }

    (intersection, union)
}