//! For generating pseudo-random numbers.

use std::borrow::Borrow;

use super::natural::*;

/// An object for generating pseudo-random numbers.
pub struct Rng {
    rng: malachite_base::num::random::RandomPrimitiveInts<u64>,
}

impl Rng {
    /// Constructor.
    pub fn new(seed: u128) -> Self {
        let mut bytes = [0u8; 32]; // Initialize an array of 32 zeroed bytes
        bytes[16..].copy_from_slice(&seed.to_be_bytes());
        let seed = malachite_base::random::Seed::from_bytes(bytes);
        Self {
            rng: malachite_base::num::random::random_primitive_ints(seed),
        }
    }

    /// Return a natural number randomly and uniformly from the interval [0, self]
    pub fn random_below(&mut self, n: impl Borrow<Natural>) -> Natural {
        Natural::from_malachite(malachite_nz::natural::random::get_random_natural_less_than(
            &mut self.rng,
            (n.borrow() + Natural::ONE).to_malachite_ref(),
        ))
    }
}
