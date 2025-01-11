use std::collections::HashMap;

use algebraeon_sets::structure::*;
use malachite_base::num::basic::traits::{One, Two, Zero};
use malachite_nz::natural::Natural;

use crate::structure::structure::*;

pub mod factor;
pub mod functions;
pub mod primes;

impl SemiRingStructure for CannonicalStructure<Natural> {
    fn zero(&self) -> Self::Set {
        Natural::ZERO
    }
    fn one(&self) -> Self::Set {
        Natural::ONE
    }
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }
    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
    }
}

pub fn nat_to_usize(n: &Natural) -> Result<usize, ()> {
    let limbs = n.to_limbs_asc();
    if limbs.len() == 0 {
        Ok(0)
    } else if limbs.len() == 1 {
        let n = limbs[0];
        if Natural::from(n) > Natural::from(usize::MAX) {
            Err(())
        } else {
            Ok(n as usize)
        }
    } else {
        Err(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nat_to_usize() {
        assert_eq!(nat_to_usize(&Natural::from(0u8)).unwrap(), 0);
        assert_eq!(nat_to_usize(&Natural::from(1u8)).unwrap(), 1);
        assert_eq!(nat_to_usize(&Natural::from(2u8)).unwrap(), 2);
    }
}
