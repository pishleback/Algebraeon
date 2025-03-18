use crate::structure::structure::*;
use algebraeon_sets::structure::*;

pub mod factor;
pub mod functions;
pub mod primes;

use algebraeon_nzq::integer::*;
use algebraeon_nzq::natural::*;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nat_to_usize() {
        assert_eq!(nat_to_usize(&Natural::from(0u8)).unwrap(), 0);
        assert_eq!(nat_to_usize(&Natural::from(1u8)).unwrap(), 1);
        assert_eq!(nat_to_usize(&Natural::from(2u8)).unwrap(), 2);

        // use malachite_base::num::arithmetic::traits::ModPow;

        // let a = malachite_nz::natural::Natural::from(1u8);
        // let c = a.mod_pow(a);
    }
}
