use std::collections::HashMap;

use malachite_base::num::basic::traits::{One, Two, Zero};
use malachite_nz::natural::Natural;

pub mod factor;
pub mod functions;
pub mod primes;

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
