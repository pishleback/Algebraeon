use crate::structure::*;
use algebraeon_nzq::{Integer, Natural, gcd, lcm, traits::Abs};
use algebraeon_sets::structure::CannonicalStructure;

impl IdealStructure for CannonicalStructure<Integer> {
    type Ideal = Natural;
}

impl IdealArithmeticStructure for CannonicalStructure<Integer> {
    fn principal_ideal(&self, a: &Self::Set) -> Self::Ideal {
        a.abs()
    }

    fn ideal_equal(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool {
        a == b
    }

    fn ideal_contains(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool {
        b % a == Natural::ZERO
    }

    fn ideal_intersection(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        lcm(a.clone(), b.clone())
    }

    fn ideal_add(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        gcd(a.clone(), b.clone())
    }

    fn ideal_mul(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        a * b
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_ideals() {
        assert!(Integer::ideal_equal(
            &Natural::from(3u32),
            &Natural::from(3u32)
        ));

        assert!(!Integer::ideal_equal(
            &Natural::from(2u32),
            &Natural::from(3u32)
        ));

        assert!(Integer::ideal_contains(
            &Natural::from(2u32),
            &Natural::from(6u32)
        ));

        assert!(!Integer::ideal_contains(
            &Natural::from(6u32),
            &Natural::from(2u32)
        ));

        assert!(!Integer::ideal_contains(
            &Natural::from(5u32),
            &Natural::from(7u32)
        ));

        assert!(Integer::ideal_equal(
            &Integer::ideal_add(&Natural::from(6u32), &Natural::from(15u32)),
            &Natural::from(3u32)
        ));

        assert!(Integer::ideal_equal(
            &Integer::ideal_intersection(&Natural::from(6u32), &Natural::from(15u32)),
            &Natural::from(30u32)
        ));

        assert!(Integer::ideal_equal(
            &Integer::ideal_mul(&Natural::from(6u32), &Natural::from(15u32)),
            &Natural::from(90u32)
        ));

        assert_eq!(Integer::generated_ideal(vec![-15, 6]), 3u32.into());
    }
}
