use crate::{
    number::natural::factorization::{factor, primes::is_prime},
    structure::*,
};
use algebraeon_nzq::{traits::Abs, *};
use algebraeon_sets::structure::MetaType;

impl IdealStructure for IntegerCanonicalStructure {
    type Ideal = Natural;
}

impl IdealArithmeticStructure for IntegerCanonicalStructure {
    fn principal_ideal(&self, a: &Self::Set) -> Self::Ideal {
        a.abs()
    }

    fn ideal_equal(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool {
        a == b
    }

    fn ideal_contains(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool {
        b % a == Natural::ZERO
    }

    fn ideal_intersect(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        lcm(a.clone(), b.clone())
    }

    fn ideal_add(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        gcd(a.clone(), b.clone())
    }

    fn ideal_mul(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        a * b
    }
}

impl PrincipalIdealDomainStructure for IntegerCanonicalStructure {
    fn ideal_generator(&self, ideal: &Self::Ideal) -> Self::Set {
        Integer::from(ideal)
    }
}

impl DedekindDomainStructure for IntegerCanonicalStructure {}

impl FactorableIdealsStructure for IntegerCanonicalStructure {
    fn factor_ideal(&self, ideal: &Self::Ideal) -> Option<DedekindDomainIdealFactorization<Self>> {
        let f = factor(ideal.clone())?;
        Some(DedekindDomainIdealFactorization::from_factor_powers(
            Integer::structure(),
            f.into_factor_powers()
                .into_iter()
                .map(|(n, k)| (DedekindDomainPrimeIdeal::from_ideal_unchecked(n), k.into()))
                .collect(),
        ))
    }
}

impl DedekindDomainPrimeIdeal<IntegerCanonicalStructure> {
    pub fn try_from_nat(n: Natural) -> Result<Self, ()> {
        if is_prime(&n) {
            Ok(DedekindDomainPrimeIdeal::from_ideal_unchecked(n))
        } else {
            Err(())
        }
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
            &Integer::ideal_intersect(&Natural::from(6u32), &Natural::from(15u32)),
            &Natural::from(30u32)
        ));

        assert!(Integer::ideal_equal(
            &Integer::ideal_mul(&Natural::from(6u32), &Natural::from(15u32)),
            &Natural::from(90u32)
        ));

        assert_eq!(Integer::generated_ideal(vec![-15, 6]), 3u32.into());
    }

    #[test]
    fn factor_integer_ideal() {
        let f = Integer::factor_ideal(&Integer::from(0).principal_ideal());
        println!("{:?}", f);
        assert!(f.is_none());

        let f = Integer::factor_ideal(&Integer::from(18).principal_ideal());
        println!("{:?}", f);
        assert!(f.is_some());
    }
}
