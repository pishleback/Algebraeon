use crate::{
    rings::natural::factorization::{
        NaturalCanonicalFactorizationStructure, factor, primes::is_prime,
    },
    structure::*,
};
use algebraeon_nzq::{traits::Abs, *};
use algebraeon_sets::structure::{BorrowedStructure, MetaType, SetSignature, Signature};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntegerIdealsStructure<B: BorrowedStructure<IntegerCanonicalStructure>> {
    integers: B,
}

impl CanonicalIdealsSignature for IntegerCanonicalStructure {
    type Ideals<SelfB: BorrowedStructure<IntegerCanonicalStructure>> =
        IntegerIdealsStructure<SelfB>;

    fn ideals<'a>(&'a self) -> Self::Ideals<&'a Self> {
        IntegerIdealsStructure { integers: self }
    }

    fn into_ideals(self) -> Self::Ideals<Self> {
        IntegerIdealsStructure { integers: self }
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> Signature for IntegerIdealsStructure<B> {}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> SetSignature for IntegerIdealsStructure<B> {
    type Set = Natural;
    fn is_element(&self, _x: &Self::Set) -> bool {
        true
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> IdealsSignature<IntegerCanonicalStructure, B>
    for IntegerIdealsStructure<B>
{
    fn ring(&self) -> &IntegerCanonicalStructure {
        self.integers.borrow()
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>>
    IdealsArithmeticSignature<IntegerCanonicalStructure, B> for IntegerIdealsStructure<B>
{
    fn principal_ideal(&self, a: &Integer) -> Self::Set {
        a.abs()
    }

    fn ideal_equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        a == b
    }

    fn ideal_contains(&self, a: &Self::Set, b: &Self::Set) -> bool {
        b % a == Natural::ZERO
    }

    fn ideal_intersect(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        lcm(a.clone(), b.clone())
    }

    fn ideal_add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        gcd(a.clone(), b.clone())
    }

    fn ideal_mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a * b
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>>
    PrincipalIdealsSignature<IntegerCanonicalStructure, B> for IntegerIdealsStructure<B>
{
    fn ideal_generator(&self, ideal: &Natural) -> Integer {
        Integer::from(ideal)
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>>
    DedekindDomainIdealsSignature<IntegerCanonicalStructure, B> for IntegerIdealsStructure<B>
{
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>>
    FactorableIdealsSignature<IntegerCanonicalStructure, B> for IntegerIdealsStructure<B>
{
    fn factor_ideal(
        &self,
        ideal: &Self::Set,
    ) -> Option<DedekindDomainIdealFactorization<Self::Set>> {
        let f = factor(ideal.clone())?;
        Some(
            self.factorizations().new_powers(
                Natural::structure()
                    .factorizations()
                    .into_powers(f)
                    .into_iter()
                    .map(|(n, k)| (DedekindDomainPrimeIdeal::from_ideal_unchecked(n), k.into()))
                    .collect(),
            ),
        )
    }
}

impl DedekindDomainPrimeIdeal<Natural> {
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
        assert!(Integer::ideals().ideal_equal(&Natural::from(3u32), &Natural::from(3u32)));

        assert!(!Integer::ideals().ideal_equal(&Natural::from(2u32), &Natural::from(3u32)));

        assert!(Integer::ideals().ideal_contains(&Natural::from(2u32), &Natural::from(6u32)));

        assert!(!Integer::ideals().ideal_contains(&Natural::from(6u32), &Natural::from(2u32)));

        assert!(!Integer::ideals().ideal_contains(&Natural::from(5u32), &Natural::from(7u32)));

        assert!(Integer::ideals().ideal_equal(
            &Integer::ideals().ideal_add(&Natural::from(6u32), &Natural::from(15u32)),
            &Natural::from(3u32)
        ));

        assert!(Integer::ideals().ideal_equal(
            &Integer::ideals().ideal_intersect(&Natural::from(6u32), &Natural::from(15u32)),
            &Natural::from(30u32)
        ));

        assert!(Integer::ideals().ideal_equal(
            &Integer::ideals().ideal_mul(&Natural::from(6u32), &Natural::from(15u32)),
            &Natural::from(90u32)
        ));

        assert_eq!(Integer::ideals().generated_ideal(vec![-15, 6]), 3u32.into());
    }

    #[test]
    fn factor_integer_ideal() {
        let f =
            Integer::ideals().factor_ideal(&Integer::ideals().principal_ideal(&Integer::from(0)));
        println!("{:?}", f);
        assert!(f.is_none());

        let f =
            Integer::ideals().factor_ideal(&Integer::ideals().principal_ideal(&Integer::from(18)));
        println!("{:?}", f);
        assert!(f.is_some());
    }
}
