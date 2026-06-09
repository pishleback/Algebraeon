use crate::structure::*;
use algebraeon_nzq::{traits::Abs, *};
use algebraeon_structures::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntegerIdealsStructure<B: BorrowedStructure<IntegerCanonicalStructure>> {
    integers: B,
}

impl RingToIdealsSignature for IntegerCanonicalStructure {
    type Ideals<SelfB: BorrowedStructure<IntegerCanonicalStructure>> =
        IntegerIdealsStructure<SelfB>;

    fn ideals(&self) -> Self::Ideals<&Self> {
        IntegerIdealsStructure { integers: self }
    }

    fn into_ideals(self) -> Self::Ideals<Self> {
        IntegerIdealsStructure { integers: self }
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> Signature for IntegerIdealsStructure<B> {}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> SetSignature for IntegerIdealsStructure<B> {
    type Elem = Natural;
    fn validate_element(&self, _x: &Self::Elem) -> Result<(), String> {
        Ok(())
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> IdealsSignature<IntegerCanonicalStructure, B>
    for IntegerIdealsStructure<B>
{
    fn ring(&self) -> &IntegerCanonicalStructure {
        self.integers.borrow()
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> EqSignature for IntegerIdealsStructure<B> {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        a == b
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> RinglikeSpecializationSignature
    for IntegerIdealsStructure<B>
{
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> ZeroSignature for IntegerIdealsStructure<B> {
    fn zero(&self) -> Self::Elem {
        Natural::ZERO
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> AdditionSignature
    for IntegerIdealsStructure<B>
{
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        gcd(a.clone(), b.clone())
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> AdditiveMonoidSignature
    for IntegerIdealsStructure<B>
{
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> TryNegateSignature
    for IntegerIdealsStructure<B>
{
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        if self.is_zero(a) {
            Some(self.zero())
        } else {
            None
        }
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> OneSignature for IntegerIdealsStructure<B> {
    fn one(&self) -> Self::Elem {
        Natural::ONE
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> MultiplicationSignature
    for IntegerIdealsStructure<B>
{
    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        a * b
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> CommutativeMultiplicationSignature
    for IntegerIdealsStructure<B>
{
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> MultiplicativeMonoidSignature
    for IntegerIdealsStructure<B>
{
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> MultiplicativeAbsorptionMonoidSignature
    for IntegerIdealsStructure<B>
{
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> LeftDistributiveMultiplicationOverAddition
    for IntegerIdealsStructure<B>
{
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> RightDistributiveMultiplicationOverAddition
    for IntegerIdealsStructure<B>
{
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> SemiRingSignature
    for IntegerIdealsStructure<B>
{
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>>
    IdealsArithmeticSignature<IntegerCanonicalStructure, B> for IntegerIdealsStructure<B>
{
    fn principal_ideal(&self, a: &Integer) -> Self::Elem {
        Abs::abs(a)
    }

    fn contains_ideal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        b % a == Natural::ZERO
    }

    fn intersect(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        lcm(a.clone(), b.clone())
    }

    fn quotient(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        if b == &Natural::ZERO {
            Natural::ONE
        } else {
            a / gcd(a.clone(), b.clone())
        }
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

impl<B: BorrowedStructure<IntegerCanonicalStructure>> TryReciprocalSignature
    for IntegerIdealsStructure<B>
{
    fn try_reciprocal(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.factorization_exponents().try_reciprocal(a)
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> FavoriteAssociateSignature
    for IntegerIdealsStructure<B>
{
    fn factor_fav_assoc(&self, a: &Self::Elem) -> (Self::Elem, Self::Elem) {
        self.factorization_exponents().factor_fav_assoc(a)
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> UniqueFactorizationMonoidSignature
    for IntegerIdealsStructure<B>
{
    type FactoredExponent = NaturalCanonicalStructure;

    fn factorization_exponents(&self) -> &Self::FactoredExponent {
        Natural::structure_ref()
    }

    fn into_factorization_exponents(self) -> Self::FactoredExponent {
        Natural::structure()
    }

    fn factorization_pow(&self, a: &Self::Elem, k: &Natural) -> Self::Elem {
        self.factorization_exponents().nat_pow(a, k)
    }

    fn try_is_irreducible(&self, a: &Self::Elem) -> Option<bool> {
        self.factorization_exponents().try_is_irreducible(a)
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> FactoringMonoidSignature
    for IntegerIdealsStructure<B>
{
    fn factor_unchecked(&self, ideal: &Natural) -> Factored<Natural, Natural> {
        Natural::structure().factor_unchecked(ideal)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_ideals() {
        let ideals = Integer::structure().into_ideals();

        assert!(ideals.equal(&Natural::from(3u32), &Natural::from(3u32)));

        assert!(!ideals.equal(&Natural::from(2u32), &Natural::from(3u32)));

        assert!(ideals.contains_ideal(&Natural::from(2u32), &Natural::from(6u32)));

        assert!(!ideals.contains_ideal(&Natural::from(6u32), &Natural::from(2u32)));

        assert!(!ideals.contains_ideal(&Natural::from(5u32), &Natural::from(7u32)));

        assert!(ideals.equal(
            &ideals.add(&Natural::from(6u32), &Natural::from(15u32)),
            &Natural::from(3u32)
        ));

        assert!(
            ideals.equal(
                &Integer::structure()
                    .ideals()
                    .intersect(&Natural::from(6u32), &Natural::from(15u32)),
                &Natural::from(30u32)
            )
        );

        assert!(
            ideals.equal(
                &Integer::structure()
                    .ideals()
                    .mul(&Natural::from(6u32), &Natural::from(15u32)),
                &Natural::from(90u32)
            )
        );

        assert_eq!(ideals.generated_ideal(vec![-15, 6]), 3u32.into());
    }

    #[test]
    fn factor_integer_ideal() {
        let ideals = Integer::structure().into_ideals();

        let f = ideals.factor(
            &Integer::structure()
                .ideals()
                .principal_ideal(&Integer::from(0)),
        );
        println!("{:?}", f);
        assert!(f.is_zero());

        let f = ideals.factor(
            &Integer::structure()
                .ideals()
                .principal_ideal(&Integer::from(18)),
        );
        println!("{:?}", f);
        assert!(!f.is_zero());
    }
}
