use crate::structure::{
    CharZeroFieldSignature, CharZeroRingSignature, DedekindDomainSignature,
    FiniteDimensionalFieldExtension, MetaGreatestCommonDivisor, RingHomomorphism,
};
use algebraeon_nzq::{Integer, Rational, RationalCanonicalStructure, traits::Fraction};
use algebraeon_sets::structure::{BorrowedStructure, InjectiveFunction};

/// An algebraic number field is a field of characteristic zero such that
/// the inclusion of its rational subfield is finite dimensional
pub trait AlgebraicNumberFieldSignature: CharZeroFieldSignature {
    type RingOfIntegers: DedekindDomainSignature + CharZeroRingSignature;
    type RingOfIntegersInclusion<B: BorrowedStructure<Self>>: AlgebraicIntegerRingInAlgebraicNumberField<Self>;

    fn finite_dimensional_rational_extension<'a>(
        &'a self,
    ) -> impl FiniteDimensionalFieldExtension<RationalCanonicalStructure, Self>;
    fn into_finite_dimensional_rational_extension(
        self,
    ) -> impl FiniteDimensionalFieldExtension<RationalCanonicalStructure, Self>;

    fn ring_of_integer_extension<'a>(&'a self) -> Self::RingOfIntegersInclusion<&'a Self>;
    fn into_ring_of_integer_extension(self) -> Self::RingOfIntegersInclusion<Self>;

    fn compute_integral_basis_and_discriminant(&self) -> (Vec<Self::Set>, Integer);

    fn is_algebraic_integer(&self, a: &Self::Set) -> bool;

    // This is the LCM of the denominators of the coefficients of the minimal polynomial of a,
    // and thus it may well be >1 even when the element a is an algebraic integer.
    fn min_poly_denominator_lcm(&self, a: &Self::Set) -> Integer {
        Integer::lcm_list(
            self.finite_dimensional_rational_extension()
                .min_poly(a)
                .coeffs()
                .into_iter()
                .map(|c| Integer::from(c.denominator()))
                .collect(),
        )
    }

    /// return a scalar multiple of $a$ which is an algebraic integer
    /// need not return $a$ itself when $a$ is already an algebraic integer
    fn integral_multiple(&self, a: &Self::Set) -> Self::Set {
        let m = self.min_poly_denominator_lcm(a);
        let b = self.mul(&self.from_rat(&Rational::from(m)).unwrap(), a);
        debug_assert!(self.is_algebraic_integer(&b));
        b
    }
}

pub trait AlgebraicIntegerRingInAlgebraicNumberField<ANF: AlgebraicNumberFieldSignature>:
    RingHomomorphism<ANF::RingOfIntegers, ANF> + InjectiveFunction<ANF::RingOfIntegers, ANF>
{
}
