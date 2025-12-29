use crate::structure::{
    CharZeroFieldSignature, FiniteDimensionalFieldExtension, FiniteRankFreeRingExtension,
    MetaGreatestCommonDivisor,
};
use algebraeon_nzq::{Integer, Rational, RationalCanonicalStructure, traits::Fraction};
use algebraeon_sets::structure::{BorrowedStructure, FiniteSetSignature};

/// An algebraic number field is a field of characteristic zero such that
/// the inclusion of its rational subfield is finite dimensional
pub trait AlgebraicNumberFieldSignature: CharZeroFieldSignature {
    type Basis: FiniteSetSignature;
    type RationalInclusion<B: BorrowedStructure<Self>>: FiniteDimensionalFieldExtension<RationalCanonicalStructure, Self>;

    fn finite_dimensional_rational_extension<'a>(&'a self) -> Self::RationalInclusion<&'a Self>;
    fn into_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self>;

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

    fn n(&self) -> usize {
        self.finite_dimensional_rational_extension().degree()
    }

    fn is_algebraic_integer(&self, a: &Self::Set) -> bool;

    /// return a scalar multiple of $a$ which is an algebraic integer
    /// need not return $a$ itself when $a$ is already an algebraic integer
    fn integral_multiple(&self, a: &Self::Set) -> Self::Set {
        let m = self.min_poly_denominator_lcm(a);
        let b = self.mul(&self.try_from_rat(&Rational::from(m)).unwrap(), a);
        debug_assert!(self.is_algebraic_integer(&b));
        b
    }
}
