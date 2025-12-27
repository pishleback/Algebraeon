use crate::{
    algebraic_number_field::ring_of_integer_extensions::RingOfIntegersExtension,
    integer::ideal::IntegerIdealsStructure,
    structure::{
        CharZeroFieldSignature, CharZeroRingSignature, DedekindDomainIdealsSignature,
        DedekindDomainSignature, FiniteDimensionalFieldExtension, FiniteRankFreeRingExtension,
        MetaGreatestCommonDivisor, RingHomomorphism, RingToIdealsSignature,
    },
};
use algebraeon_nzq::{
    Integer, IntegerCanonicalStructure, Rational, RationalCanonicalStructure, traits::Fraction,
};
use algebraeon_sets::structure::{
    BorrowedStructure, FiniteSetSignature, InjectiveFunction, MetaType,
};

/// An algebraic number field is a field of characteristic zero such that
/// the inclusion of its rational subfield is finite dimensional
pub trait AlgebraicNumberFieldSignature: CharZeroFieldSignature {
    type Basis: FiniteSetSignature;
    type RingOfIntegers: DedekindDomainSignature + CharZeroRingSignature;
    type RingOfIntegersInclusion: AlgebraicIntegerRingInAlgebraicNumberField<Self>;
    type RationalInclusion<B: BorrowedStructure<Self>>: FiniteDimensionalFieldExtension<RationalCanonicalStructure, Self>;

    fn finite_dimensional_rational_extension<'a>(&'a self) -> Self::RationalInclusion<&'a Self>;
    fn into_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self>;

    fn into_ring_of_integers_extension(self) -> Self::RingOfIntegersInclusion;

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

    fn absolute_degree(&self) -> usize {
        self.finite_dimensional_rational_extension().degree()
    }

    /// return a scalar multiple of $a$ which is an algebraic integer
    /// need not return $a$ itself when $a$ is already an algebraic integer
    fn integral_multiple(&self, a: &Self::Set) -> Self::Set {
        let m = self.min_poly_denominator_lcm(a);
        let b = self.mul(&self.try_from_rat(&Rational::from(m)).unwrap(), a);
        debug_assert!(self.is_algebraic_integer(&b));
        b
    }
}

pub trait AlgebraicIntegerRingInAlgebraicNumberField<
    ANF: AlgebraicNumberFieldSignature<RingOfIntegersInclusion = Self>,
>:
    RingHomomorphism<ANF::RingOfIntegers, ANF> + InjectiveFunction<ANF::RingOfIntegers, ANF>
{
    fn discriminant(&self) -> Integer;

    fn zq_extension<'a>(
        &'a self,
    ) -> RingOfIntegersExtension<
        ANF,
        &'a Self,
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        &'a ANF::RingOfIntegers,
        <ANF::RingOfIntegers as RingToIdealsSignature>::Ideals<&'a ANF::RingOfIntegers>,
    >
    where
        ANF::RingOfIntegers: RingToIdealsSignature,
        <ANF::RingOfIntegers as RingToIdealsSignature>::Ideals<&'a ANF::RingOfIntegers>:
            DedekindDomainIdealsSignature<ANF::RingOfIntegers, &'a ANF::RingOfIntegers>,
    {
        let ideals_z = Integer::structure().into_ideals();
        let ideals_r = self.domain().ideals();
        RingOfIntegersExtension::new_integer_extension(self, ideals_z, ideals_r)
    }
}
