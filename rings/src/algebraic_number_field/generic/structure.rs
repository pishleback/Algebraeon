use std::marker::PhantomData;

use crate::{
    algebraic_number_field::{OrderWithBasis, RingOfIntegersExtension},
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
    BorrowedStructure, FiniteSetSignature, Function, InjectiveFunction, MetaType, Morphism,
    SetSignature,
};

/// An algebraic number field is a field of characteristic zero such that
/// the inclusion of its rational subfield is finite dimensional
pub trait AlgebraicNumberFieldSignature: CharZeroFieldSignature {
    type Basis: FiniteSetSignature;
    type RationalInclusion<B: BorrowedStructure<Self>>: FiniteDimensionalFieldExtension<RationalCanonicalStructure, Self>;

    fn inbound_finite_dimensional_rational_extension<'a>(
        &'a self,
    ) -> Self::RationalInclusion<&'a Self>;
    fn into_inbound_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self>;

    /// The dimension of this algebraic number field as a vector space over the rational numbers
    fn n(&self) -> usize {
        self.inbound_finite_dimensional_rational_extension()
            .degree()
    }

    /// An element which generates this algebraic number field when adjoined to the rational numbers
    /// Such an element always exists by the primitive element theorem
    fn generator(&self) -> Self::Set;

    /// Determine whether an element is integral over the integers i.e. is it a root of a monic integer polynomial
    fn is_algebraic_integer(&self, a: &Self::Set) -> bool;

    /// The discriminant of this algebraic number field i.e. the discriminant of its ring of integers
    /// Implementations should not compute this by constructing the ring of integers, as the constructor for a maximal OrderWithBasis calls this function to validate its input
    fn discriminant(&self) -> Integer;

    /// A list of self.n() elements which generate the ring of integers as a Z-module
    fn integral_basis(&self) -> Vec<Self::Set>;

    fn ring_of_integers<'a>(&'a self) -> OrderWithBasis<Self, &'a Self, true> {
        OrderWithBasis::new_maximal_unchecked(self, self.integral_basis())
    }
    fn into_ring_of_integers(self) -> OrderWithBasis<Self, Self, true> {
        let basis = self.integral_basis();
        OrderWithBasis::new_maximal_unchecked(self, basis)
    }

    /// The LCM of the denominators of the coefficients of the minimal polynomial of a.
    ///
    /// It may well be >1 even when the element a is an algebraic integer.
    fn min_poly_denominator_lcm(&self, a: &Self::Set) -> Integer {
        Integer::lcm_list(
            self.inbound_finite_dimensional_rational_extension()
                .min_poly(a)
                .coeffs()
                .into_iter()
                .map(|c| Integer::from(c.denominator()))
                .collect(),
        )
    }

    /// A scalar multiple of $a$ which is an algebraic integer.
    ///
    /// It need not return $a$ itself when $a$ is already an algebraic integer.
    fn integral_multiple(&self, a: &Self::Set) -> Self::Set {
        let m = self.min_poly_denominator_lcm(a);
        let b = self.mul(&self.try_from_rat(&Rational::from(m)).unwrap(), a);
        debug_assert!(self.is_algebraic_integer(&b));
        b
    }
}

pub trait AlgebraicIntegerRingSignature<K: AlgebraicNumberFieldSignature>:
    DedekindDomainSignature + CharZeroRingSignature
{
    fn anf(&self) -> &K;

    fn into_outbound_roi_to_anf_inclusion(
        self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<K, Self, Self> {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(self)
    }

    fn outbound_roi_to_anf_inclusion<'a>(
        &'a self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<K, Self, &'a Self> {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(self)
    }
}

#[derive(Debug, Clone)]
pub struct RingOfIntegersToAlgebraicNumberFieldInclusion<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> {
    _roi: PhantomData<R>,
    roi: RB,
    _anf: PhantomData<K>,
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
{
    pub fn from_ring_of_integers(roi: RB) -> Self {
        Self {
            _roi: PhantomData,
            _anf: PhantomData,
            roi,
        }
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> Morphism<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
{
    fn domain(&self) -> &R {
        self.roi.borrow()
    }

    fn range(&self) -> &K {
        self.roi.borrow().anf()
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> Function<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
where
    Self: AlgebraicIntegerRingInAlgebraicNumberFieldSignature<
            AlgebraicNumberField = K,
            RingOfIntegers = R,
        >,
{
    fn image(&self, x: &<R as SetSignature>::Set) -> <K as SetSignature>::Set {
        self.roi_to_anf(x)
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> InjectiveFunction<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
where
    Self: AlgebraicIntegerRingInAlgebraicNumberFieldSignature<
            AlgebraicNumberField = K,
            RingOfIntegers = R,
        >,
{
    fn try_preimage(&self, x: &<K as SetSignature>::Set) -> Option<<R as SetSignature>::Set> {
        self.try_anf_to_roi(x)
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> RingHomomorphism<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
where
    Self: AlgebraicIntegerRingInAlgebraicNumberFieldSignature<
            AlgebraicNumberField = K,
            RingOfIntegers = R,
        >,
{
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
{
    pub fn zq_extension<'a>(
        &'a self,
    ) -> RingOfIntegersExtension<
        K,
        R,
        &'a R,
        RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, &'a R>,
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        <R as RingToIdealsSignature>::Ideals<&'a R>,
    >
    where
        R: RingToIdealsSignature,
        <R as RingToIdealsSignature>::Ideals<&'a R>: DedekindDomainIdealsSignature<R, &'a R>,
    {
        let ideals_z = Integer::structure().into_ideals();
        let ideals_r = self.domain().ideals();
        RingOfIntegersExtension::new_integer_extension(
            RingOfIntegersToAlgebraicNumberFieldInclusion {
                _roi: PhantomData,
                _anf: PhantomData,
                roi: self.domain(),
            },
            ideals_z,
            ideals_r,
        )
    }
}

pub trait AlgebraicIntegerRingInAlgebraicNumberFieldSignature:
    RingHomomorphism<Self::RingOfIntegers, Self::AlgebraicNumberField>
    + InjectiveFunction<Self::RingOfIntegers, Self::AlgebraicNumberField>
{
    type AlgebraicNumberField: AlgebraicNumberFieldSignature;
    type RingOfIntegers: AlgebraicIntegerRingSignature<Self::AlgebraicNumberField>;

    fn roi(&self) -> &Self::RingOfIntegers {
        self.domain()
    }

    fn anf(&self) -> &Self::AlgebraicNumberField {
        self.range()
    }

    fn discriminant(&self) -> Integer;

    fn roi_to_anf(
        &self,
        x: &<Self::RingOfIntegers as SetSignature>::Set,
    ) -> <Self::AlgebraicNumberField as SetSignature>::Set;

    fn try_anf_to_roi(
        &self,
        y: &<Self::AlgebraicNumberField as SetSignature>::Set,
    ) -> Option<<Self::RingOfIntegers as SetSignature>::Set>;
}
