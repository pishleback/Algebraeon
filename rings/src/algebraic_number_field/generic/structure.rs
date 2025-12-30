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

    fn n(&self) -> usize {
        self.inbound_finite_dimensional_rational_extension()
            .degree()
    }

    fn generator(&self) -> Self::Set;

    fn is_algebraic_integer(&self, a: &Self::Set) -> bool;

    fn maximal_order<'a>(&'a self) -> OrderWithBasis<Self, &'a Self, true>;
    fn into_maximal_order(self) -> OrderWithBasis<Self, Self, true>;

    fn discriminant(&self) -> Integer;

    // This is the LCM of the denominators of the coefficients of the minimal polynomial of a,
    // and thus it may well be >1 even when the element a is an algebraic integer.
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

    /// return a scalar multiple of $a$ which is an algebraic integer
    /// need not return $a$ itself when $a$ is already an algebraic integer
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
