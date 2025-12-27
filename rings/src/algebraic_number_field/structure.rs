use std::marker::PhantomData;

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
    BorrowedStructure, FiniteSetSignature, Function, InjectiveFunction, MetaType, Morphism,
    SetSignature,
};

pub trait AlgebraicIntegerRingSignature: DedekindDomainSignature + CharZeroRingSignature {
    type AlgebraicNumberField: AlgebraicNumberFieldSignature<RingOfIntegers = Self>;

    fn anf(&self) -> Self::AlgebraicNumberField;

    fn anf_inclusion<'a>(
        &'a self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<
        Self,
        &'a Self,
        Self::AlgebraicNumberField,
        Self::AlgebraicNumberField,
    > {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(self)
    }
    fn into_anf_inclusion(
        self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<
        Self,
        Self,
        Self::AlgebraicNumberField,
        Self::AlgebraicNumberField,
    > {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(self)
    }
}

/// An algebraic number field is a field of characteristic zero such that
/// the inclusion of its rational subfield is finite dimensional
pub trait AlgebraicNumberFieldSignature: CharZeroFieldSignature {
    type Basis: FiniteSetSignature;
    type RingOfIntegers: AlgebraicIntegerRingSignature<AlgebraicNumberField = Self>;
    type RationalInclusion<B: BorrowedStructure<Self>>: FiniteDimensionalFieldExtension<RationalCanonicalStructure, Self>;

    fn roi(&self) -> Self::RingOfIntegers;

    fn finite_dimensional_rational_extension<'a>(&'a self) -> Self::RationalInclusion<&'a Self>;
    fn into_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self>;

    fn ring_of_integers_extension<'a>(
        &'a self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<
        Self::RingOfIntegers,
        Self::RingOfIntegers,
        Self,
        &'a Self,
    > {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_algebraic_number_field(self)
    }
    fn into_ring_of_integers_extension(
        self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<
        Self::RingOfIntegers,
        Self::RingOfIntegers,
        Self,
        Self,
    > {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_algebraic_number_field(self)
    }

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

#[derive(Debug, Clone)]
pub struct RingOfIntegersToAlgebraicNumberFieldInclusion<
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    RB: BorrowedStructure<R>,
    K: AlgebraicNumberFieldSignature<RingOfIntegers = R>,
    KB: BorrowedStructure<K>,
> {
    _roi: PhantomData<R>,
    roi: RB,
    _anf: PhantomData<K>,
    anf: KB,
}

impl<
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    K: AlgebraicNumberFieldSignature<RingOfIntegers = R>,
    KB: BorrowedStructure<K>,
> RingOfIntegersToAlgebraicNumberFieldInclusion<R, R, K, KB>
{
    pub fn from_algebraic_number_field(anf: KB) -> Self {
        Self {
            _roi: PhantomData,
            _anf: PhantomData,
            roi: anf.borrow().roi(),
            anf,
        }
    }
}

impl<
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    RB: BorrowedStructure<R>,
    K: AlgebraicNumberFieldSignature<RingOfIntegers = R>,
> RingOfIntegersToAlgebraicNumberFieldInclusion<R, RB, K, K>
{
    pub fn from_ring_of_integers(roi: RB) -> Self {
        Self {
            _roi: PhantomData,
            _anf: PhantomData,
            anf: roi.borrow().anf(),
            roi,
        }
    }
}

impl<
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    RB: BorrowedStructure<R>,
    K: AlgebraicNumberFieldSignature<RingOfIntegers = R>,
    KB: BorrowedStructure<K>,
> Morphism<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<R, RB, K, KB>
{
    fn domain(&self) -> &R {
        self.roi.borrow()
    }

    fn range(&self) -> &K {
        self.anf.borrow()
    }
}

impl<
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    RB: BorrowedStructure<R>,
    K: AlgebraicNumberFieldSignature<RingOfIntegers = R>,
    KB: BorrowedStructure<K>,
> Function<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<R, RB, K, KB>
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
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    RB: BorrowedStructure<R>,
    K: AlgebraicNumberFieldSignature<RingOfIntegers = R>,
    KB: BorrowedStructure<K>,
> InjectiveFunction<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<R, RB, K, KB>
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
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    RB: BorrowedStructure<R>,
    K: AlgebraicNumberFieldSignature<RingOfIntegers = R>,
    KB: BorrowedStructure<K>,
> RingHomomorphism<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<R, RB, K, KB>
where
    Self: AlgebraicIntegerRingInAlgebraicNumberFieldSignature<
            AlgebraicNumberField = K,
            RingOfIntegers = R,
        >,
{
}

impl<
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    RB: BorrowedStructure<R>,
    K: AlgebraicNumberFieldSignature<RingOfIntegers = R>,
    KB: BorrowedStructure<K>,
> RingOfIntegersToAlgebraicNumberFieldInclusion<R, RB, K, KB>
{
    pub fn zq_extension<'a>(
        &'a self,
    ) -> RingOfIntegersExtension<
        R,
        &'a R,
        K,
        &'a K,
        RingOfIntegersToAlgebraicNumberFieldInclusion<R, &'a R, K, &'a K>,
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
                roi: self.domain(),
                _anf: PhantomData,
                anf: self.range(),
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
    type AlgebraicNumberField: AlgebraicNumberFieldSignature<RingOfIntegers = Self::RingOfIntegers>;
    type RingOfIntegers: AlgebraicIntegerRingSignature<
        AlgebraicNumberField = Self::AlgebraicNumberField,
    >;

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
        x: &<Self::AlgebraicNumberField as SetSignature>::Set,
    ) -> Option<<Self::RingOfIntegers as SetSignature>::Set>;
}
