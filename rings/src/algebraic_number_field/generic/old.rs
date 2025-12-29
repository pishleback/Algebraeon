use crate::{
    algebraic_number_field::{AlgebraicNumberFieldSignature, RingOfIntegersExtension},
    integer::ideal::IntegerIdealsStructure,
    structure::{
        CharZeroRingSignature, DedekindDomainIdealsSignature, DedekindDomainSignature,
        RingHomomorphism, RingToIdealsSignature,
    },
};
use algebraeon_nzq::{Integer, IntegerCanonicalStructure};
use algebraeon_sets::structure::{
    BorrowedStructure, Function, InjectiveFunction, MetaType, Morphism, SetSignature,
};
use std::marker::PhantomData;

pub trait AlgebraicIntegerRingSignature: DedekindDomainSignature + CharZeroRingSignature {
    type AlgebraicNumberField: AlgebraicNumberFieldWithRingOfIntegersSignature<
        RingOfIntegers = Self,
    >;

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

pub trait AlgebraicNumberFieldWithRingOfIntegersSignature: AlgebraicNumberFieldSignature {
    type RingOfIntegers: AlgebraicIntegerRingSignature<AlgebraicNumberField = Self>;

    fn roi(&self) -> Self::RingOfIntegers;

    fn roi_inclusion<'a>(
        &'a self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<
        Self::RingOfIntegers,
        Self::RingOfIntegers,
        Self,
        &'a Self,
    > {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_algebraic_number_field(self)
    }
    fn into_roi_inclusion(
        self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<
        Self::RingOfIntegers,
        Self::RingOfIntegers,
        Self,
        Self,
    > {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_algebraic_number_field(self)
    }
}

#[derive(Debug, Clone)]
pub struct RingOfIntegersToAlgebraicNumberFieldInclusion<
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    RB: BorrowedStructure<R>,
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
> {
    _roi: PhantomData<R>,
    roi: RB,
    _anf: PhantomData<K>,
    anf: KB,
}

impl<
    R: AlgebraicIntegerRingSignature<AlgebraicNumberField = K>,
    K: AlgebraicNumberFieldWithRingOfIntegersSignature<RingOfIntegers = R>,
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
    K: AlgebraicNumberFieldSignature,
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
    K: AlgebraicNumberFieldSignature,
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
    K: AlgebraicNumberFieldSignature,
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
    K: AlgebraicNumberFieldSignature,
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
    K: AlgebraicNumberFieldSignature,
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
    K: AlgebraicNumberFieldSignature,
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
    type AlgebraicNumberField: AlgebraicNumberFieldWithRingOfIntegersSignature<
        RingOfIntegers = Self::RingOfIntegers,
    >;
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
