use algebraeon_nzq::{Integer, Natural, Rational, RationalCanonicalStructure};
use algebraeon_sets::structure::{
    BorrowedStructure, Function, InjectiveFunction, MetaType, Morphism,
};
use algebraeon_sets::structure::{
    CanonicalStructure, CountableSetSignature, EqSignature, FiniteSetSignature, SetSignature,
    Signature,
};

use crate::algebraic_number_field::quadratic_ring_of_integers::QuadraticRingOfIntegersStructure;
use crate::algebraic_number_field::structure::AlgebraicIntegerRingInAlgebraicNumberField;
use crate::structure::{
    FinitelyGeneratedModuleSignature, FreeModuleSignature, PrincipalRationalSubfieldInclusion,
    RingHomomorphism, RingHomomorphismRangeModuleStructure, SemiModuleSignature,
};
use crate::{
    algebraic_number_field::structure::AlgebraicNumberFieldSignature,
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, CharZeroFieldSignature,
        CharZeroRingSignature, CharacteristicSignature, FieldSignature, IntegralDomainSignature,
        RingSignature, SemiRingSignature, SemiRingUnitsSignature,
    },
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, CanonicalStructure)]
#[canonical_structure(eq)]
pub enum QuadraticNumberFieldBasis {
    Real,
    Algebraic,
}

impl CountableSetSignature for QuadraticNumberFieldBasisCanonicalStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
        vec![
            QuadraticNumberFieldBasis::Real,
            QuadraticNumberFieldBasis::Algebraic,
        ]
        .into_iter()
    }
}

impl FiniteSetSignature for QuadraticNumberFieldBasisCanonicalStructure {
    fn size(&self) -> usize {
        2
    }
}

#[derive(Debug, Clone)]
pub struct QuadraticNumberFieldElement {
    integer_part: Rational,
    algebraic_part: Rational,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuadraticNumberFieldStructure {
    /// A squarefree integer
    d: Integer,
}

impl Signature for QuadraticNumberFieldStructure {}

impl SetSignature for QuadraticNumberFieldStructure {
    type Set = QuadraticNumberFieldElement;

    fn is_element(&self, _: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl EqSignature for QuadraticNumberFieldStructure {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        todo!()
    }
}

impl AdditiveMonoidSignature for QuadraticNumberFieldStructure {
    fn zero(&self) -> Self::Set {
        todo!()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        todo!()
    }
}

impl AdditiveGroupSignature for QuadraticNumberFieldStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        todo!()
    }
}

impl SemiRingSignature for QuadraticNumberFieldStructure {
    fn one(&self) -> Self::Set {
        todo!()
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        todo!()
    }
}

impl RingSignature for QuadraticNumberFieldStructure {}

impl SemiRingUnitsSignature for QuadraticNumberFieldStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, crate::structure::RingDivisionError> {
        todo!()
    }
}

impl IntegralDomainSignature for QuadraticNumberFieldStructure {
    fn div(
        &self,
        a: &Self::Set,
        b: &Self::Set,
    ) -> Result<Self::Set, crate::structure::RingDivisionError> {
        todo!()
    }
}

impl CharacteristicSignature for QuadraticNumberFieldStructure {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl CharZeroRingSignature for QuadraticNumberFieldStructure {
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer> {
        todo!()
    }
}

impl FieldSignature for QuadraticNumberFieldStructure {}

impl CharZeroFieldSignature for QuadraticNumberFieldStructure {
    fn try_to_rat(&self, x: &Self::Set) -> Option<Rational> {
        todo!()
    }
}

impl SemiModuleSignature<RationalCanonicalStructure> for QuadraticNumberFieldStructure {
    fn ring(&self) -> &RationalCanonicalStructure {
        todo!()
    }

    fn scalar_mul(&self, a: &Self::Set, x: &Rational) -> Self::Set {
        todo!()
    }
}

impl FreeModuleSignature<RationalCanonicalStructure> for QuadraticNumberFieldStructure {
    type Basis = QuadraticNumberFieldBasisCanonicalStructure;

    fn basis_set(&self) -> impl std::borrow::Borrow<Self::Basis> {
        QuadraticNumberFieldBasis::structure()
    }

    fn to_component<'a>(
        &self,
        b: &QuadraticNumberFieldBasis,
        v: &'a QuadraticNumberFieldElement,
    ) -> std::borrow::Cow<'a, Rational> {
        todo!()
    }

    fn from_component(&self, b: &QuadraticNumberFieldBasis, r: &Rational) -> Self::Set {
        todo!()
    }
}

static_assertions::const_assert!(
    impls::impls!(QuadraticNumberFieldBasisCanonicalStructure : FiniteSetSignature)
);

static_assertions::const_assert!(
    impls::impls!(QuadraticNumberFieldStructure : FinitelyGeneratedModuleSignature<RationalCanonicalStructure>)
);

static_assertions::const_assert!(
    impls::impls!(<QuadraticNumberFieldStructure as FreeModuleSignature<RationalCanonicalStructure>>::Basis : FiniteSetSignature)
);

impl<'h, B: BorrowedStructure<QuadraticNumberFieldStructure>>
    FreeModuleSignature<RationalCanonicalStructure>
    for RingHomomorphismRangeModuleStructure<
        'h,
        RationalCanonicalStructure,
        QuadraticNumberFieldStructure,
        PrincipalRationalSubfieldInclusion<QuadraticNumberFieldStructure, B>,
    >
{
    type Basis = QuadraticNumberFieldBasisCanonicalStructure;

    fn basis_set(&self) -> impl std::borrow::Borrow<Self::Basis> {
        QuadraticNumberFieldBasis::structure()
    }

    fn to_component<'a>(
        &self,
        b: &QuadraticNumberFieldBasis,
        v: &'a QuadraticNumberFieldElement,
    ) -> std::borrow::Cow<'a, Rational> {
        todo!()
    }

    fn from_component(&self, b: &QuadraticNumberFieldBasis, r: &Rational) -> Self::Set {
        todo!()
    }
}

#[derive(Debug, Clone)]
pub struct QuadraticNumberFieldRingOfIntegersInclusionStructure {}

impl Morphism<QuadraticRingOfIntegersStructure, QuadraticNumberFieldStructure>
    for QuadraticNumberFieldRingOfIntegersInclusionStructure
{
    fn domain(&self) -> &QuadraticRingOfIntegersStructure {
        todo!()
    }

    fn range(&self) -> &QuadraticNumberFieldStructure {
        todo!()
    }
}

impl Function<QuadraticRingOfIntegersStructure, QuadraticNumberFieldStructure>
    for QuadraticNumberFieldRingOfIntegersInclusionStructure
{
    fn image(
        &self,
        x: &<QuadraticRingOfIntegersStructure as SetSignature>::Set,
    ) -> <QuadraticNumberFieldStructure as SetSignature>::Set {
        todo!()
    }
}

impl InjectiveFunction<QuadraticRingOfIntegersStructure, QuadraticNumberFieldStructure>
    for QuadraticNumberFieldRingOfIntegersInclusionStructure
{
    fn try_preimage(
        &self,
        y: &<QuadraticNumberFieldStructure as SetSignature>::Set,
    ) -> Option<<QuadraticRingOfIntegersStructure as SetSignature>::Set> {
        todo!()
    }
}

impl RingHomomorphism<QuadraticRingOfIntegersStructure, QuadraticNumberFieldStructure>
    for QuadraticNumberFieldRingOfIntegersInclusionStructure
{
}

impl AlgebraicIntegerRingInAlgebraicNumberField<QuadraticNumberFieldStructure>
    for QuadraticNumberFieldRingOfIntegersInclusionStructure
{
    fn discriminant(&self) -> Integer {
        todo!()
    }
}

impl AlgebraicNumberFieldSignature for QuadraticNumberFieldStructure {
    type Basis = QuadraticNumberFieldBasisCanonicalStructure;
    type RingOfIntegers = QuadraticRingOfIntegersStructure;
    type RingOfIntegersInclusion = QuadraticNumberFieldRingOfIntegersInclusionStructure;
    type RationalInclusion<B: BorrowedStructure<Self>> =
        PrincipalRationalSubfieldInclusion<Self, B>;

    fn finite_dimensional_rational_extension<'a>(&'a self) -> Self::RationalInclusion<&'a Self> {
        todo!()
    }

    fn into_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self> {
        todo!()
    }

    fn into_ring_of_integers_extension(self) -> Self::RingOfIntegersInclusion {
        todo!()
    }

    fn is_algebraic_integer(&self, a: &Self::Set) -> bool {
        todo!()
    }
}
