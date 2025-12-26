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
    AdditiveMonoidEqSignature, FinitelyGeneratedModuleSignature, FreeModuleSignature,
    MetaAdditiveMonoidEq, MetaCharZeroRing, PrincipalRationalSubfieldInclusion, RingDivisionError,
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
    /// Represents `rational_part + algebraic_part * sqrt(d)`
    rational_part: Rational,
    algebraic_part: Rational,
}

impl QuadraticNumberFieldElement {
    const ZERO: Self = Self {
        rational_part: Rational::ZERO,
        algebraic_part: Rational::ZERO,
    };
    const ONE: Self = Self {
        rational_part: Rational::ONE,
        algebraic_part: Rational::ZERO,
    };
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
        a.rational_part == b.rational_part && a.algebraic_part == b.algebraic_part
    }
}

impl AdditiveMonoidSignature for QuadraticNumberFieldStructure {
    fn zero(&self) -> Self::Set {
        Self::Set::ZERO
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        QuadraticNumberFieldElement {
            rational_part: &a.rational_part + &b.rational_part,
            algebraic_part: &a.algebraic_part + &b.algebraic_part,
        }
    }
}

impl AdditiveGroupSignature for QuadraticNumberFieldStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        QuadraticNumberFieldElement {
            rational_part: -&a.rational_part,
            algebraic_part: -&a.algebraic_part,
        }
    }
}

impl SemiRingSignature for QuadraticNumberFieldStructure {
    fn one(&self) -> Self::Set {
        Self::Set::ONE
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        // (x + y sqrd(d))(z + w sqrt(d)) = (xz + dyw) + (xw + yz) sqrt(d)
        QuadraticNumberFieldElement {
            rational_part: &a.rational_part * &b.rational_part
                + Rational::from(&self.d) * &a.algebraic_part * &b.algebraic_part,
            algebraic_part: &a.rational_part * &b.algebraic_part
                + &a.algebraic_part * &b.rational_part,
        }
    }
}

impl RingSignature for QuadraticNumberFieldStructure {}

impl SemiRingUnitsSignature for QuadraticNumberFieldStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        // (x + y sqrt(d))^{-1} = (a - b sqrt(d)) / (x^2 + dy^2)
        debug_assert!(!self.d.is_zero()); // it's squarefree in particular non-zero
        let d = &a.rational_part * &a.rational_part
            + Rational::from(&self.d) * &a.algebraic_part * &a.algebraic_part;
        debug_assert_eq!(d == Rational::ZERO, self.is_zero(a));
        if d == Rational::ZERO {
            return Err(RingDivisionError::DivideByZero);
        }
        Ok(QuadraticNumberFieldElement {
            rational_part: &a.rational_part / &d,
            algebraic_part: -&a.algebraic_part / d,
        })
    }
}

impl IntegralDomainSignature for QuadraticNumberFieldStructure {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        Ok(self.mul(a, &self.inv(b)?))
    }
}

impl CharacteristicSignature for QuadraticNumberFieldStructure {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl CharZeroRingSignature for QuadraticNumberFieldStructure {
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer> {
        if x.algebraic_part == Rational::ZERO {
            x.rational_part.try_to_int()
        } else {
            None
        }
    }
}

impl FieldSignature for QuadraticNumberFieldStructure {}

impl CharZeroFieldSignature for QuadraticNumberFieldStructure {
    fn try_to_rat(&self, x: &Self::Set) -> Option<Rational> {
        if x.algebraic_part == Rational::ZERO {
            Some(x.rational_part.clone())
        } else {
            None
        }
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
