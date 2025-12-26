use crate::algebraic_number_field::quadratic_ring_of_integers::QuadraticRingOfIntegersStructure;
use crate::algebraic_number_field::structure::AlgebraicIntegerRingInAlgebraicNumberField;
use crate::structure::{
    AdditiveMonoidEqSignature, FreeModuleSignature, MetaAdditiveMonoidEq, MetaCharZeroRing,
    MetaFactorableSignature, PrincipalRationalSubfieldInclusion, RingDivisionError,
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
use algebraeon_nzq::{Integer, Natural, Rational, RationalCanonicalStructure};
use algebraeon_sets::structure::{
    BorrowedStructure, Function, InjectiveFunction, MetaType, Morphism,
};
use algebraeon_sets::structure::{
    CanonicalStructure, CountableSetSignature, EqSignature, FiniteSetSignature, SetSignature,
    Signature,
};
use std::borrow::Cow;

#[derive(Debug, Clone, Copy, PartialEq, Eq, CanonicalStructure)]
#[canonical_structure(eq)]
pub enum QuadraticNumberFieldBasis {
    Rational,
    Algebraic,
}

impl CountableSetSignature for QuadraticNumberFieldBasisCanonicalStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
        vec![
            QuadraticNumberFieldBasis::Rational,
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
pub struct QuadraticNumberFieldStructure<D: BorrowedStructure<Integer>> {
    /// A squarefree integer
    d: D,
}

impl QuadraticNumberFieldStructure<Integer> {
    pub fn new_unchecked(d: Integer) -> Self {
        debug_assert!(d.is_squarefree());
        Self { d }
    }
}

impl<'a> QuadraticNumberFieldStructure<&'a Integer> {
    pub fn new_unchecked_ref(d: &'a Integer) -> Self {
        debug_assert!(d.is_squarefree());
        Self { d }
    }
}

impl<D: BorrowedStructure<Integer>> QuadraticNumberFieldStructure<D> {
    pub fn d(&self) -> &Integer {
        self.d.borrow()
    }

    pub fn roi<'d>(&'d self) -> QuadraticRingOfIntegersStructure<&'d Integer> {
        QuadraticRingOfIntegersStructure::new_unchecked_ref(self.d())
    }
}

impl<D: BorrowedStructure<Integer>> Signature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedStructure<Integer>> SetSignature for QuadraticNumberFieldStructure<D> {
    type Set = QuadraticNumberFieldElement;

    fn is_element(&self, _: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl<D: BorrowedStructure<Integer>> EqSignature for QuadraticNumberFieldStructure<D> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        a.rational_part == b.rational_part && a.algebraic_part == b.algebraic_part
    }
}

impl<D: BorrowedStructure<Integer>> AdditiveMonoidSignature for QuadraticNumberFieldStructure<D> {
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

impl<D: BorrowedStructure<Integer>> AdditiveGroupSignature for QuadraticNumberFieldStructure<D> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        QuadraticNumberFieldElement {
            rational_part: -&a.rational_part,
            algebraic_part: -&a.algebraic_part,
        }
    }
}

impl<D: BorrowedStructure<Integer>> SemiRingSignature for QuadraticNumberFieldStructure<D> {
    fn one(&self) -> Self::Set {
        Self::Set::ONE
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        // (x + y sqrd(d))(z + w sqrt(d)) = (xz + dyw) + (xw + yz) sqrt(d)
        QuadraticNumberFieldElement {
            rational_part: &a.rational_part * &b.rational_part
                + Rational::from(self.d()) * &a.algebraic_part * &b.algebraic_part,
            algebraic_part: &a.rational_part * &b.algebraic_part
                + &a.algebraic_part * &b.rational_part,
        }
    }
}

impl<D: BorrowedStructure<Integer>> RingSignature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedStructure<Integer>> SemiRingUnitsSignature for QuadraticNumberFieldStructure<D> {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        // (x + y sqrt(d))^{-1} = (a - b sqrt(d)) / (x^2 + dy^2)
        debug_assert!(!self.d().is_zero()); // it's squarefree in particular non-zero
        let d = &a.rational_part * &a.rational_part
            + Rational::from(self.d()) * &a.algebraic_part * &a.algebraic_part;
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

impl<D: BorrowedStructure<Integer>> IntegralDomainSignature for QuadraticNumberFieldStructure<D> {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        Ok(self.mul(a, &self.inv(b)?))
    }
}

impl<D: BorrowedStructure<Integer>> CharacteristicSignature for QuadraticNumberFieldStructure<D> {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl<D: BorrowedStructure<Integer>> CharZeroRingSignature for QuadraticNumberFieldStructure<D> {
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer> {
        if x.algebraic_part == Rational::ZERO {
            x.rational_part.try_to_int()
        } else {
            None
        }
    }
}

impl<D: BorrowedStructure<Integer>> FieldSignature for QuadraticNumberFieldStructure<D> {}

impl<D: BorrowedStructure<Integer>> CharZeroFieldSignature for QuadraticNumberFieldStructure<D> {
    fn try_to_rat(&self, x: &Self::Set) -> Option<Rational> {
        if x.algebraic_part == Rational::ZERO {
            Some(x.rational_part.clone())
        } else {
            None
        }
    }
}

impl<D: BorrowedStructure<Integer>> SemiModuleSignature<RationalCanonicalStructure>
    for QuadraticNumberFieldStructure<D>
{
    fn ring(&self) -> &RationalCanonicalStructure {
        Rational::structure_ref()
    }

    fn scalar_mul(&self, a: &Self::Set, x: &Rational) -> Self::Set {
        QuadraticNumberFieldElement {
            rational_part: x * &a.rational_part,
            algebraic_part: x * &a.algebraic_part,
        }
    }
}

impl<'h, D: BorrowedStructure<Integer>, B: BorrowedStructure<QuadraticNumberFieldStructure<D>>>
    FreeModuleSignature<RationalCanonicalStructure>
    for RingHomomorphismRangeModuleStructure<
        'h,
        RationalCanonicalStructure,
        QuadraticNumberFieldStructure<D>,
        PrincipalRationalSubfieldInclusion<QuadraticNumberFieldStructure<D>, B>,
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
        match b {
            QuadraticNumberFieldBasis::Rational => Cow::Borrowed(&v.rational_part),
            QuadraticNumberFieldBasis::Algebraic => Cow::Borrowed(&v.algebraic_part),
        }
    }

    fn from_component(&self, b: &QuadraticNumberFieldBasis, r: &Rational) -> Self::Set {
        match b {
            QuadraticNumberFieldBasis::Rational => QuadraticNumberFieldElement {
                rational_part: r.clone(),
                algebraic_part: Rational::ZERO,
            },
            QuadraticNumberFieldBasis::Algebraic => QuadraticNumberFieldElement {
                rational_part: Rational::ZERO,
                algebraic_part: r.clone(),
            },
        }
    }
}

#[derive(Debug, Clone)]
pub struct QuadraticNumberFieldRingOfIntegersInclusionStructure<D: BorrowedStructure<Integer>> {
    qanf: QuadraticNumberFieldStructure<D>,
    qroi: QuadraticRingOfIntegersStructure<D>,
}

impl QuadraticNumberFieldRingOfIntegersInclusionStructure<Integer> {
    pub fn new_unchecked(&self, d: Integer) -> Self {
        debug_assert!(d.is_squarefree());
        Self {
            qanf: QuadraticNumberFieldStructure::new_unchecked(d.clone()),
            qroi: QuadraticRingOfIntegersStructure::new_unchecked(d),
        }
    }
}

impl<'d> QuadraticNumberFieldRingOfIntegersInclusionStructure<&'d Integer> {
    pub fn new_unchecked(&self, d: &'d Integer) -> Self {
        debug_assert!(d.is_squarefree());
        Self {
            qanf: QuadraticNumberFieldStructure::new_unchecked_ref(d),
            qroi: QuadraticRingOfIntegersStructure::new_unchecked_ref(d),
        }
    }
}

impl<D: BorrowedStructure<Integer>>
    Morphism<QuadraticRingOfIntegersStructure<D>, QuadraticNumberFieldStructure<D>>
    for QuadraticNumberFieldRingOfIntegersInclusionStructure<D>
{
    fn domain(&self) -> &QuadraticRingOfIntegersStructure<D> {
        &self.qroi
    }

    fn range(&self) -> &QuadraticNumberFieldStructure<D> {
        &self.qanf
    }
}

impl<D: BorrowedStructure<Integer>>
    Function<QuadraticRingOfIntegersStructure<D>, QuadraticNumberFieldStructure<D>>
    for QuadraticNumberFieldRingOfIntegersInclusionStructure<D>
{
    fn image(&self, x: &QuadraticNumberFieldElement) -> QuadraticNumberFieldElement {
        x.clone()
    }
}

impl<D: BorrowedStructure<Integer>>
    InjectiveFunction<QuadraticRingOfIntegersStructure<D>, QuadraticNumberFieldStructure<D>>
    for QuadraticNumberFieldRingOfIntegersInclusionStructure<D>
{
    fn try_preimage(&self, y: &QuadraticNumberFieldElement) -> Option<QuadraticNumberFieldElement> {
        todo!()
    }
}

impl<D: BorrowedStructure<Integer>>
    RingHomomorphism<QuadraticRingOfIntegersStructure<D>, QuadraticNumberFieldStructure<D>>
    for QuadraticNumberFieldRingOfIntegersInclusionStructure<D>
{
}

impl<D: BorrowedStructure<Integer>>
    AlgebraicIntegerRingInAlgebraicNumberField<QuadraticNumberFieldStructure<D>>
    for QuadraticNumberFieldRingOfIntegersInclusionStructure<D>
{
    fn discriminant(&self) -> Integer {
        todo!()
    }
}

impl<D: BorrowedStructure<Integer>> AlgebraicNumberFieldSignature
    for QuadraticNumberFieldStructure<D>
{
    type Basis = QuadraticNumberFieldBasisCanonicalStructure;
    type RingOfIntegers = QuadraticRingOfIntegersStructure<D>;
    type RingOfIntegersInclusion = QuadraticNumberFieldRingOfIntegersInclusionStructure<D>;
    type RationalInclusion<B: BorrowedStructure<Self>> =
        PrincipalRationalSubfieldInclusion<Self, B>;

    fn finite_dimensional_rational_extension<'a>(&'a self) -> Self::RationalInclusion<&'a Self> {
        PrincipalRationalSubfieldInclusion::new(self)
    }

    fn into_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self> {
        PrincipalRationalSubfieldInclusion::new(self)
    }

    fn into_ring_of_integers_extension(self) -> Self::RingOfIntegersInclusion {
        todo!()
    }

    fn is_algebraic_integer(&self, a: &Self::Set) -> bool {
        todo!()
    }
}
