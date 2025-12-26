use crate::{
    algebraic_number_field::quadratic_number_field::QuadraticNumberFieldElement,
    structure::{
        AdditiveGroupSignature, AdditiveMonoidSignature, CharZeroRingSignature,
        CharacteristicSignature, DedekindDomainSignature, IntegralDomainSignature, RingSignature,
        SemiRingSignature, SemiRingUnitsSignature,
    },
};
use algebraeon_nzq::{Integer, Natural};
use algebraeon_sets::structure::{EqSignature, SetSignature, Signature};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QuadraticRingOfIntegersStructure {
    /// A squarefree integer
    d: Integer,
}

impl Signature for QuadraticRingOfIntegersStructure {}

impl SetSignature for QuadraticRingOfIntegersStructure {
    type Set = QuadraticNumberFieldElement;

    fn is_element(&self, _: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl EqSignature for QuadraticRingOfIntegersStructure {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        todo!()
    }
}

impl AdditiveMonoidSignature for QuadraticRingOfIntegersStructure {
    fn zero(&self) -> Self::Set {
        todo!()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        todo!()
    }
}

impl AdditiveGroupSignature for QuadraticRingOfIntegersStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        todo!()
    }
}

impl SemiRingSignature for QuadraticRingOfIntegersStructure {
    fn one(&self) -> Self::Set {
        todo!()
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        todo!()
    }
}

impl RingSignature for QuadraticRingOfIntegersStructure {}

impl SemiRingUnitsSignature for QuadraticRingOfIntegersStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, crate::structure::RingDivisionError> {
        todo!()
    }
}

impl IntegralDomainSignature for QuadraticRingOfIntegersStructure {
    fn div(
        &self,
        a: &Self::Set,
        b: &Self::Set,
    ) -> Result<Self::Set, crate::structure::RingDivisionError> {
        todo!()
    }
}

impl CharacteristicSignature for QuadraticRingOfIntegersStructure {
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl CharZeroRingSignature for QuadraticRingOfIntegersStructure {
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer> {
        todo!()
    }
}

impl DedekindDomainSignature for QuadraticRingOfIntegersStructure {}
