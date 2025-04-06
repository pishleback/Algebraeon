use super::number_field::AlgebraicNumberFieldStructure;
use crate::{
    polynomial::Polynomial,
    structure::{
        FiniteUnitsStructure, IntegralDomainStructure, RingStructure, SemiRingStructure,
        UnitsStructure,
    },
};
use algebraeon_nzq::{Integer, Rational};
use algebraeon_sets::structure::{EqStructure, PartialEqStructure, Structure};

#[derive(Debug, Clone)]
pub struct RingOfIntegersWithIntegralBasisStructure {
    algebraic_number_field: AlgebraicNumberFieldStructure,
    integral_basis: Vec<Polynomial<Rational>>,
    discriminant: Integer,
}

impl PartialEq for RingOfIntegersWithIntegralBasisStructure {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self, other)
    }
}
impl Eq for RingOfIntegersWithIntegralBasisStructure {}

impl RingOfIntegersWithIntegralBasisStructure {
    pub fn new(
        algebraic_number_field: AlgebraicNumberFieldStructure,
        integral_basis: Vec<Polynomial<Rational>>,
        discriminant: Integer,
    ) -> Self {
        debug_assert_eq!(
            algebraic_number_field.discriminant(&integral_basis),
            discriminant
        );
        let (_, true_discriminant) =
            algebraic_number_field.compute_integral_basis_and_discriminant();
        debug_assert_eq!(discriminant, true_discriminant);
        for a in &integral_basis {
            debug_assert!(algebraic_number_field.is_algebraic_integer(a))
        }
        Self {
            algebraic_number_field,
            integral_basis,
            discriminant,
        }
    }

    pub fn degree(&self) -> usize {
        self.algebraic_number_field.degree()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingOfIntegersWithIntegralBasisElement {
    coefficients: Vec<Integer>,
}

impl RingOfIntegersWithIntegralBasisStructure {
    pub fn roi_to_anf(&self, elem: RingOfIntegersWithIntegralBasisElement) -> Polynomial<Rational> {
        todo!()
    }

    pub fn anf_to_roi(
        &self,
        elem: Polynomial<Rational>,
    ) -> Option<RingOfIntegersWithIntegralBasisElement> {
        todo!()
    }
}

impl Structure for RingOfIntegersWithIntegralBasisStructure {
    type Set = RingOfIntegersWithIntegralBasisElement;
}

impl PartialEqStructure for RingOfIntegersWithIntegralBasisStructure {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        a == b
    }
}

impl EqStructure for RingOfIntegersWithIntegralBasisStructure {}

impl SemiRingStructure for RingOfIntegersWithIntegralBasisStructure {
    fn zero(&self) -> Self::Set {
        let coefficients = vec![Integer::ZERO; self.degree()];
        RingOfIntegersWithIntegralBasisElement { coefficients }
    }

    fn one(&self) -> Self::Set {
        todo!()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let n = self.degree();
        debug_assert_eq!(a.coefficients.len(), n);
        debug_assert_eq!(b.coefficients.len(), n);
        let coefficients = (0..n)
            .map(|i| &a.coefficients[i] + &b.coefficients[i])
            .collect();
        RingOfIntegersWithIntegralBasisElement { coefficients }
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        todo!()
    }
}

impl RingStructure for RingOfIntegersWithIntegralBasisStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        todo!()
    }
}

impl UnitsStructure for RingOfIntegersWithIntegralBasisStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, crate::structure::RingDivisionError> {
        todo!()
    }
}

impl IntegralDomainStructure for RingOfIntegersWithIntegralBasisStructure {
    fn div(
        &self,
        a: &Self::Set,
        b: &Self::Set,
    ) -> Result<Self::Set, crate::structure::RingDivisionError> {
        todo!()
    }
}

impl FiniteUnitsStructure for RingOfIntegersWithIntegralBasisStructure {
    fn all_units(&self) -> Vec<Self::Set> {
        todo!()
    }
}
