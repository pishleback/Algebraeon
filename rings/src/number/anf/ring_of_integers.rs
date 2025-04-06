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
        debug_assert_eq!(elem.coefficients.len(), self.degree());
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        polynomial::PolynomialStructure,
        structure::{IntoErgonomic, MetaSemiRing},
    };

    #[test]
    fn ring_of_integer_arithmetic() {
        let x = Polynomial::<Rational>::var().into_ergonomic();

        // Take the integral basis (0 + x, 1/2 + 1/2x)
        let a = Polynomial::<Rational>::from_coeffs(vec![Rational::ZERO, Rational::ONE]);
        let b = Polynomial::<Rational>::from_coeffs(vec![Rational::ONE_HALF, Rational::ONE_HALF]);

        let anf = (x.pow(2) + 7).into_verbose().algebraic_number_field();
        let roi = RingOfIntegersWithIntegralBasisStructure::new(
            anf.clone(),
            vec![a.clone(), b.clone()],
            Integer::from(-7),
        );

        {
            assert_eq!(
                roi.roi_to_anf(RingOfIntegersWithIntegralBasisElement {
                    coefficients: vec![Integer::from(1), Integer::from(4)]
                }),
                (2 + 3 * x).into_verbose()
            );
        }

        {
            assert!(
                roi.anf_to_roi(Polynomial::<Rational>::from_coeffs(vec![
                    Rational::ONE_HALF,
                    Rational::ONE,
                ]))
                .is_none()
            );

            let c = roi
                .anf_to_roi(Polynomial::<Rational>::from_coeffs(vec![
                    Rational::from(2),
                    Rational::from(3),
                ]))
                .unwrap()
                .coefficients;
            assert_eq!(c.len(), 2);
            assert_eq!(c[0], Integer::from(1));
            assert_eq!(c[1], Integer::from(4));
        }

        {
            // 0 = 0 * (0+x) + 0 * (1/2 + 1/2x)
            let zero = roi.zero().coefficients;
            assert_eq!(zero.len(), 2);
            assert_eq!(zero[0], Integer::ZERO);
            assert_eq!(zero[1], Integer::ZERO);
        }

        {
            // 1 = -1 * (0+x) + 2 * (1/2 + 1/2x)
            let one = roi.one().coefficients;
            assert_eq!(one.len(), 2);
            assert_eq!(one[0], Integer::from(-1));
            assert_eq!(one[1], Integer::from(2));
        }
    }
}
