use super::number_field::AlgebraicNumberFieldStructure;
use crate::polynomial::Polynomial;
use algebraeon_nzq::{Integer, Rational};

#[derive(Debug, Clone)]
pub struct RingOfIntegersStructure {
    algebraic_number_field: AlgebraicNumberFieldStructure,
    integral_basis: Vec<Polynomial<Rational>>,
    discriminant: Integer,
}

impl RingOfIntegersStructure {
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
}
