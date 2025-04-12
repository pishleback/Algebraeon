use crate::{polynomial::*, structure::*};

pub mod conway_finite_fields;
pub mod conway_polynomials;
pub mod extension;
pub mod modulo;
pub mod polynomial;
pub mod quaternary_field;

impl<FS: FiniteFieldStructure> FactorableStructure for PolynomialStructure<FS> {
    fn factor(&self, p: &Self::Set) -> Option<crate::structure::Factored<Self>> {
        Some(
            self.factorize_monic(p)?
                .factorize_squarefree()
                .factorize_distinct_degree()
                .factorize_cantor_zassenhaus(),
        )
    }
}
