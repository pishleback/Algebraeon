use crate::{polynomial::*, structure::*};
use std::sync::Arc;
pub mod conway_finite_fields;
pub mod conway_polynomials;
pub mod extension;
pub mod polynomial;
pub mod quaternary_field;
use algebraeon_structures::*;

impl<FS: FiniteFieldSignature> FactoringMonoidSignature for PolynomialStructure<FS> {
    fn factor_unchecked(self: &Arc<Self>, p: &Self::Elem) -> Factored<Self::Elem, Natural> {
        if let Some(p) = self.factorize_monic(p) {
            p.factorize_squarefree()
                .factorize_distinct_degree()
                .factorize_cantor_zassenhaus()
        } else {
            Factored::Zero
        }
    }
}
