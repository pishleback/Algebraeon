use crate::{polynomial::*, structure::*};
use algebraeon_sets::structure::BorrowedStructure;
pub mod conway_finite_fields;
pub mod conway_polynomials;
pub mod extension;
pub mod modulo;
pub mod polynomial;
pub mod quaternary_field;

impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> FactorableSignature
    for PolynomialStructure<FS, FSB>
{
    fn factor(&self, p: &Self::Set) -> Option<crate::structure::FactoredRingElement<Self>> {
        Some(
            self.factorize_monic(p)?
                .factorize_squarefree()
                .factorize_distinct_degree()
                .factorize_cantor_zassenhaus(),
        )
    }
}
