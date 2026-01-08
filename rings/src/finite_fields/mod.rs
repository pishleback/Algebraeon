use crate::{polynomial::*, structure::*};
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::*;
pub mod conway_finite_fields;
pub mod conway_polynomials;
pub mod extension;
pub mod modulo;
pub mod polynomial;
pub mod quaternary_field;

impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> FactoringMonoidSignature
    for PolynomialStructure<FS, FSB>
{
    fn factor_unchecked(&self, p: &Self::Set) -> Factored<Self::Set, Natural> {
        self.factorize_monic(p)?
            .factorize_squarefree()
            .factorize_distinct_degree()
            .factorize_cantor_zassenhaus()
    }
}
