use crate::{polynomial::polynomial::*, structure::structure::*};

pub mod extension;
pub mod modulo;
pub mod polynomial;
pub mod quaternary_field;

impl<FS: FiniteFieldStructure> UniqueFactorizationStructure for PolynomialStructure<FS> {
    fn factor(&self, p: &Self::Set) -> Option<crate::structure::factorization::Factored<Self>> {
        Some(
            self.factorize_monic(p)?
                .factorize_squarefree()
                .factorize_distinct_degree()
                .factorize_cantor_zassenhaus(),
        )
    }
}
