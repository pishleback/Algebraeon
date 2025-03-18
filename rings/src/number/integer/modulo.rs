use crate::{polynomial::polynomial::*, structure::quotient::*};

use super::*;

impl FiniteUnitsStructure for QuotientStructure<CannonicalStructure<Integer>, true> {
    fn all_units(&self) -> Vec<Self::Set> {
        let mut units = vec![];
        let mut u = Integer::from(1);
        while u < self.modulus().unsigned_abs_ref() {
            units.push(u.clone());
            u += Integer::ONE;
        }
        units
    }
}

impl FiniteFieldStructure for QuotientStructure<CannonicalStructure<Integer>, true> {
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (self.modulus().unsigned_abs_ref(), Natural::ONE)
    }
}

impl UniqueFactorizationStructure
    for PolynomialStructure<QuotientStructure<CannonicalStructure<Integer>, true>>
{
    fn factor(&self, p: &Self::Set) -> Option<crate::structure::factorization::Factored<Self>> {
        Some(
            self.factorize_monic(p)?
                .factorize_squarefree()
                .factorize_distinct_degree()
                .factorize_cantor_zassenhaus(),
        )
    }
}
