use crate::{polynomial::polynomial::*, ring_structure::quotient::*};

use super::*;

impl FiniteUnitsStructure for QuotientStructure<CannonicalStructure<Integer>, true> {
    fn all_units(&self) -> Vec<Self::Set> {
        let mut units = vec![];
        let mut u = Integer::from(1);
        while u < self.modulus().unsigned_abs() {
            units.push(u.clone());
            u += Integer::ONE;
        }
        units
    }
}

impl FiniteFieldStructure for QuotientStructure<CannonicalStructure<Integer>, true> {
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (self.modulus().unsigned_abs(), Natural::ONE)
    }
}

impl UniqueFactorizationStructure
    for PolynomialStructure<QuotientStructure<CannonicalStructure<Integer>, true>>
{
    fn factor(
        &self,
        p: &Self::Set,
    ) -> Option<crate::ring_structure::factorization::Factored<Self>> {
        Some(
            self.factorize_monic(p)?
                .factorize_squarefree()
                .factorize_berlekamps(),
        )
    }
}
