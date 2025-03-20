use crate::structure::quotient::*;

use super::*;

impl FiniteUnitsStructure for QuotientStructure<CannonicalStructure<Integer>, true> {
    fn all_units(&self) -> Vec<Self::Set> {
        let mut units = vec![];
        let mut u = Integer::from(1);
        while u < self.modulus().abs() {
            units.push(u.clone());
            u += Integer::ONE;
        }
        units
    }
}

impl FiniteFieldStructure for QuotientStructure<CannonicalStructure<Integer>, true> {
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (self.modulus().abs(), Natural::ONE)
    }
}
