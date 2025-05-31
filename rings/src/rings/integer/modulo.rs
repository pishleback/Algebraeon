use crate::rings::quotient::QuotientStructure;

use super::*;

impl FiniteUnitsSignature for QuotientStructure<IntegerCanonicalStructure, true> {
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

impl FiniteFieldSignature for QuotientStructure<IntegerCanonicalStructure, true> {
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (self.modulus().abs(), Natural::ONE)
    }
}

impl<const IS_FIELD: bool> CountableSetSignature
    for QuotientStructure<IntegerCanonicalStructure, IS_FIELD>
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> {
        (0usize..)
            .map(Integer::from)
            .take_while(|n| n < self.modulus())
    }
}

impl<const IS_FIELD: bool> FiniteSetSignature
    for QuotientStructure<IntegerCanonicalStructure, IS_FIELD>
{
    fn size(&self) -> usize {
        self.modulus().abs().try_into().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_elements() {
        assert_eq!(
            QuotientStructure::new_ring(Integer::structure(), 26.into())
                .list_all_elements()
                .len(),
            26
        );
    }
}
