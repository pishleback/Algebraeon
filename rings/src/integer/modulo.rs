use super::*;

impl<
    B: BorrowedStructure<IntegerCanonicalStructure>,
    BE: BorrowedStructure<EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, true>>,
> CountableSetSignature
    for MultiplicativeMonoidUnitsStructure<
        EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, true>,
        BE,
    >
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
        self.list_all_elements().into_iter()
    }
}

impl<
    B: BorrowedStructure<IntegerCanonicalStructure>,
    BE: BorrowedStructure<EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, true>>,
> FiniteSetSignature
    for MultiplicativeMonoidUnitsStructure<
        EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, true>,
        BE,
    >
{
    fn list_all_elements(&self) -> Vec<Self::Set> {
        let mut units = vec![];
        let mut u = Integer::from(1);
        while u < Abs::abs(self.monoid().modulus()) {
            units.push(u.clone());
            u += Integer::ONE;
        }
        units
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> FiniteFieldSignature
    for EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, true>
{
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (Abs::abs(self.modulus()), Natural::ONE)
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>, const IS_FIELD: bool> CountableSetSignature
    for EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, IS_FIELD>
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> + Clone {
        (0usize..)
            .map(Integer::from)
            .take_while(|n| n < self.modulus())
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>, const IS_FIELD: bool> FiniteSetSignature
    for EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, IS_FIELD>
{
    fn size(&self) -> usize {
        Abs::abs(self.modulus()).try_into().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn count_elements() {
        assert_eq!(
            Integer::structure()
                .into_quotient_ring(26.into())
                .unwrap()
                .list_all_elements()
                .len(),
            26
        );
    }
}
