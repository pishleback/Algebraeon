use algebraeon_structures::{
    BorrowedStructure, FiniteSetSignature, MetaType, Natural, OptionalStructure, SetSignature,
    Signature,
};

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> FiniteSetSignature
    for OptionalStructure<Set, SetB>
{
    fn list_all_elements(&self) -> Vec<Self::Elem> {
        self.generate_all_elements().collect()
    }

    fn size(&self) -> Natural {
        self.set().size() + Natural::ONE
    }

    fn generate_random_elements(&self, seed: u64) -> impl Iterator<Item = Self::Elem> {
        let rng = StdRng::seed_from_u64(seed);
        FiniteSetRandomElementGenerator::<Self, StdRng> {
            all_elements: self.list_all_elements(),
            rng,
        }
    }
}
