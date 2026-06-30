use crate::*;
use algebraeon_macros::{signature_meta_trait, skip_meta};
use std::collections::HashMap;

#[signature_meta_trait]
pub trait PermutationsSignature<Set: SetSignature>: GroupSignature {
    #[skip_meta]
    fn set(&self) -> &Set;

    /// Error if a == b
    #[allow(clippy::result_unit_err)]
    fn new_swap(&self, a: Set::Elem, b: Set::Elem) -> Result<Self::Elem, ()> {
        self.new_cycle(vec![a, b])
    }

    /// Error if any two elements of cycle are equal
    #[allow(clippy::result_unit_err)]
    fn new_cycle(&self, cycle: Vec<Set::Elem>) -> Result<Self::Elem, ()>;

    /// Error if any cycle is invalid
    #[allow(clippy::result_unit_err)]
    fn new_cycles(&self, cycles: Vec<Vec<Set::Elem>>) -> Result<Self::Elem, ()> {
        Ok(self.compose_list(
            cycles
                .into_iter()
                .map(|cycle| self.new_cycle(cycle).unwrap())
                .collect(),
        ))
    }

    /// Error if not valid permutation
    #[allow(clippy::result_unit_err)]
    fn new_perm(
        &self,
        pairs: Vec<(impl Borrow<Set::Elem>, impl Borrow<Set::Elem>)>,
    ) -> Result<Self::Elem, ()>;

    fn image(&self, perm: &Self::Elem, elem: &Set::Elem) -> Set::Elem;

    fn preimage(&self, perm: &Self::Elem, elem: &Set::Elem) -> Set::Elem;

    /// All elements not fixed by this permutation
    fn support(&self, perm: Self::Elem) -> Vec<Set::Elem>;

    /// How many elements are moved by the permutation
    fn support_size(&self, perm: &Self::Elem) -> usize;

    /// The disjoint cycle decomposition
    /// Cycles map appear in any order and the elements of each cycle may come in any order
    fn disjoint_cycles(&self, perm: &Self::Elem) -> Vec<Vec<Set::Elem>>;

    fn cycle_shape(&self, perm: &Self::Elem) -> HashMap<usize, usize> {
        let mut cycle_counts = HashMap::new();
        for cycle in self.disjoint_cycles(perm) {
            let k = cycle.len();
            *cycle_counts.entry(k).or_insert(0) += 1;
        }
        cycle_counts
    }

    fn is_even(&self, perm: &Self::Elem) -> bool {
        let mut is_even = true;
        for cycle in self.disjoint_cycles(perm) {
            if cycle.len() % 2 == 0 {
                is_even = !is_even;
            }
        }
        is_even
    }

    fn is_odd(&self, perm: &Self::Elem) -> bool {
        !self.is_even(perm)
    }
}
