use super::examples::c2::C2;
use crate::composition_table::group::MetaGenerateFiniteSubgroupTableSignature;
use algebraeon_macros::CanonicalStructure;
use algebraeon_structures::*;
use itertools::Itertools;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;

impl std::fmt::Display for Permutation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut cycles = self.disjoint_cycles();
        cycles.retain(|cycle| cycle.len() != 1);

        if cycles.is_empty() {
            f.write_str("()")?;
        }

        #[allow(clippy::redundant_closure_for_method_calls)]
        let string = cycles
            .iter()
            .map(|cycle| {
                "(".to_owned()
                    + &cycle
                        .cyc
                        .iter()
                        .map(|x| x.to_string())
                        .collect::<Vec<String>>()
                        .join(" ")
                    + ")"
            })
            .collect::<Vec<String>>()
            .join(" ");

        f.write_str(string.as_str())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_composition() {
        let a = Permutation::new(vec![0, 2, 1, 3]).unwrap();
        let b = Permutation::new(vec![1, 2, 0, 3]).unwrap();
        let c = Permutation::new(vec![2, 1, 0]).unwrap();

        println!("a = {}", a);
        println!("b = {}", b);
        println!("ab = {}", Permutation::compose(&a, &b));
        println!("c = {}", c);

        assert_eq!(Permutation::compose(&a, &b), c);
    }

    // #[test]
    // fn test_dcn() {
    //     let a = Permutation::new(vec![0, 2, 1, 5, 3, 4]).unwrap();
    //     assert_eq!(
    //         a.disjoint_cycles(),
    //         vec![Cycle::new(vec![1, 2]), Cycle::new(vec![3, 4, 5])]
    //     );
    // }

    #[test]
    fn test_sign() {
        let a = Permutation::new(vec![0, 2, 1]).unwrap();
        let b = Permutation::new(vec![1, 2, 0]).unwrap();
        let c = Permutation::new(vec![2, 1, 0]).unwrap();

        println!("a = {}", a);
        println!("b = {}", b);
        println!("c = {}", c);

        assert_eq!(a.sign(), C2::Flip);
        assert_eq!(b.sign(), C2::Identity);
        assert_eq!(c.sign(), C2::Flip);
    }

    #[test]
    fn test_all_permutations() {
        assert_eq!(Permutation::all_permutations(0).collect_vec().len(), 1);
        assert_eq!(Permutation::all_permutations(1).collect_vec().len(), 1);
        assert_eq!(Permutation::all_permutations(2).collect_vec().len(), 2);
        assert_eq!(Permutation::all_permutations(3).collect_vec().len(), 6);
        assert_eq!(Permutation::all_permutations(4).collect_vec().len(), 24);
        assert_eq!(Permutation::all_permutations(5).collect_vec().len(), 120);
    }
}
