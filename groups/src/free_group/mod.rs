pub mod todd_coxeter;

use malachite_base::num::basic::traits::{One, Zero};
use malachite_nz::integer::Integer;
use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::Hash;

mod underlying_set {
    use super::*;
    pub trait UnderlyingSet: Debug + Clone + PartialEq + Eq + Hash {}
}
use underlying_set::*;
impl<T: Debug + Clone + PartialEq + Eq + Hash> UnderlyingSet for T {}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FreeGroupTerm<S: UnderlyingSet> {
    generator: S,
    power: Integer,
}
impl<S: UnderlyingSet> FreeGroupTerm<S> {
    fn inverse(self) -> Self {
        Self {
            generator: self.generator,
            power: -self.power,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum FreeGroupAtom<S: UnderlyingSet> {
    Identity(S),
    Inverse(S),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FreeGroup<S: UnderlyingSet> {
    reduced_word: Vec<FreeGroupTerm<S>>,
}

impl<S: UnderlyingSet> FreeGroup<S> {
    pub fn generator(s: S) -> Self {
        Self {
            reduced_word: vec![FreeGroupTerm {
                generator: s,
                power: Integer::ONE,
            }],
        }
    }

    fn into_atoms(&self) -> Vec<FreeGroupAtom<S>> {
        let mut atoms = vec![];
        for term in &self.reduced_word {
            debug_assert_ne!(term.power, Integer::ZERO);
            let mut i = Integer::ZERO;
            while i < term.power {
                atoms.push(FreeGroupAtom::Identity(term.generator.clone()));
                i += Integer::ONE;
            }
            while i > term.power {
                atoms.push(FreeGroupAtom::Inverse(term.generator.clone()));
                i -= Integer::ONE;
            }
            debug_assert_eq!(i, term.power);
        }
        atoms
    }
}

impl<S: UnderlyingSet> super::group::Group for FreeGroup<S> {
    fn identity() -> Self {
        Self {
            reduced_word: vec![],
        }
    }

    fn inverse(self) -> Self {
        Self {
            reduced_word: self
                .reduced_word
                .into_iter()
                .rev()
                .map(|term| term.inverse())
                .collect(),
        }
    }

    fn compose_mut(&mut self, other: &Self) {
        for other_term in &other.reduced_word {
            match self.reduced_word.last_mut() {
                Some(self_term) => {
                    if self_term.generator == other_term.generator {
                        if &self_term.power + &other_term.power == 0 {
                            self.reduced_word.pop().unwrap();
                        } else {
                            self_term.power += &other_term.power;
                        }
                    } else {
                        self.reduced_word.push(other_term.clone());
                    }
                }
                None => {
                    self.reduced_word.push(other_term.clone());
                }
            }
        }
    }
}

impl<S: UnderlyingSet> FreeGroup<S> {
    pub fn generated_finite_group(
        gens: Vec<S>,
        relations: Vec<Self>,
    ) -> super::composition_table::group::Group {
        let gen_to_idx: HashMap<&S, usize> = gens.iter().enumerate().map(|(s, i)| (i, s)).collect();
        todd_coxeter::enumerate_group(
            gens.len(),
            relations
                .iter()
                .map(|relation| {
                    relation
                        .into_atoms()
                        .into_iter()
                        .map(|a| match a {
                            FreeGroupAtom::Identity(g) => 2 * gen_to_idx.get(&g).unwrap(),
                            FreeGroupAtom::Inverse(g) => 2 * gen_to_idx.get(&g).unwrap() + 1,
                        })
                        .collect()
                })
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::natural::Natural;

    use crate::group::Group;

    use super::*;

    #[test]
    fn test_int_pow() {
        let a = FreeGroup::generator("a");
        assert_eq!(
            a.int_pow(&Integer::from(-100)),
            FreeGroup::compose_list(vec![
                a.int_pow(&Integer::from(-200)),
                a.int_pow(&Integer::from(100))
            ])
        );
    }

    #[test]
    fn test_generate_finite_group() {
        let a = FreeGroup::generator("a");
        let b = FreeGroup::generator("b");
        let c = FreeGroup::generator("c");
        let g = FreeGroup::generated_finite_group(
            vec!["a", "b", "c"],
            vec![
                FreeGroup::compose_list(vec![&a, &a]),
                FreeGroup::compose_list(vec![&b, &b]),
                FreeGroup::compose_list(vec![&c, &c]),
                FreeGroup::compose_list(vec![&a, &b]).nat_pow(&Natural::from(3u32)),
                FreeGroup::compose_list(vec![&b, &c]).nat_pow(&Natural::from(5u32)),
                FreeGroup::compose_list(vec![&a, &c]).nat_pow(&Natural::from(2u32)),
            ],
        );
        assert_eq!(g.size(), 120);
    }
}
