#![allow(dead_code)]

use std::{
    borrow::Borrow,
    collections::{BTreeSet, HashSet},
    rc::Rc,
};

use super::function::*;

pub trait SetT: Borrow<Set> + Clone {}
impl<T: Borrow<Set> + Clone> SetT for T {}

#[derive(Debug)]
pub struct Set {
    n: usize,
}

impl Set {
    pub fn new(n: usize) -> Self {
        Set { n }
    }

    pub fn size(&self) -> usize {
        self.n
    }

    pub fn elems(&self) -> std::ops::Range<usize> {
        0..self.n
    }

    pub fn contains(&self, x: usize) -> bool {
        x < self.n
    }
}

#[derive(Debug)]
pub struct Subset<'a> {
    set: &'a Set,
    elems: BTreeSet<usize>,
}

impl<'a> Subset<'a> {
    pub fn check_state(&self) -> Result<(), &'static str> {
        for x in &self.elems {
            if !self.set.contains(*x) {
                return Err("subset contains invalid set element");
            }
        }
        Ok(())
    }

    pub fn natural_inclusion(&self) -> Function<Rc<Set>, &'a Set> {
        Function::new_unchecked(
            Rc::new(Set {
                n: self.elems.len(),
            }),
            self.set,
            self.elems.clone().into_iter().collect(),
        )
    }
}

#[derive(Debug)]
pub struct Partition<'a> {
    set: &'a Set,
    parts: Vec<BTreeSet<usize>>,
    proj: Vec<usize>,
}

impl<'a> Partition<'a> {
    pub fn check_state(&self) -> Result<(), &'static str> {
        //is a partition
        let mut accounted_elems: HashSet<usize> = HashSet::new();
        for part in &self.parts {
            for x in part {
                if accounted_elems.contains(x) {
                    return Err("partition contains duplicate elements");
                }
                accounted_elems.insert(*x);
            }
        }
        if accounted_elems.len() != self.set.n {
            return Err("partition is missing some elements");
        }

        //lookup is correct
        if !(self.proj.len() == self.set.n) {
            return Err("partition lookup has the wrong length");
        }
        for (elem, part_idx) in self.proj.iter().enumerate() {
            if !(*part_idx < self.parts.len()) {
                return Err("partition lookup index is bigger partition size");
            }
            if !self.parts[*part_idx].contains(&elem) {
                return Err("partition lookup points to wrong partition set");
            }
        }
        Ok(())
    }

    pub fn natural_projection(&self) -> (Function<&'a Set, Rc<Set>>, &Vec<BTreeSet<usize>>) {
        (
            Function::new_unchecked(
                self.set,
                Rc::new(Set {
                    n: self.parts.len(),
                }),
                self.proj.clone(),
            ),
            &self.parts,
        )
    }
}

struct DisjointUnion<FirstT: SetT, SecondT: SetT> {
    first: FirstT,
    second: SecondT,
}

impl<FirstT: SetT, SecondT: SetT> DisjointUnion<FirstT, SecondT> {
    pub fn to_set(
        &self,
    ) -> (
        Rc<Set>,
        Function<FirstT, Rc<Set>>,
        Function<SecondT, Rc<Set>>,
    ) {
        let a = self.first.borrow().size();
        let b = self.second.borrow().size();
        let n = a + b;
        let set = Rc::new(Set { n });
        (
            set.clone(),
            Function::new_unchecked(
                self.first.clone(),
                set.clone(),
                (0..a).map(|x| 0 + x).collect(),
            ),
            Function::new_unchecked(
                self.second.clone(),
                set.clone(),
                (0..b).map(|x| a + x).collect(),
            ),
        )
    }

    pub fn size(&self) -> usize {
        self.first.borrow().size() + self.second.borrow().size()
    }
}

struct CartesianProduct<FirstT: SetT, SecondT: SetT> {
    first: FirstT,
    second: SecondT,
}

impl<FirstT: SetT, SecondT: SetT> CartesianProduct<FirstT, SecondT> {
    pub fn to_set(
        &self,
    ) -> (
        Rc<Set>,
        Function<Rc<Set>, FirstT>,
        Function<Rc<Set>, SecondT>,
    ) {
        let a = self.first.borrow().size();
        let b = self.second.borrow().size();
        let n = a * b;
        let set = Rc::new(Set { n });
        (
            set.clone(),
            Function::new_unchecked(
                set.clone(),
                self.first.clone(),
                (0..n).map(|x| x % a).collect(),
            ),
            Function::new_unchecked(
                set.clone(),
                self.second.clone(),
                (0..n).map(|x| x / a).collect(),
            ),
        )
    }

    pub fn size(&self) -> usize {
        self.first.borrow().size() * self.second.borrow().size()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_subset_state() {
        let set = Set { n: 0 };
        //empty subset of empty set
        let subset = Subset {
            set: &set,
            elems: vec![].into_iter().collect(),
        };
        assert!(subset.check_state().is_ok());
        //bad subset of empty set
        let subset = Subset {
            set: &set,
            elems: vec![0].into_iter().collect(),
        };
        assert!(subset.check_state().is_err());

        let set = Set { n: 7 };
        //typical
        let subset = Subset {
            set: &set,
            elems: vec![1, 2, 3].into_iter().collect(),
        };
        assert!(subset.check_state().is_ok());
        //empty
        let subset = Subset {
            set: &set,
            elems: vec![].into_iter().collect(),
        };
        assert!(subset.check_state().is_ok());
        //full
        let subset = Subset {
            set: &set,
            elems: vec![0, 1, 2, 3, 4, 5, 6].into_iter().collect(),
        };
        assert!(subset.check_state().is_ok());
        //partially invalid
        let subset = Subset {
            set: &set,
            elems: vec![5, 6, 7, 8].into_iter().collect(),
        };
        assert!(subset.check_state().is_err());
        //fully invalid
        let subset = Subset {
            set: &set,
            elems: vec![7, 8, 9].into_iter().collect(),
        };
        assert!(subset.check_state().is_err());
    }

    #[test]
    pub fn test_subset_natural_inclusion() {
        let set = Set { n: 7 };
        let subset = Subset {
            set: &set,
            elems: vec![2, 3, 4].into_iter().collect(),
        };
        let f = subset.natural_inclusion();
        f.check_state().unwrap();
        assert_eq!(f.domain().size(), 3);
        assert_eq!(f.range().size(), 7);
        assert!(std::ptr::eq(&set, f.range()));
        assert!(f.is_injective());
    }

    #[test]
    pub fn test_partition_state() {
        let set = Set { n: 7 };
        //typical
        let partition = Partition {
            set: &set,
            parts: vec![
                vec![0, 1].into_iter().collect(),
                vec![2, 3].into_iter().collect(),
                vec![4, 5, 6].into_iter().collect(),
            ]
            .into_iter()
            .collect(),
            proj: vec![0, 0, 1, 1, 2, 2, 2],
        };
        assert!(partition.check_state().is_ok());

        //not enough
        let partition = Partition {
            set: &set,
            parts: vec![
                vec![0, 1].into_iter().collect(),
                vec![2, 3].into_iter().collect(),
                vec![5, 6].into_iter().collect(),
            ]
            .into_iter()
            .collect(),
            proj: vec![0, 0, 1, 1, 2, 2, 2],
        };
        assert!(partition.check_state().is_err());

        //too many
        let partition = Partition {
            set: &set,
            parts: vec![
                vec![0, 1].into_iter().collect(),
                vec![2, 3, 4].into_iter().collect(),
                vec![4, 5, 6].into_iter().collect(),
            ]
            .into_iter()
            .collect(),
            proj: vec![0, 0, 1, 1, 2, 2, 2],
        };
        assert!(partition.check_state().is_err());

        //bad lookup
        let partition = Partition {
            set: &set,
            parts: vec![
                vec![0, 1].into_iter().collect(),
                vec![2, 3].into_iter().collect(),
                vec![4, 5, 6].into_iter().collect(),
            ]
            .into_iter()
            .collect(),
            proj: vec![0, 0, 1, 1, 1, 2, 2],
        };
        assert!(partition.check_state().is_err());
    }

    #[test]
    pub fn test_partition_natural_projection() {
        let set = Set { n: 7 };
        let partition = Partition {
            set: &set,
            parts: vec![
                vec![0, 1].into_iter().collect(),
                vec![2, 3].into_iter().collect(),
                vec![4, 5, 6].into_iter().collect(),
            ]
            .into_iter()
            .collect(),
            proj: vec![0, 0, 1, 1, 2, 2, 2],
        };
        let (f, _parts) = partition.natural_projection();
        f.check_state().unwrap();
        println!("{:?}", f);
        assert_eq!(f.domain().size(), 7);
        assert_eq!(f.range().size(), 3);
        assert!(std::ptr::eq(&set, f.domain()));
        assert!(f.is_surjective());
    }

    #[test]
    pub fn test_disjoint_union() {
        let a_set = Set { n: 3 };
        let b_set = Set { n: 5 };
        let union = DisjointUnion {
            first: &a_set,
            second: &b_set,
        };
        assert_eq!(union.size(), a_set.size() + b_set.size());
        let (union_set, inc_a, inc_b) = union.to_set();
        assert_eq!(union_set.size(), union.size());
        assert!(Rc::ptr_eq(&union_set, &inc_a.range()));
        assert!(Rc::ptr_eq(&union_set, &inc_b.range()));
        inc_a.check_state().unwrap();
        inc_b.check_state().unwrap();
        assert!(inc_a.is_injective());
        assert!(inc_b.is_injective());
    }

    #[test]
    pub fn test_cartesian_product() {
        let a_set = Set { n: 3 };
        let b_set = Set { n: 5 };
        let prod = CartesianProduct {
            first: &a_set,
            second: &b_set,
        };
        assert_eq!(prod.size(), a_set.size() * b_set.size());
        let (prod_set, proj_a, proj_b) = prod.to_set();
        assert_eq!(prod_set.size(), prod.size());
        assert!(Rc::ptr_eq(&prod_set, &proj_a.domain()));
        assert!(Rc::ptr_eq(&prod_set, &proj_b.domain()));
        proj_a.check_state().unwrap();
        proj_b.check_state().unwrap();
        assert!(proj_a.is_surjective());
        assert!(proj_b.is_surjective());
    }
}