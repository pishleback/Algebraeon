use std::collections::{BTreeSet, HashSet};
use std::rc::Rc;

use super::function::*;
use super::set::*;

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
        if accounted_elems.len() != self.set.size() {
            return Err("partition is missing some elements");
        }

        //lookup is correct
        if !(self.proj.len() == self.set.size()) {
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
                Rc::new(Set::new(self.parts.len())),
                self.proj.clone(),
            ),
            &self.parts,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_partition_state() {
        let set = Set::new(7);
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
        let set = Set::new(7);
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
}
