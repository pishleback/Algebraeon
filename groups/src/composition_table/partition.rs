use std::collections::{BTreeSet, HashSet};

use super::group::*;
use super::subset::*;

#[derive(Clone, Debug)]
pub struct SetPartition {
    classes: Vec<BTreeSet<usize>>, //vector of conjugacy classes
    lookup: Vec<usize>,            //for each element, the index of its part class
}

impl SetPartition {
    pub fn check_state(&self, group: &FiniteGroup) -> Result<(), &'static str> {
        //is a partition
        let mut accounted_elems: HashSet<usize> = HashSet::new();
        for class in &self.classes {
            if class.len() == 0 {
                return Err("partition contains an empty set");
            }
            for x in class {
                if accounted_elems.contains(x) {
                    return Err("partition contains duplicate elements");
                }
                accounted_elems.insert(*x);
            }
        }
        if accounted_elems.len() != group.size() {
            return Err("partition is missing some elements");
        }

        //lookup is correct
        if !(self.lookup.len() == group.size()) {
            return Err("partition lookup has the wrong length");
        }
        for (elem, class_idx) in self.lookup.iter().enumerate() {
            if !(*class_idx < self.classes.len()) {
                return Err("partition lookup index is bigger partition size");
            }
            if !self.classes[*class_idx].contains(&elem) {
                return Err("partition lookup points to wrong partition set");
            }
        }
        Ok(())
    }

    pub fn new_unchecked(classes: Vec<BTreeSet<usize>>, lookup: Vec<usize>) -> Self {
        Self { classes, lookup }
    }

    pub fn project(&self, x: usize) -> &BTreeSet<usize> {
        &self.classes[self.lookup[x]]
    }
}

pub struct GroupPartition<'a> {
    pub group: &'a FiniteGroup,
    pub state: SetPartition,
    // is_left_cosets: Option<bool>,
    // is_right_cosets: Option<bool>,
}

impl<'a> PartialEq for GroupPartition<'a> {
    fn eq(&self, other: &Self) -> bool {
        let grp = self.group;
        if !std::ptr::eq(grp, other.group) {
            return false;
        }
        let n = self.size();
        if n != other.size() {
            return false;
        }
        //check equality by comparing lookup lists
        //equal iff lookups are equal up to permutation
        //build up a permutation f from indicies of self.lookup to indicies of other.lookup
        let mut f: Vec<Option<usize>> = vec![None; n];
        for i in 0..grp.size() {
            match f[self.state.lookup[i]] {
                Some(expected_other_lookup_i) => {
                    if expected_other_lookup_i != other.state.lookup[i] {
                        return false;
                    }
                }
                None => {
                    f[self.state.lookup[i]] = Some(other.state.lookup[i]);
                }
            }
        }
        true
    }
}

impl<'a> Eq for GroupPartition<'a> {}

impl<'a> GroupPartition<'a> {
    pub fn check_state(&self) -> Result<(), &'static str> {
        match self.state.check_state(self.group) {
            Err(msg) => {
                return Err(msg);
            }
            Ok(()) => {}
        };

        Ok(())
    }

    pub fn size(&self) -> usize {
        self.state.classes.len()
    }

    pub fn is_left_cosets(&self) -> bool {
        let ident_class = Subset::new_unchecked(
            self.group,
            self.state.classes[self.state.lookup[self.group.ident()]].clone(),
        );
        match ident_class.to_subgroup() {
            Some(ident_subgroup) => &ident_subgroup.left_cosets() == self,
            None => false,
        }
    }

    pub fn is_right_cosets(&self) -> bool {
        let ident_class = Subset::new_unchecked(
            self.group,
            self.state.classes[self.state.lookup[self.group.ident()]].clone(),
        );
        match ident_class.to_subgroup() {
            Some(ident_subgroup) => &ident_subgroup.right_cosets() == self,
            None => false,
        }
    }

    pub fn is_congruence(&self) -> bool {
        self.is_left_cosets() && self.is_right_cosets()
    }

    pub fn from_subsets(group: &'a FiniteGroup, subsets: Vec<BTreeSet<usize>>) -> Result<Self, ()> {
        let classes = subsets;
        let mut lookup = vec![0; group.size()];
        for (class_idx, class) in classes.iter().enumerate() {
            for x in class {
                if !(*x < group.size()) {
                    return Err(());
                }
                lookup[*x] = class_idx;
            }
        }
        let partition = Self {
            group,
            state: SetPartition { classes, lookup },
        };
        match partition.check_state() {
            Ok(_) => {}
            Err(_) => {
                return Err(());
            }
        }
        Ok(partition)
    }
}

pub struct Congruence<'a> {
    pub partition: GroupPartition<'a>,
}

impl<'a> Congruence<'a> {
    pub fn check_state(&self) -> Result<(), &'static str> {
        match self.partition.check_state() {
            Ok(()) => {}
            Err(msg) => {
                return Err(msg);
            }
        }

        match self.partition.is_congruence() {
            false => return Err("congruence is not a congruence"),
            true => {}
        }

        Ok(())
    }

    pub fn size(&self) -> usize {
        self.partition.size()
    }

    pub fn quotient_group(&self) -> FiniteGroup {
        let n = self.size();

        FiniteGroup::new_unchecked(
            n,
            self.partition.state.lookup[self.partition.group.ident()],
            {
                let mut inv = vec![0; n];
                for i in 0..n {
                    inv[i] = self.partition.state.lookup[self
                        .partition
                        .group
                        .inv(*self.partition.state.classes[i].iter().next().unwrap())];
                }
                inv
            },
            {
                let mut mul = vec![vec![0; n]; n];
                for i in 0..n {
                    for j in 0..n {
                        mul[i][j] = self.partition.state.lookup[self.partition.group.mul(
                            *self.partition.state.classes[i].iter().next().unwrap(),
                            *self.partition.state.classes[j].iter().next().unwrap(),
                        )];
                    }
                }
                mul
            },
            None,
            None,
        )
    }
}

#[cfg(test)]
mod partition_tests {
    use super::*;

    #[test]
    fn partition_check_bad_state() {
        let grp = examples::cyclic_group_structure(6);

        //elements too big
        let p = SetPartition {
            classes: vec![
                vec![0, 1, 2, 3].into_iter().collect(),
                vec![4, 5, 6].into_iter().collect(),
            ],
            lookup: vec![0, 0, 0, 0, 1, 1, 1],
        };
        match p.check_state(&grp) {
            Ok(()) => assert!(false),
            Err(_) => {}
        }

        //not a covering set
        let p = SetPartition {
            classes: vec![
                vec![0, 2].into_iter().collect(),
                vec![3, 5].into_iter().collect(),
            ],
            lookup: vec![0, 0, 0, 1, 1, 1],
        };
        match p.check_state(&grp) {
            Ok(()) => assert!(false),
            Err(_) => {}
        }

        //not disjoint
        let p = SetPartition {
            classes: vec![
                vec![0, 1, 2, 3].into_iter().collect(),
                vec![2, 3, 4, 5].into_iter().collect(),
            ],
            lookup: vec![0, 0, 0, 0, 1, 1],
        };
        match p.check_state(&grp) {
            Ok(()) => assert!(false),
            Err(_) => {}
        }

        //lookup values too big
        let p = SetPartition {
            classes: vec![
                vec![0, 1, 2].into_iter().collect(),
                vec![3, 4, 5].into_iter().collect(),
            ],
            lookup: vec![0, 0, 0, 1, 1, 2],
        };
        match p.check_state(&grp) {
            Ok(()) => assert!(false),
            Err(_) => {}
        }

        //incorrect lookup values
        let p = SetPartition {
            classes: vec![
                vec![0, 1, 2].into_iter().collect(),
                vec![3, 4, 5].into_iter().collect(),
            ],
            lookup: vec![0, 0, 1, 1, 1, 1],
        };
        match p.check_state(&grp) {
            Ok(()) => assert!(false),
            Err(_) => {}
        }
    }

    #[test]
    fn congruence_check_state() {
        let grp = examples::symmetric_group_structure(4);
        for (sg, _gens) in grp.subgroups() {
            match sg.to_normal_subgroup() {
                None => {
                    let cong: Congruence<'_> = Congruence {
                        partition: sg.left_cosets(),
                    };
                    match cong.check_state() {
                        Ok(_) => panic!(),
                        Err(_) => {}
                    }
                    let cong: Congruence<'_> = Congruence {
                        partition: sg.right_cosets(),
                    };
                    match cong.check_state() {
                        Ok(_) => panic!(),
                        Err(_) => {}
                    }
                }
                Some(nsg) => {
                    let cong = nsg.cosets();
                    match cong.check_state() {
                        Ok(_) => {}
                        Err(_) => panic!(),
                    }
                }
            }
        }
    }
}
