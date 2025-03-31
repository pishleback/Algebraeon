use algebraeon_sets::combinatorics::Partition;

use super::group::*;
use super::subset::*;
use std::collections::HashSet;

#[derive(Debug)]
pub struct GroupPartition<'a> {
    pub group: &'a FiniteGroupMultiplicationTable,
    pub partition: Partition,
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
            match f[self.partition.project(i)] {
                Some(expected_other_lookup_i) => {
                    if expected_other_lookup_i != other.partition.project(i) {
                        return false;
                    }
                }
                None => {
                    f[self.partition.project(i)] = Some(other.partition.project(i));
                }
            }
        }
        true
    }
}

impl<'a> Eq for GroupPartition<'a> {}

impl<'a> GroupPartition<'a> {
    pub fn check_state(&self) -> Result<(), &'static str> {
        if self.partition.num_elements() != self.group.size() {
            return Err("Partition is of set of the wrong size");
        }
        Ok(())
    }

    pub fn size(&self) -> usize {
        self.partition.size()
    }

    pub fn is_left_cosets(&self) -> bool {
        let ident_class = Subset::new_unchecked(
            self.group,
            self.partition.class_containing(self.group.ident()).clone(),
        );
        match ident_class.to_subgroup() {
            Some(ident_subgroup) => &ident_subgroup.left_cosets() == self,
            None => false,
        }
    }

    pub fn is_right_cosets(&self) -> bool {
        let ident_class = Subset::new_unchecked(
            self.group,
            self.partition.class_containing(self.group.ident()).clone(),
        );
        match ident_class.to_subgroup() {
            Some(ident_subgroup) => &ident_subgroup.right_cosets() == self,
            None => false,
        }
    }

    pub fn is_congruence(&self) -> bool {
        self.is_left_cosets() && self.is_right_cosets()
    }

    pub fn from_subsets(
        group: &'a FiniteGroupMultiplicationTable,
        subsets: Vec<HashSet<usize>>,
    ) -> Result<Self, ()> {
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
            partition: Partition::new_unchecked(classes, lookup),
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

    pub fn quotient_group(&self) -> FiniteGroupMultiplicationTable {
        let n = self.size();
        self.partition.group.check_state().unwrap();
        FiniteGroupMultiplicationTable::new_unchecked(
            n,
            self.partition
                .partition
                .project(self.partition.group.ident()),
            {
                let mut inv = vec![0; n];
                for i in 0..n {
                    inv[i] = self.partition.partition.project(
                        self.partition
                            .group
                            .inv(*self.partition.partition.get_class(i).iter().next().unwrap()),
                    );
                }
                inv
            },
            {
                let mut mul = vec![vec![0; n]; n];
                for i in 0..n {
                    for j in 0..n {
                        mul[i][j] = self.partition.partition.project(self.partition.group.mul(
                            *self.partition.partition.get_class(i).iter().next().unwrap(),
                            *self.partition.partition.get_class(j).iter().next().unwrap(),
                        ));
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
