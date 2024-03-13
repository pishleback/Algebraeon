use super::group::*;
use super::partition::*;
use super::subgroup::*;

pub struct NormalSubgroup<'a> {
    subgroup: Subgroup<'a>,
}

impl<'a> NormalSubgroup<'a> {
    pub fn check_state(&self) -> Result<(), &'static str> {
        match self.subgroup.check_state() {
            Ok(_) => {}
            Err(msg) => {
                return Err(msg);
            }
        }

        if !(self.subgroup.is_normal_subgroup()) {
            return Err("normal subgroup is not norpppmal");
        }

        Ok(())
    }

    pub fn new_unchecked(subgroup: Subgroup<'a>) -> Self {
        Self { subgroup }
    }

    pub fn size(&self) -> usize {
        self.subgroup.size()
    }

    pub fn subgroup(&self) -> &Subgroup<'a> {
        &self.subgroup
    }

    pub fn to_subgroup(self) -> Subgroup<'a> {
        self.subgroup
    }

    pub fn cosets(&self) -> Congruence {
        Congruence {
            partition: self.subgroup.left_cosets(),
        } //should be the same as right cosets
    }

    pub fn quotient_group(&self) -> Group {
        self.cosets().quotient_group()
    }
}

#[cfg(test)]
mod normal_subgroup_tests {
    use super::super::super::sets::permutations::*;
    use super::super::subset::*;
    use super::*;

    #[test]
    fn normal_subgroup_state() {
        //are the states of subgroups, normal subgroups, and generating sets correct when produced by grp.subgroups() and grp.normal_subgroups()?
        for mut grp in vec![
            symmetric_group_structure(4).0,
            dihedral_group_structure(12),
            cyclic_group_structure(8),
        ] {
            grp.cache_conjugacy_classes();
            let nsgs = grp.normal_subgroups();
            for (nsg, gens) in nsgs {
                nsg.check_state().unwrap();
                gens.check_state().unwrap();
                // let left_cosets = nsg.left_cosets().unwrap();
                // match left_cosets.is_right_cosets {
                //     Some(flag) => {
                //         assert!(flag)
                //     }
                //     None => {
                //         assert!(false)
                //     }
                // }
            }
        }
    }

    #[test]
    fn normal_subgroup_counts() {
        for (mut grp, num_sgs) in vec![
            (cyclic_group_structure(12), 6),  //Abelian so same as # of subgroups
            (cyclic_group_structure(13), 2),  //Abelian so same as # of subgroups
            (dihedral_group_structure(1), 2), //Abelian so same as # of subgroups
            (dihedral_group_structure(12), 9),
            (dihedral_group_structure(13), 3),   //{e}, C13, D13
            (symmetric_group_structure(1).0, 1), //trivial group
            (symmetric_group_structure(2).0, 2), //cyclic of order 2
            (symmetric_group_structure(3).0, 3), //{e}, C3, S3
            (symmetric_group_structure(4).0, 4), //{e}, V4, A4, S4
            (symmetric_group_structure(5).0, 3), //{e}, A5, S5
        ] {
            grp.cache_conjugacy_classes();
            assert_eq!(grp.normal_subgroups().len(), num_sgs);
        }
    }

    #[test]
    fn normal_subgroup_cosets() {
        let (grp, _perms, elems) = symmetric_group_structure(3);
        let nsg = NormalSubgroup {
            subgroup: Subgroup {
                subset: Subset::new_unchecked(
                    &grp,
                    vec![
                        elems[&Permutation::new(vec![0, 1, 2]).unwrap()],
                        elems[&Permutation::new(vec![1, 2, 0]).unwrap()],
                        elems[&Permutation::new(vec![2, 0, 1]).unwrap()],
                    ]
                    .into_iter()
                    .collect(),
                ),
            },
        };
        nsg.check_state().unwrap();
        let left_cosets = nsg.subgroup.left_cosets();
        let right_cosets = nsg.subgroup.right_cosets();
        assert!(left_cosets == right_cosets);
        let cosets = nsg.cosets();
        cosets.check_state().unwrap();
        assert!(cosets.partition == left_cosets);
        assert!(cosets.partition == right_cosets);
    }

    #[test]
    fn normal_subgroup_quotient() {
        let (mut grp, _perms, _elems) = symmetric_group_structure(4);
        grp.cache_conjugacy_classes();
        for (nsg, _gens) in grp.normal_subgroups() {
            let qgrp = nsg.quotient_group();
            qgrp.check_state().unwrap();
            assert_eq!(qgrp.size() * nsg.size(), grp.size());
        }
    }
}
