
use std::collections::BTreeSet;
use std::collections::HashSet;

use super::partition::*;
use super::normal_subgroup::*;
use super::subset::*;
use super::group::*;
use std::hash::Hash;

#[derive(Clone)]
pub struct Subgroup<'a> {
    pub subset: Subset<'a>,
}

impl<'a> PartialEq for Subgroup<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.subset == other.subset
    }
}

impl<'a> Eq for Subgroup<'a> {}

impl<'a> Hash for Subgroup<'a> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.subset.hash(state);
    }
}

impl<'a> Subgroup<'a> {
    pub fn check_state(&self) -> Result<(), &'static str> {
        match self.subset.check_state() {
            Ok(_) => {}
            Err(msg) => {
                return Err(msg);
            }
        }

        if !self.subset.is_subgroup() {
            return Err("subgroup is not closed under composition");
        }

        Ok(())
    }

    pub fn size(&self) -> usize {
        self.subset.size()
    }

    pub fn to_group(&self) -> Group {
        let sg_elems: Vec<usize> = self.subset.elems.clone().into_iter().collect();
        let k = sg_elems.len();
        let mut group_to_subgroup: Vec<Option<usize>> = vec![None; self.subset.group.size()];
        for (i, x) in sg_elems.iter().enumerate() {
            group_to_subgroup[*x] = Some(i);
        }
        //TODO: add a test that this group has valid structure
        Group {
            n: self.size(),
            ident: group_to_subgroup[self.subset.group.ident].unwrap(),
            inv: (0..k)
                .map(|x| group_to_subgroup[self.subset.group.inv[sg_elems[x]]].unwrap())
                .collect(),
            mul: (0..k)
                .map(|x| {
                    (0..k)
                        .map(|y| {
                            group_to_subgroup[self.subset.group.mul[sg_elems[x]][sg_elems[y]]]
                                .unwrap()
                        })
                        .collect()
                })
                .collect(),
            conjugacy_classes: None,
            is_abelian: None,
            is_simple: None,
        }
    }

    pub fn left_cosets(&self) -> Partition {
        let mut cosets: Vec<BTreeSet<usize>> = vec![];
        let mut coset_lookup = vec![0; self.subset.group.size()];
        let mut missing: HashSet<usize> = self.subset.group.elems().collect();
        while missing.len() > 0 {
            let g = missing.iter().next().unwrap();
            let coset = self.subset.left_mul(*g).unwrap().elems;
            for h in &coset {
                missing.remove(h);
                coset_lookup[*h] = cosets.len();
            }
            cosets.push(coset);
        }
        Partition {
            group: self.subset.group,
            state: PartitionState {
                classes: cosets,
                lookup: coset_lookup,
            },
        }
    }

    pub fn right_cosets(&self) -> Partition {
        let mut cosets: Vec<BTreeSet<usize>> = vec![];
        let mut coset_lookup = vec![0; self.subset.group.size()];
        let mut missing: HashSet<usize> = self.subset.group.elems().collect();
        while missing.len() > 0 {
            let g = missing.iter().next().unwrap();
            let coset = self.subset.right_mul(*g).unwrap().elems;
            for h in &coset {
                missing.remove(h);
                coset_lookup[*h] = cosets.len();
            }
            cosets.push(coset);
        }
        Partition {
            group: self.subset.group,
            state: PartitionState {
                classes: cosets,
                lookup: coset_lookup,
            },
        }
    }

    pub fn is_normal_subgroup(&self) -> bool {
        for x in &self.subset.elems {
            for g in self.subset.group.elems() {
                if !self.subset.elems.contains(
                    &self.subset.group.mul[self.subset.group.mul[g][*x]][self.subset.group.inv[g]],
                )
                //gxg^{-1}
                {
                    return false;
                }
            }
        }
        true
    }

    pub fn to_normal_subgroup(&self) -> Option<NormalSubgroup> {
        match self.is_normal_subgroup() {
            true => Some(NormalSubgroup {
                subgroup: self.clone(),
            }),
            false => None,
        }
    }
}

#[cfg(test)]
mod subgroup_tests {
    use super::*;
    use super::super::super::permutations::*;
    
    #[test]
    fn subgroup_counts() {
        for (grp, num_sgs) in vec![
            (cyclic_group_structure(12), 6),
            (cyclic_group_structure(13), 2),
            (dihedral_group_structure(1), 2),
            (dihedral_group_structure(12), 34),
            (dihedral_group_structure(13), 16),
            (symmetric_group_structure(1).0, 1),
            (symmetric_group_structure(2).0, 2),
            (symmetric_group_structure(3).0, 6),
            (symmetric_group_structure(4).0, 30),
            (symmetric_group_structure(5).0, 156),
        ] {
            assert_eq!(grp.subgroups().len(), num_sgs);
        }
    }

    #[test]
    fn subgroup_state() {
        //are the states of subgroups, normal subgroups, and generating sets correct when produced by grp.subgroups() and grp.normal_subgroups()?
        for grp in vec![
            symmetric_group_structure(4).0,
            dihedral_group_structure(12),
            cyclic_group_structure(8),
        ] {
            let sgs = grp.subgroups();
            for (sg, gens) in sgs {
                // sg.is_normal_subgroup(); //cache is to check if the coset partition has the right is_congruence value
                sg.check_state().unwrap();
                gens.check_state().unwrap();
                // sg.left_cosets().unwrap().check_state().unwrap();
                // sg.right_cosets().unwrap().check_state().unwrap();
            }
        }
    }


    #[test]
    fn subgroup_left_cosets() {
        let (grp, _perms, elems) = symmetric_group_structure(3);
        let sg = Subgroup {
            subset: Subset {
                group: &grp,
                elems: vec![
                    elems[&Permutation::new(vec![0, 1, 2]).unwrap()],
                    elems[&Permutation::new(vec![1, 0, 2]).unwrap()],
                ]
                .into_iter()
                .collect(),
            },
        };
        sg.check_state().unwrap();
        let cosets = sg.left_cosets();
        cosets.check_state().unwrap();
        assert_eq!(cosets.size(), 3);
        assert!(
            cosets
                == Partition::from_subsets(
                    &grp,
                    vec![
                        vec![
                            elems[&Permutation::new(vec![0, 1, 2]).unwrap()],
                            elems[&Permutation::new(vec![1, 0, 2]).unwrap()],
                        ]
                        .into_iter()
                        .collect(),
                        vec![
                            elems[&Permutation::new(vec![1, 2, 0]).unwrap()],
                            elems[&Permutation::new(vec![2, 1, 0]).unwrap()],
                        ]
                        .into_iter()
                        .collect(),
                        vec![
                            elems[&Permutation::new(vec![2, 0, 1]).unwrap()],
                            elems[&Permutation::new(vec![0, 2, 1]).unwrap()],
                        ]
                        .into_iter()
                        .collect()
                    ]
                )
                .unwrap()
        );
    }

    #[test]
    fn subgroup_right_cosets() {
        let (grp, _perms, elems) = symmetric_group_structure(3);
        let sg = Subgroup {
            subset: Subset {
                group: &grp,
                elems: vec![
                    elems[&Permutation::new(vec![0, 1, 2]).unwrap()],
                    elems[&Permutation::new(vec![1, 0, 2]).unwrap()],
                ]
                .into_iter()
                .collect(),
            },
        };
        sg.check_state().unwrap();
        let cosets = sg.right_cosets();
        cosets.check_state().unwrap();
        assert_eq!(cosets.size(), 3);
        assert!(
            cosets
                == Partition::from_subsets(
                    &grp,
                    vec![
                        vec![
                            elems[&Permutation::new(vec![0, 1, 2]).unwrap()],
                            elems[&Permutation::new(vec![1, 0, 2]).unwrap()],
                        ]
                        .into_iter()
                        .collect(),
                        vec![
                            elems[&Permutation::new(vec![1, 2, 0]).unwrap()],
                            elems[&Permutation::new(vec![0, 2, 1]).unwrap()],
                        ]
                        .into_iter()
                        .collect(),
                        vec![
                            elems[&Permutation::new(vec![2, 0, 1]).unwrap()],
                            elems[&Permutation::new(vec![2, 1, 0]).unwrap()],
                        ]
                        .into_iter()
                        .collect()
                    ]
                )
                .unwrap()
        );
    }
}