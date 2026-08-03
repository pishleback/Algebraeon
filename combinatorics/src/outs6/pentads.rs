use crate::outs6::*;
use algebraeon_macros::signature_meta_trait;
use algebraeon_sets::sets::{
    FiniteSetToFinitelySupportedPermutationsStructure, FinitelySupportedPermutation,
};
use algebraeon_structures::*;
use std::{cmp::Ordering, sync::Arc};

#[derive(Debug, Clone)]
pub struct Pentad<Elem> {
    // must have syntheme_1 < syntheme_2 < syntheme_3 < syntheme_4 < syntheme_5 and all disjoint
    pub(crate) synthemes: [Syntheme<Elem>; 5],
}

/// The 15-element set of duads on a 6-element set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PentadsStructure<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> {
    set: Arc<Set>,
}

impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> PentadsStructure<Set> {
    pub fn new(set: Arc<Set>) -> Arc<Self> {
        debug_assert_eq!(set.size(), Natural::from(6usize));
        Self { set }.into()
    }

    pub fn set(&self) -> &Arc<Set> {
        &self.set
    }
}

pub trait SetToPentadsSignature:
    ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature
{
    fn pentads(self: &Arc<Self>) -> Arc<PentadsStructure<Self>> {
        PentadsStructure::new(self.clone())
    }
}
impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> SetToPentadsSignature
    for Set
{
}

impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> Signature
    for PentadsStructure<Set>
{
}

impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> SetSignature
    for PentadsStructure<Set>
{
    type Elem = Pentad<Set::Elem>;

    fn validate_element(self: &Arc<Self>, p: &Self::Elem) -> Result<(), String> {
        let synthemes = self.set().synthemes();
        for s in &p.synthemes {
            synthemes.validate_element(s)?;
        }
        if !synthemes.is_sorted(&p.synthemes) {
            return Err("synthemes are not sorted".to_string());
        }
        for i in 0..5 {
            for j in (i + 1)..5 {
                if !synthemes
                    .overlap(&p.synthemes[i], &p.synthemes[j])
                    .is_disjoint()
                {
                    return Err("synthemes are not disjoint".to_string());
                }
            }
        }
        Ok(())
    }
}

impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> EqSignature
    for PentadsStructure<Set>
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.cmp(a, b).is_eq()
    }
}

impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> PartialOrdSignature
    for PentadsStructure<Set>
{
    fn partial_cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> OrdSignature
    for PentadsStructure<Set>
{
    fn cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        Natural::cmp(
            &self.element_to_enumeration(a),
            &self.element_to_enumeration(b),
        )
    }
}

impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> CountableSetSignature
    for PentadsStructure<Set>
{
    fn generate_all_elements(self: Arc<Self>) -> impl Iterator<Item = Self::Elem> {
        self.list_all_elements_ordered().into_iter()
    }
}

impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> FiniteSetSignature
    for PentadsStructure<Set>
{
    fn size(self: &Arc<Self>) -> Natural {
        Natural::from(6usize)
    }
}

impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> ConstSizeFiniteSetSignature<6>
    for PentadsStructure<Set>
{
}

impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> OrderedFiniteSetSignature
    for PentadsStructure<Set>
{
    fn list_all_elements_ordered(self: &Arc<Self>) -> Vec<Self::Elem> {
        // the ordering here is arbitrary but must be the same every time
        let synthemes_set = self.set().synthemes();

        let p = |i: usize| {
            self.set()
                .enumeration_to_element(&Natural::from(i))
                .unwrap()
        };

        // An arbitrary pentad from which we obtain the others by applying each syntheme permutation
        let root_pentad = Pentad {
            synthemes: [
                Syntheme {
                    duads: [
                        Duad {
                            points: [p(0), p(1)],
                        },
                        Duad {
                            points: [p(3), p(4)],
                        },
                        Duad {
                            points: [p(2), p(5)],
                        },
                    ],
                },
                Syntheme {
                    duads: [
                        Duad {
                            points: [p(0), p(2)],
                        },
                        Duad {
                            points: [p(1), p(3)],
                        },
                        Duad {
                            points: [p(4), p(5)],
                        },
                    ],
                },
                Syntheme {
                    duads: [
                        Duad {
                            points: [p(0), p(3)],
                        },
                        Duad {
                            points: [p(2), p(4)],
                        },
                        Duad {
                            points: [p(1), p(5)],
                        },
                    ],
                },
                Syntheme {
                    duads: [
                        Duad {
                            points: [p(1), p(2)],
                        },
                        Duad {
                            points: [p(0), p(4)],
                        },
                        Duad {
                            points: [p(3), p(5)],
                        },
                    ],
                },
                Syntheme {
                    duads: [
                        Duad {
                            points: [p(2), p(3)],
                        },
                        Duad {
                            points: [p(1), p(4)],
                        },
                        Duad {
                            points: [p(0), p(5)],
                        },
                    ],
                },
            ],
        };
        debug_assert!(self.is_element(&root_pentad));

        let mut pentads = vec![];
        for syntheme in &root_pentad.synthemes {
            pentads.push(
                self.set()
                    .permutations()
                    .pentad_image(&synthemes_set.to_permutation(syntheme), &root_pentad),
            );
        }
        pentads.push(root_pentad);
        debug_assert_eq!(pentads.len(), 6);
        pentads
    }

    fn element_to_enumeration(self: &Arc<Self>, elem: &Self::Elem) -> Natural {
        assert!(self.validate_element(elem).is_ok());
        // found by printing the pentads produced by self.list_all_elements_ordered() and extracting sufficient information to enumerate them
        let x: usize = self
            .set()
            .element_to_enumeration(&elem.synthemes[1].duads[2].points[0])
            .try_into()
            .unwrap();
        let y: usize = self
            .set()
            .element_to_enumeration(&elem.synthemes[3].duads[2].points[0])
            .try_into()
            .unwrap();
        Natural::from(match (x, y) {
            (3, 1) => 0,
            (4, 1) => 1,
            (3, 2) => 2,
            (1, 3) => 3,
            (1, 2) => 4,
            (4, 3) => 5,
            _ => unreachable!(),
        } as usize)
    }

    fn enumeration_to_element(self: &Arc<Self>, num: &Natural) -> Option<Self::Elem> {
        if num < &Natural::from(6usize) {
            let num: usize = num.try_into().unwrap();
            Some(
                self.list_all_elements_ordered()
                    .into_iter()
                    .nth(num)
                    .unwrap(),
            )
        } else {
            None
        }
    }
}

impl<Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature> PentadsStructure<Set> {
    pub fn pentad(
        self: &Arc<Self>,
        synthemes: [Syntheme<Set::Elem>; 5],
    ) -> Result<Pentad<Set::Elem>, &'static str> {
        let synthemes_set = self.set().synthemes();
        for i in 0..5 {
            for j in (i + 1)..5 {
                if !synthemes_set
                    .overlap(&synthemes[i], &synthemes[j])
                    .is_disjoint()
                {
                    return Err("not disjoint");
                }
            }
        }
        let pentad = Pentad {
            synthemes: synthemes_set.sort(synthemes.into()).try_into().unwrap(),
        };
        debug_assert!(self.is_element(&pentad));
        Ok(pentad)
    }
}

#[signature_meta_trait]
pub trait SetPermutationAsPentadPermutation<
    Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature,
>: PermutationsSignature<Set>
{
    fn pentad_image(
        self: &Arc<Self>,
        set_perm: &Self::Elem,
        pentad: &Pentad<Set::Elem>,
    ) -> Pentad<Set::Elem> {
        let set = self.set();
        debug_assert_eq!(set.size(), Natural::from(6usize));
        let pentads = set.pentads();
        pentads
            .pentad(
                pentad
                    .synthemes
                    .iter()
                    .map(|s| self.syntheme_image(set_perm, s))
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap(),
            )
            .unwrap()
    }

    fn pentad_action(
        self: &Arc<Self>,
        set_perm: &Self::Elem,
    ) -> FinitelySupportedPermutation<Pentad<Set::Elem>> {
        let set = self.set();
        debug_assert_eq!(set.size(), Natural::from(6usize));
        let pentads = set.pentads();
        let pentads_perms = pentads.permutations();
        pentads_perms
            .new_perm(
                pentads
                    .list_all_elements()
                    .into_iter()
                    .map(|from| {
                        let to = self.pentad_image(set_perm, &from);
                        (from, to)
                    })
                    .collect(),
            )
            .unwrap()
    }
}
impl<
    Set: ConstSizeFiniteSetSignature<6> + OrderedFiniteSetSignature,
    SetPerms: PermutationsSignature<Set>,
> SetPermutationAsPentadPermutation<Set> for SetPerms
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_structures::MetaType;
    use std::collections::HashMap;

    #[test]
    fn test_enumeration() {
        algebraeon_structures::assert_enumerated_ord_finite_set!(
            i32::structure()
                .const_size_finite_subset([1, 2, 3, 4, 5, 6])
                .pentads(),
            6
        );
    }

    #[test]
    fn test_permutation() {
        let set = i32::structure().const_size_finite_subset([1, 2, 3, 4, 5, 6]);
        let set_perms = set.permutations();
        let pentads = set.pentads();
        let pentad_perms = pentads.permutations();
        assert_eq!(
            pentad_perms
                .cycle_shape(&set_perms.pentad_action(&set_perms.new_cycle(vec![1, 2]).unwrap())),
            HashMap::from([(2, 3)])
        );
    }
}
