use crate::outs6::*;
use algebraeon_macros::signature_meta_trait;
use algebraeon_sets::sets::{
    FiniteSetToFinitelySupportedPermutationsStructure, FinitelySupportedPermutation,
};
use algebraeon_structures::*;
use std::{cmp::Ordering, marker::PhantomData};

#[derive(Debug, Clone)]
pub struct Pentad<Elem> {
    // must have syntheme_1 < syntheme_2 < syntheme_3 < syntheme_4 < syntheme_5 and all disjoint
    pub(crate) synthemes: [Syntheme<Elem>; 5],
}

/// The 15-element set of duads on a 6-element set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PentadsStructure<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>
    PentadsStructure<Set, SetB>
{
    pub fn try_new(set: SetB) -> Option<Self> {
        if set.borrow().size() == Natural::from(6usize) {
            Some(Self {
                _set: PhantomData,
                set,
            })
        } else {
            None
        }
    }

    pub fn set(&self) -> &Set {
        self.set.borrow()
    }
}

pub trait SetToPentadsSignature: EnumeratedOrdFiniteSetSignature {
    fn pentads(&self) -> Option<PentadsStructure<Self, &Self>> {
        PentadsStructure::try_new(self)
    }

    fn into_pentads(self) -> Option<PentadsStructure<Self, Self>> {
        PentadsStructure::try_new(self)
    }
}
impl<Set: EnumeratedOrdFiniteSetSignature> SetToPentadsSignature for Set {}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> Signature
    for PentadsStructure<Set, SetB>
{
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> SetSignature
    for PentadsStructure<Set, SetB>
{
    type Elem = Pentad<Set::Elem>;

    fn validate_element(&self, p: &Self::Elem) -> Result<(), String> {
        let synthemes = self.set().synthemes().unwrap();
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

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> EqSignature
    for PentadsStructure<Set, SetB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.cmp(a, b).is_eq()
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> PartialOrdSignature
    for PentadsStructure<Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> OrdSignature
    for PentadsStructure<Set, SetB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        Natural::cmp(
            &self.element_to_enumeration(a),
            &self.element_to_enumeration(b),
        )
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> CountableSetSignature
    for PentadsStructure<Set, SetB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.list_all_elements_ordered().into_iter()
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> FiniteSetSignature
    for PentadsStructure<Set, SetB>
{
    fn size(&self) -> Natural {
        Natural::from(6usize)
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>
    EnumeratedOrdFiniteSetSignature for PentadsStructure<Set, SetB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        // the ordering here is arbitrary but must be the same every time
        let synthemes_set = self.set().synthemes().unwrap();

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

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
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

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
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

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>
    PentadsStructure<Set, SetB>
{
    pub fn pentad(
        &self,
        synthemes: [Syntheme<Set::Elem>; 5],
    ) -> Result<Pentad<Set::Elem>, &'static str> {
        let synthemes_set = self.set().synthemes().unwrap();
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
pub trait SetPermutationAsPentadPermutation<Set: EnumeratedOrdFiniteSetSignature>:
    PermutationsSignature<Set>
{
    fn pentad_image(&self, set_perm: &Self::Elem, pentad: &Pentad<Set::Elem>) -> Pentad<Set::Elem> {
        let set = self.set();
        debug_assert_eq!(set.size(), Natural::from(6usize));
        let pentads = set.pentads().unwrap();
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
        &self,
        set_perm: &Self::Elem,
    ) -> FinitelySupportedPermutation<Pentad<Set::Elem>> {
        let set = self.set();
        debug_assert_eq!(set.size(), Natural::from(6usize));
        let pentads = set.pentads().unwrap();
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
impl<Set: EnumeratedOrdFiniteSetSignature, SetPerms: PermutationsSignature<Set>>
    SetPermutationAsPentadPermutation<Set> for SetPerms
{
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use algebraeon_sets::sets::SetToFiniteSubsetByOrdSignature;
    use algebraeon_structures::MetaType;

    #[test]
    fn test_enumeration() {
        let set = i32::structure().into_finite_subset(vec![1, 2, 3, 4, 5, 6]);
        let pentads_set = set.pentads().unwrap();
        let pentads = pentads_set.list_all_elements_ordered();
        assert_eq!(pentads.len(), 6);
        assert_eq!(pentads_set.size(), Natural::from(6usize));

        debug_assert_eq!(pentads.len(), 6);
        // pentads are all valid
        for p in &pentads {
            println!("{:?}", p);
            assert!(pentads_set.validate_element(p).is_ok());
        }

        // synthemes are all distinct
        for i in 0..6 {
            for j in (i + 1)..6 {
                let si = &pentads[i];
                let sj = &pentads[j];
                assert!(pentads_set.cmp(si, sj).is_lt());
            }
        }

        // enumeration is correct
        for (i, s) in pentads.iter().enumerate() {
            assert_eq!(Natural::from(i), pentads_set.element_to_enumeration(s));
            assert!(
                pentads_set.equal(
                    &pentads_set
                        .enumeration_to_element(&Natural::from(i))
                        .unwrap(),
                    s
                )
            );
        }
        assert!(
            pentads_set
                .enumeration_to_element(&Natural::from(6usize))
                .is_none()
        );
    }

    #[test]
    fn test_permutation() {
        let set = i32::structure().into_finite_subset(vec![1, 2, 3, 4, 5, 6]);
        let set_perms = set.permutations();
        let pentads = set.pentads().unwrap();
        let pentad_perms = pentads.permutations();
        assert_eq!(
            pentad_perms
                .cycle_shape(&set_perms.pentad_action(&set_perms.new_cycle(vec![1, 2]).unwrap())),
            HashMap::from([(2, 3)])
        );
    }
}
