use crate::outs6::*;
use algebraeon_macros::signature_meta_trait;
use algebraeon_sets::sets::{
    FiniteSetToFinitelySupportedPermutationsStructure, FinitelySupportedPermutation,
    SetToFiniteSubsetsByOrdSignature,
};
use algebraeon_structures::*;
use std::{cmp::Ordering, marker::PhantomData};

#[derive(Debug, Clone)]
pub struct Syntheme<Elem> {
    // must have duad_1 < duad_2 < duad_3 and all disjoint
    pub(crate) duads: [Duad<Elem>; 3],
}

/// The 15-element set of duads on a 6-element set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SynthemesStructure<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> SynthemesStructure<Set, SetB>
{
    pub fn new(set: SetB) -> Self {
        debug_assert_eq!(set.borrow().size(), Natural::from(6usize));
        Self {
            _set: PhantomData,
            set,
        }
    }

    pub fn set(&self) -> &Set {
        self.set.borrow()
    }
}

pub trait SetToSynthemesSignature:
    FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature
{
    fn synthemes(&self) -> SynthemesStructure<Self, &Self> {
        SynthemesStructure::new(self)
    }

    fn into_synthemes(self) -> SynthemesStructure<Self, Self> {
        SynthemesStructure::new(self)
    }
}
impl<Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature> SetToSynthemesSignature
    for Set
{
}

impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> Signature for SynthemesStructure<Set, SetB>
{
}

impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> SetSignature for SynthemesStructure<Set, SetB>
{
    type Elem = Syntheme<Set::Elem>;

    fn validate_element(&self, s: &Self::Elem) -> Result<(), String> {
        let duads = self.set().duads();
        for duad in &s.duads {
            duads.validate_element(duad)?;
        }
        if !duads.is_sorted(&s.duads) {
            return Err("duads are not sorted".to_string());
        }
        if !self
            .set()
            .finite_subsets()
            .is_disjoint(&(&s.duads[0]).into(), &(&s.duads[1]).into())
            || !self
                .set()
                .finite_subsets()
                .is_disjoint(&(&s.duads[0]).into(), &(&s.duads[2]).into())
            || !self
                .set()
                .finite_subsets()
                .is_disjoint(&(&s.duads[1]).into(), &(&s.duads[2]).into())
        {
            return Err("duads are not disjoint".to_string());
        }
        Ok(())
    }
}

impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> EqSignature for SynthemesStructure<Set, SetB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        let duads = self.set().duads();
        (0..3).all(|i| duads.equal(&a.duads[i], &b.duads[i]))
    }
}

impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> PartialOrdSignature for SynthemesStructure<Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> OrdSignature for SynthemesStructure<Set, SetB>
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

impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> CountableSetSignature for SynthemesStructure<Set, SetB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        (0usize..15).map(move |i| self.enumeration_to_element(&Natural::from(i)).unwrap())
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> FiniteSetSignature for SynthemesStructure<Set, SetB>
{
    fn size(&self) -> Natural {
        Natural::from(15usize)
    }
}

impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> FiniteSetSizedSignature<15> for SynthemesStructure<Set, SetB>
{
}

impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> EnumeratedOrdFiniteSetSignature for SynthemesStructure<Set, SetB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.list_all_elements()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        assert!(self.validate_element(elem).is_ok());
        for (i, s) in self.generate_all_elements().enumerate() {
            if self.equal(&s, elem) {
                return Natural::from(i);
            }
        }
        unreachable!()
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        let p = |i: usize| {
            self.set()
                .enumeration_to_element(&Natural::from(i))
                .unwrap()
        };
        if *num == Natural::from(0usize) {
            Some(Syntheme {
                duads: [
                    Duad {
                        points: [p(0), p(1)],
                    },
                    Duad {
                        points: [p(2), p(3)],
                    },
                    Duad {
                        points: [p(4), p(5)],
                    },
                ],
            })
        } else if *num == Natural::from(1usize) {
            Some(Syntheme {
                duads: [
                    Duad {
                        points: [p(0), p(1)],
                    },
                    Duad {
                        points: [p(2), p(4)],
                    },
                    Duad {
                        points: [p(3), p(5)],
                    },
                ],
            })
        } else if *num == Natural::from(2usize) {
            Some(Syntheme {
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
            })
        } else if *num == Natural::from(3usize) {
            Some(Syntheme {
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
            })
        } else if *num == Natural::from(4usize) {
            Some(Syntheme {
                duads: [
                    Duad {
                        points: [p(0), p(2)],
                    },
                    Duad {
                        points: [p(1), p(4)],
                    },
                    Duad {
                        points: [p(3), p(5)],
                    },
                ],
            })
        } else if *num == Natural::from(5usize) {
            Some(Syntheme {
                duads: [
                    Duad {
                        points: [p(0), p(2)],
                    },
                    Duad {
                        points: [p(3), p(4)],
                    },
                    Duad {
                        points: [p(1), p(5)],
                    },
                ],
            })
        } else if *num == Natural::from(6usize) {
            Some(Syntheme {
                duads: [
                    Duad {
                        points: [p(1), p(2)],
                    },
                    Duad {
                        points: [p(0), p(3)],
                    },
                    Duad {
                        points: [p(4), p(5)],
                    },
                ],
            })
        } else if *num == Natural::from(7usize) {
            Some(Syntheme {
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
            })
        } else if *num == Natural::from(8usize) {
            Some(Syntheme {
                duads: [
                    Duad {
                        points: [p(0), p(3)],
                    },
                    Duad {
                        points: [p(1), p(4)],
                    },
                    Duad {
                        points: [p(2), p(5)],
                    },
                ],
            })
        } else if *num == Natural::from(9usize) {
            Some(Syntheme {
                duads: [
                    Duad {
                        points: [p(2), p(3)],
                    },
                    Duad {
                        points: [p(0), p(4)],
                    },
                    Duad {
                        points: [p(1), p(5)],
                    },
                ],
            })
        } else if *num == Natural::from(10usize) {
            Some(Syntheme {
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
            })
        } else if *num == Natural::from(11usize) {
            Some(Syntheme {
                duads: [
                    Duad {
                        points: [p(1), p(3)],
                    },
                    Duad {
                        points: [p(0), p(4)],
                    },
                    Duad {
                        points: [p(2), p(5)],
                    },
                ],
            })
        } else if *num == Natural::from(12usize) {
            Some(Syntheme {
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
            })
        } else if *num == Natural::from(13usize) {
            Some(Syntheme {
                duads: [
                    Duad {
                        points: [p(1), p(3)],
                    },
                    Duad {
                        points: [p(2), p(4)],
                    },
                    Duad {
                        points: [p(0), p(5)],
                    },
                ],
            })
        } else if *num == Natural::from(14usize) {
            Some(Syntheme {
                duads: [
                    Duad {
                        points: [p(1), p(2)],
                    },
                    Duad {
                        points: [p(3), p(4)],
                    },
                    Duad {
                        points: [p(0), p(5)],
                    },
                ],
            })
        } else {
            None
        }
    }
}

pub enum SynthemeOverlapResult<Set: EnumeratedOrdFiniteSetSignature> {
    Equal,
    Disjoint,
    UniqueCommonDuad(Duad<Set::Elem>),
}

impl<Set: EnumeratedOrdFiniteSetSignature> SynthemeOverlapResult<Set> {
    pub fn is_equal(&self) -> bool {
        match self {
            SynthemeOverlapResult::Equal => true,
            SynthemeOverlapResult::Disjoint | SynthemeOverlapResult::UniqueCommonDuad(_) => false,
        }
    }

    pub fn is_disjoint(&self) -> bool {
        match self {
            SynthemeOverlapResult::Disjoint => true,
            SynthemeOverlapResult::Equal | SynthemeOverlapResult::UniqueCommonDuad(_) => false,
        }
    }

    pub fn unwrap_common_duad(&self) -> &Duad<Set::Elem> {
        match self {
            SynthemeOverlapResult::Equal | SynthemeOverlapResult::Disjoint => panic!(),
            SynthemeOverlapResult::UniqueCommonDuad(duad) => duad,
        }
    }

    pub fn into_unwrap_common_duad(self) -> Duad<Set::Elem> {
        match self {
            SynthemeOverlapResult::Equal | SynthemeOverlapResult::Disjoint => panic!(),
            SynthemeOverlapResult::UniqueCommonDuad(duad) => duad,
        }
    }
}

impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
> SynthemesStructure<Set, SetB>
{
    pub fn syntheme(
        &self,
        duads: [Duad<Set::Elem>; 3],
    ) -> Result<Syntheme<Set::Elem>, &'static str> {
        let duads_set = self.set().duads();
        let sorted_duads: [_; 3] = duads_set.sort(duads.into()).try_into().unwrap();
        if duads_set
            .overlap(&sorted_duads[0], &sorted_duads[1])
            .is_disjoint()
            && duads_set
                .overlap(&sorted_duads[0], &sorted_duads[2])
                .is_disjoint()
            && duads_set
                .overlap(&sorted_duads[1], &sorted_duads[2])
                .is_disjoint()
        {
            let syntheme = Syntheme {
                duads: sorted_duads,
            };
            debug_assert!(self.is_element(&syntheme));
            Ok(syntheme)
        } else {
            Err("duads are not disjoint")
        }
    }

    pub fn overlap(
        &self,
        s1: &Syntheme<Set::Elem>,
        s2: &Syntheme<Set::Elem>,
    ) -> SynthemeOverlapResult<Set> {
        let duads_set = self.set().duads();
        let mut common = vec![];
        for item in
            duads_set.merge_sorted_and_unique(s1.duads.iter().collect(), s2.duads.iter().collect())
        {
            match item {
                MergedUniqueSource::First(_) | MergedUniqueSource::Second(_) => {}
                MergedUniqueSource::Both(d1, d2) => {
                    debug_assert!(duads_set.equal(d1, d2));
                    common.push(d1);
                }
            }
        }
        match common.len() {
            0 => SynthemeOverlapResult::Disjoint,
            1 => {
                SynthemeOverlapResult::UniqueCommonDuad(common.into_iter().next().unwrap().clone())
            }
            3 => SynthemeOverlapResult::Equal,
            _ => {
                unreachable!()
            }
        }
    }

    /// Convert a syntheme into a 2^3-cycle
    pub fn to_permutation(
        &self,
        syntheme: &Syntheme<Set::Elem>,
    ) -> FinitelySupportedPermutation<Set::Elem> {
        let p1 = &syntheme.duads[0].points[0];
        let p2 = &syntheme.duads[0].points[1];
        let p3 = &syntheme.duads[1].points[0];
        let p4 = &syntheme.duads[1].points[1];
        let p5 = &syntheme.duads[2].points[0];
        let p6 = &syntheme.duads[2].points[1];
        self.set()
            .permutations()
            .new_perm(vec![
                (p1.clone(), p2.clone()),
                (p2.clone(), p1.clone()),
                (p3.clone(), p4.clone()),
                (p4.clone(), p3.clone()),
                (p5.clone(), p6.clone()),
                (p6.clone(), p5.clone()),
            ])
            .unwrap()
    }

    /// Convert a 2^3-cycle into a syntheme, or return None if the permutation is not a 2^3-cycle
    pub fn try_from_permutation(
        &self,
        swap: &FinitelySupportedPermutation<Set::Elem>,
    ) -> Option<Syntheme<Set::Elem>> {
        let disjoint_cycles = self.set().permutations().disjoint_cycles(swap);
        if disjoint_cycles.len() != 3 {
            return None;
        }
        for cycle in &disjoint_cycles {
            if cycle.len() != 2 {
                return None;
            }
        }
        let duads = self.set().duads();
        Some(
            self.syntheme(
                disjoint_cycles
                    .into_iter()
                    .map(|cycle| {
                        let mut cycle = cycle.into_iter();
                        duads
                            .duad(cycle.next().unwrap(), cycle.next().unwrap())
                            .unwrap()
                    })
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap(),
            )
            .unwrap(),
        )
    }
}

#[signature_meta_trait]
pub trait SetPermutationAsSynthemePermutation<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
>: PermutationsSignature<Set>
{
    fn syntheme_image(
        &self,
        set_perm: &Self::Elem,
        syntheme: &Syntheme<Set::Elem>,
    ) -> Syntheme<Set::Elem> {
        let set = self.set();
        debug_assert_eq!(set.size(), Natural::from(6usize));
        let synthemes = set.synthemes();
        synthemes
            .syntheme([
                self.duad_image(set_perm, &syntheme.duads[0]),
                self.duad_image(set_perm, &syntheme.duads[1]),
                self.duad_image(set_perm, &syntheme.duads[2]),
            ])
            .unwrap()
    }

    fn syntheme_action(
        &self,
        set_perm: &Self::Elem,
    ) -> FinitelySupportedPermutation<Syntheme<Set::Elem>> {
        let set = self.set();
        debug_assert_eq!(set.size(), Natural::from(6usize));
        let synthemes = set.synthemes();
        let syntheme_perms = synthemes.permutations();
        syntheme_perms
            .new_perm(
                synthemes
                    .list_all_elements()
                    .into_iter()
                    .map(|from| {
                        let to = self.syntheme_image(set_perm, &from);
                        (from, to)
                    })
                    .collect(),
            )
            .unwrap()
    }
}
impl<
    Set: FiniteSetSizedSignature<6> + EnumeratedOrdFiniteSetSignature,
    SetPerms: PermutationsSignature<Set>,
> SetPermutationAsSynthemePermutation<Set> for SetPerms
{
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::*;
    use algebraeon_sets::sets::SetToFiniteSubsetByOrdSizedSignature;
    use algebraeon_structures::MetaType;

    #[test]
    fn test_enumeration() {
        let set = i32::structure().into_finite_subset_sized([1, 2, 3, 4, 5, 6]);
        let synthemes_set = set.synthemes();
        let synthemes = synthemes_set.list_all_elements_ordered();
        assert_eq!(synthemes.len(), 15);
        assert_eq!(synthemes_set.size(), Natural::from(15usize));

        // synthemes are all valid
        for s in &synthemes {
            println!("{:?}", s);
            assert!(synthemes_set.validate_element(s).is_ok());
        }

        // synthemes are all disjoint
        for i in 0..15 {
            for j in (i + 1)..15 {
                let si = &synthemes[i];
                let sj = &synthemes[j];
                assert!(synthemes_set.cmp(si, sj).is_lt());
            }
        }

        // enumeration is correct
        for (i, s) in synthemes.iter().enumerate() {
            assert_eq!(Natural::from(i), synthemes_set.element_to_enumeration(s));
            assert!(
                synthemes_set.equal(
                    &synthemes_set
                        .enumeration_to_element(&Natural::from(i))
                        .unwrap(),
                    s
                )
            );
        }
        assert!(
            synthemes_set
                .enumeration_to_element(&Natural::from(15usize))
                .is_none()
        );
    }

    #[test]
    fn test_overlap() {
        let set = i32::structure().into_finite_subset_sized([1, 2, 3, 4, 5, 6]);
        let duads_set = set.duads();
        let synthemes_set = set.synthemes();

        assert!(
            synthemes_set
                .overlap(
                    &synthemes_set
                        .syntheme([
                            duads_set.duad(1, 2).unwrap(),
                            duads_set.duad(3, 4).unwrap(),
                            duads_set.duad(5, 6).unwrap(),
                        ])
                        .unwrap(),
                    &synthemes_set
                        .syntheme([
                            duads_set.duad(1, 2).unwrap(),
                            duads_set.duad(3, 4).unwrap(),
                            duads_set.duad(5, 6).unwrap(),
                        ])
                        .unwrap()
                )
                .is_equal()
        );
    }

    #[test]
    fn test_permutation() {
        let set = i32::structure().into_finite_subset_sized([1, 2, 3, 4, 5, 6]);
        let set_perms = set.permutations();
        let duads = set.duads();
        let synthemes = set.synthemes();
        let syntheme_perms = synthemes.permutations();

        assert_eq!(
            syntheme_perms.cycle_shape(
                &set_perms.syntheme_action(&set_perms.new_cycle(vec![1, 2, 3]).unwrap())
            ),
            HashMap::from([(3, 5)])
        );

        assert!(
            synthemes.equal(
                &synthemes
                    .try_from_permutation(
                        &set_perms
                            .new_perm(vec![(1, 2), (2, 1), (3, 4), (4, 3), (5, 6), (6, 5)])
                            .unwrap()
                    )
                    .unwrap(),
                &synthemes
                    .syntheme([
                        duads.duad(1, 2).unwrap(),
                        duads.duad(3, 4).unwrap(),
                        duads.duad(5, 6).unwrap()
                    ])
                    .unwrap()
            )
        );

        assert!(
            &synthemes
                .try_from_permutation(&set_perms.new_cycle(vec![2, 3, 4]).unwrap())
                .is_none()
        );
    }
}
