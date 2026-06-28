use crate::outs6::*;
use algebraeon_sets::sets::SetToFiniteSubsetsByOrdSignature;
use algebraeon_structures::*;
use std::{cmp::Ordering, marker::PhantomData};

#[derive(Debug, Clone)]
pub struct Syntheme<Set: EnumeratedOrdFiniteSetSignature> {
    // must have duad_1 < duad_2 < duad_3 and all disjoint
    pub duad_1: Duad<Set>,
    pub duad_2: Duad<Set>,
    pub duad_3: Duad<Set>,
}

/// The 15-element set of duads on a 6-element set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SynthemesStructure<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>
    SynthemesStructure<Set, SetB>
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

pub trait SetToSynthemesSignature: EnumeratedOrdFiniteSetSignature {
    fn synthemes(&self) -> Option<SynthemesStructure<Self, &Self>> {
        SynthemesStructure::try_new(self)
    }

    fn into_synthemes(self) -> Option<SynthemesStructure<Self, Self>> {
        SynthemesStructure::try_new(self)
    }
}
impl<Set: EnumeratedOrdFiniteSetSignature> SetToSynthemesSignature for Set {}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> Signature
    for SynthemesStructure<Set, SetB>
{
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> SetSignature
    for SynthemesStructure<Set, SetB>
{
    type Elem = Syntheme<Set>;

    fn validate_element(&self, s: &Self::Elem) -> Result<(), String> {
        let duads = self.set().duads().unwrap();
        duads.validate_element(&s.duad_1)?;
        duads.validate_element(&s.duad_2)?;
        duads.validate_element(&s.duad_3)?;
        if !duads.is_sorted(&[&s.duad_1, &s.duad_2, &s.duad_3]) {
            return Err("duads are not sorted".to_string());
        }
        if !self
            .set()
            .finite_subsets()
            .is_disjoint(&(&s.duad_1).into(), &(&s.duad_2).into())
            || !self
                .set()
                .finite_subsets()
                .is_disjoint(&(&s.duad_1).into(), &(&s.duad_3).into())
            || !self
                .set()
                .finite_subsets()
                .is_disjoint(&(&s.duad_2).into(), &(&s.duad_3).into())
        {
            return Err("duads are not disjoint".to_string());
        }
        Ok(())
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> EqSignature
    for SynthemesStructure<Set, SetB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        let duads = self.set().duads().unwrap();
        duads.equal(&a.duad_1, &b.duad_1)
            && duads.equal(&a.duad_2, &b.duad_2)
            && duads.equal(&a.duad_3, &b.duad_3)
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> PartialOrdSignature
    for SynthemesStructure<Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> OrdSignature
    for SynthemesStructure<Set, SetB>
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
    for SynthemesStructure<Set, SetB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        (0usize..15).map(move |i| self.enumeration_to_element(&Natural::from(i)).unwrap())
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> FiniteSetSignature
    for SynthemesStructure<Set, SetB>
{
    fn size(&self) -> Natural {
        Natural::from(15usize)
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>
    EnumeratedOrdFiniteSetSignature for SynthemesStructure<Set, SetB>
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
                duad_1: Duad { p1: p(0), p2: p(1) },
                duad_2: Duad { p1: p(2), p2: p(3) },
                duad_3: Duad { p1: p(4), p2: p(5) },
            })
        } else if *num == Natural::from(1usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(0), p2: p(1) },
                duad_2: Duad { p1: p(2), p2: p(4) },
                duad_3: Duad { p1: p(3), p2: p(5) },
            })
        } else if *num == Natural::from(2usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(0), p2: p(1) },
                duad_2: Duad { p1: p(3), p2: p(4) },
                duad_3: Duad { p1: p(2), p2: p(5) },
            })
        } else if *num == Natural::from(3usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(0), p2: p(2) },
                duad_2: Duad { p1: p(1), p2: p(3) },
                duad_3: Duad { p1: p(4), p2: p(5) },
            })
        } else if *num == Natural::from(4usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(0), p2: p(2) },
                duad_2: Duad { p1: p(1), p2: p(4) },
                duad_3: Duad { p1: p(3), p2: p(5) },
            })
        } else if *num == Natural::from(5usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(0), p2: p(2) },
                duad_2: Duad { p1: p(3), p2: p(4) },
                duad_3: Duad { p1: p(1), p2: p(5) },
            })
        } else if *num == Natural::from(6usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(1), p2: p(2) },
                duad_2: Duad { p1: p(0), p2: p(3) },
                duad_3: Duad { p1: p(4), p2: p(5) },
            })
        } else if *num == Natural::from(7usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(0), p2: p(3) },
                duad_2: Duad { p1: p(2), p2: p(4) },
                duad_3: Duad { p1: p(1), p2: p(5) },
            })
        } else if *num == Natural::from(8usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(0), p2: p(3) },
                duad_2: Duad { p1: p(1), p2: p(4) },
                duad_3: Duad { p1: p(2), p2: p(5) },
            })
        } else if *num == Natural::from(9usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(2), p2: p(3) },
                duad_2: Duad { p1: p(0), p2: p(4) },
                duad_3: Duad { p1: p(1), p2: p(5) },
            })
        } else if *num == Natural::from(10usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(1), p2: p(2) },
                duad_2: Duad { p1: p(0), p2: p(4) },
                duad_3: Duad { p1: p(3), p2: p(5) },
            })
        } else if *num == Natural::from(11usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(1), p2: p(3) },
                duad_2: Duad { p1: p(0), p2: p(4) },
                duad_3: Duad { p1: p(2), p2: p(5) },
            })
        } else if *num == Natural::from(12usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(2), p2: p(3) },
                duad_2: Duad { p1: p(1), p2: p(4) },
                duad_3: Duad { p1: p(0), p2: p(5) },
            })
        } else if *num == Natural::from(13usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(1), p2: p(3) },
                duad_2: Duad { p1: p(2), p2: p(4) },
                duad_3: Duad { p1: p(0), p2: p(5) },
            })
        } else if *num == Natural::from(14usize) {
            Some(Syntheme {
                duad_1: Duad { p1: p(1), p2: p(2) },
                duad_2: Duad { p1: p(3), p2: p(4) },
                duad_3: Duad { p1: p(0), p2: p(5) },
            })
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_sets::sets::SetToFiniteSubsetByOrdSignature;
    use algebraeon_structures::MetaType;

    #[test]
    fn test_enumeration() {
        let set = i32::structure().into_finite_subset(vec![1, 2, 3, 4, 5, 6]);
        for p in set.generate_all_elements() {
            println!("{:?}", p);
        }

        let synthemes_set = set.synthemes().unwrap();
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
}
