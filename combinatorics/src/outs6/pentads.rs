use crate::outs6::*;
use algebraeon_structures::*;
use std::{cmp::Ordering, marker::PhantomData};

#[derive(Debug, Clone)]
pub struct Pentad<Set: EnumeratedOrdFiniteSetSignature> {
    // must have syntheme_1 < syntheme_2 < syntheme_3 < syntheme_4 < syntheme_5 and all disjoint
    syntheme_1: Syntheme<Set>,
    syntheme_2: Syntheme<Set>,
    syntheme_3: Syntheme<Set>,
    syntheme_4: Syntheme<Set>,
    syntheme_5: Syntheme<Set>,
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
    type Elem = Pentad<Set>;

    fn validate_element(&self, p: &Self::Elem) -> Result<(), String> {
        let synthemes = self.set().synthemes().unwrap();
        synthemes.validate_element(&p.syntheme_1)?;
        synthemes.validate_element(&p.syntheme_2)?;
        synthemes.validate_element(&p.syntheme_3)?;
        synthemes.validate_element(&p.syntheme_4)?;
        synthemes.validate_element(&p.syntheme_5)?;
        let p_synthemes = [
            &p.syntheme_1,
            &p.syntheme_2,
            &p.syntheme_3,
            &p.syntheme_4,
            &p.syntheme_5,
        ];
        if !synthemes.is_sorted(&p_synthemes) {
            return Err("synthemes are not sorted".to_string());
        }
        for i in 0..5 {
            for j in (i + 1)..5 {
                if !synthemes
                    .overlap(&p_synthemes[i], &p_synthemes[j])
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
        let synthemes = self.set().synthemes().unwrap();
        synthemes.equal(&a.syntheme_1, &b.syntheme_1)
            && synthemes.equal(&a.syntheme_2, &b.syntheme_2)
            && synthemes.equal(&a.syntheme_3, &b.syntheme_3)
            && synthemes.equal(&a.syntheme_4, &b.syntheme_4)
            && synthemes.equal(&a.syntheme_5, &b.syntheme_5)
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
        let duads_set = self.set().duads().unwrap();
        let synthemes_set = self.set().synthemes().unwrap();

        let p = |i: usize| {
            self.set()
                .enumeration_to_element(&Natural::from(i))
                .unwrap()
        };

        let mut pentads = vec![];
        todo!();

        debug_assert_eq!(pentads.len(), 6);
        pentads
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

#[cfg(test)]
mod tests {
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

        // pentads are all valid
        for s in &pentads {
            println!("{:?}", s);
            assert!(pentads_set.validate_element(s).is_ok());
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
}
