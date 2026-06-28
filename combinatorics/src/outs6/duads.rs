//! Given a 6-element set make the following definitions
//!  - A duad is a 2-element subset
//!  - A syntheme is an unordered set of 3 disjoint duads
//!  - A pentad is an unordered set of 5 disjoint synthemes
//!
//! Then there are 6 pentads, and any permutation of the 6 points induces a permutation of the 6 pentads via an outer automorphism of S6

use algebraeon_sets::sets::{
    FiniteSubsetByOrd, SetToFiniteSubsetsByOrdSignature, SetToFixedSizeFiniteSubsetsByOrdSignature,
};
use algebraeon_structures::*;
use std::{cmp::Ordering, marker::PhantomData};

/// A 2-element subset of a 6-element set
#[derive(Debug, Clone)]
pub struct Duad<Set: EnumeratedOrdFiniteSetSignature> {
    // must have p1 < p2
    p1: Set::Elem,
    p2: Set::Elem,
}

impl<Set: EnumeratedOrdFiniteSetSignature> TryFrom<FiniteSubsetByOrd<Set>> for Duad<Set> {
    type Error = &'static str;

    fn try_from(subset: FiniteSubsetByOrd<Set>) -> Result<Self, Self::Error> {
        let mut elems = subset.elems.into_iter();
        let duad = Self {
            p1: elems.next().ok_or("foo")?,
            p2: elems.next().ok_or("foo")?,
        };
        if elems.next().is_some() {
            return Err("foobar");
        }
        Ok(duad)
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature> From<Duad<Set>> for FiniteSubsetByOrd<Set> {
    fn from(val: Duad<Set>) -> Self {
        FiniteSubsetByOrd {
            elems: vec![val.p1, val.p2],
        }
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature> From<&Duad<Set>> for FiniteSubsetByOrd<Set> {
    fn from(val: &Duad<Set>) -> Self {
        FiniteSubsetByOrd {
            elems: vec![val.p1.clone(), val.p2.clone()],
        }
    }
}

/// The 15-element set of duads on a 6-element set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DuadsStructure<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> DuadsStructure<Set, SetB> {
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

pub trait SetToDuadsSignature: EnumeratedOrdFiniteSetSignature {
    fn duads(&self) -> Option<DuadsStructure<Self, &Self>> {
        DuadsStructure::try_new(self)
    }

    fn into_duads(self) -> Option<DuadsStructure<Self, Self>> {
        DuadsStructure::try_new(self)
    }
}
impl<Set: EnumeratedOrdFiniteSetSignature> SetToDuadsSignature for Set {}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> Signature
    for DuadsStructure<Set, SetB>
{
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> SetSignature
    for DuadsStructure<Set, SetB>
{
    type Elem = Duad<Set>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if !self.set().cmp(&x.p1, &x.p2).is_lt() {
            return Err("invalid duad".to_string());
        }
        Ok(())
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> EqSignature
    for DuadsStructure<Set, SetB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.set().finite_subsets().equal(&a.into(), &b.into())
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> PartialOrdSignature
    for DuadsStructure<Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        self.set()
            .finite_subsets()
            .partial_cmp(&a.into(), &b.into())
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> OrdSignature
    for DuadsStructure<Set, SetB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        self.set().finite_subsets().cmp(&a.into(), &b.into())
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> CountableSetSignature
    for DuadsStructure<Set, SetB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.set()
            .clone()
            .into_fixed_size_finite_subsets(2)
            .into_generate_all_elements()
            .map(|subset| Duad::try_from(subset).unwrap())
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>> FiniteSetSignature
    for DuadsStructure<Set, SetB>
{
    fn size(&self) -> Natural {
        debug_assert_eq!(
            self.set().fixed_size_finite_subsets(2).size(),
            Natural::from(15usize)
        );
        Natural::from(15usize)
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>
    EnumeratedOrdFiniteSetSignature for DuadsStructure<Set, SetB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.list_all_elements()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        self.set()
            .fixed_size_finite_subsets(2)
            .element_to_enumeration(&elem.into())
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        self.set()
            .fixed_size_finite_subsets(2)
            .enumeration_to_element(num)
            .map(|subset| Duad::try_from(subset).unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_sets::sets::SetToFiniteSubsetByOrdSignature;
    use algebraeon_structures::MetaType;

    #[test]
    fn test() {
        let set = i32::structure().into_finite_subset(vec![1, 2, 3, 4, 5, 6]);

        for p in set.generate_all_elements() {
            println!("{:?}", p);
        }

        let duads = set.duads().unwrap();

        for d in duads.generate_all_elements() {
            println!("{:?}", d);
        }

        println!(
            "{:?}",
            duads.enumeration_to_element(&Natural::from(15usize))
        );
    }
}
