use crate::sets::SetToFiniteSubsetsByOrdSignature;
use algebraeon_structures::*;
use std::cmp::Ordering;
use std::marker::PhantomData;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSubsetByOrd<Set: SetSignature> {
    // ordered
    pub elems: Vec<Set::Elem>,
}

impl<Set: OrdSignature> FiniteSubsetByOrd<Set> {
    pub fn size(&self) -> usize {
        self.elems.len()
    }
}

// A finite subset of a set
#[derive(Debug, Clone)]
pub struct FiniteSubsetByOrdStructure<Set: OrdSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
    // ordered
    subset: FiniteSubsetByOrd<Set>,
}

pub trait SetToFiniteSubsetByOrdSignature: OrdSignature {
    fn finite_subset(&self, elems: Vec<Self::Elem>) -> FiniteSubsetByOrdStructure<Self, &Self> {
        FiniteSubsetByOrdStructure::new(self, elems)
    }

    fn into_finite_subset(self, elems: Vec<Self::Elem>) -> FiniteSubsetByOrdStructure<Self, Self> {
        FiniteSubsetByOrdStructure::new(self, elems)
    }
}
impl<Set: OrdSignature> SetToFiniteSubsetByOrdSignature for Set {}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> PartialEq
    for FiniteSubsetByOrdStructure<Set, SetB>
{
    fn eq(&self, other: &Self) -> bool {
        let set = self.set();
        if set != other.set() {
            return false;
        }
        set.finite_subsets().equal(&self.subset, &other.subset)
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> Eq for FiniteSubsetByOrdStructure<Set, SetB> {}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> FiniteSubsetByOrdStructure<Set, SetB> {
    pub fn new(set: SetB, elems: Vec<Set::Elem>) -> Self {
        debug_assert!(set.borrow().is_sorted_and_unique(&elems));
        Self {
            _set: PhantomData,
            set,
            subset: FiniteSubsetByOrd { elems },
        }
    }

    pub fn set(&self) -> &Set {
        self.set.borrow()
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> Signature
    for FiniteSubsetByOrdStructure<Set, SetB>
{
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> SetSignature
    for FiniteSubsetByOrdStructure<Set, SetB>
{
    type Elem = Set::Elem;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if !self.set().binary_search(&self.subset.elems, x) {
            return Err("element not in finite subset".to_string());
        }
        Ok(())
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> EqSignature
    for FiniteSubsetByOrdStructure<Set, SetB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().equal(a, b)
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> PartialOrdSignature
    for FiniteSubsetByOrdStructure<Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().partial_cmp(a, b)
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> OrdSignature
    for FiniteSubsetByOrdStructure<Set, SetB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().cmp(a, b)
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> CountableSetSignature
    for FiniteSubsetByOrdStructure<Set, SetB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.subset.elems.into_iter()
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> FiniteSetSignature
    for FiniteSubsetByOrdStructure<Set, SetB>
{
    fn size(&self) -> Natural {
        Natural::from(self.subset.size())
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>
    EnumeratedOrdFiniteSetSignature for FiniteSubsetByOrdStructure<Set, SetB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.list_all_elements()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        #[cfg(debug_assertions)]
        self.validate_element(elem).unwrap();
        Natural::from(
            self.set()
                .binary_search_index(&self.subset.elems, elem)
                .unwrap(),
        )
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        let num: Result<usize, ()> = num.try_into();
        if let Ok(num) = num {
            self.subset.elems.get(num).cloned()
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::sets::SetToFiniteSubsetByOrdSignature;
    use algebraeon_structures::*;

    #[test]
    fn test() {
        let set = i32::structure().into_finite_subset(vec![1, 2, 3, 4, 5, 6]);
        assert_eq!(set.size(), Natural::from(6usize));
    }
}
