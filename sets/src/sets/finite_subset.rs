use algebraeon_structures::*;
use std::cmp::Ordering;
use std::collections::HashSet;
use std::hash::Hash;
use std::marker::PhantomData;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSubsetStructure<Set: SetSignature, SetB: BorrowedStructure<Set>>
where
    Set::Elem: MetaType + Eq + Hash,
{
    _set: PhantomData<Set>,
    set: SetB,
    elems: HashSet<Set::Elem>,
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> FiniteSubsetStructure<Set, SetB>
where
    Set::Elem: MetaType + Eq + Hash,
{
    pub fn new(set: SetB, elems: HashSet<Set::Elem>) -> Self {
        Self {
            _set: PhantomData,
            set,
            elems,
        }
    }

    pub fn set(&self) -> &Set {
        self.set.borrow()
    }
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> Signature for FiniteSubsetStructure<Set, SetB> where
    Set::Elem: MetaType + Eq + Hash
{
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> SetSignature
    for FiniteSubsetStructure<Set, SetB>
where
    Set::Elem: MetaType + Eq + Hash,
{
    type Elem = Set::Elem;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if !self.elems.contains(x) {
            return Err("element not in finite subset".to_string());
        }
        Ok(())
    }
}

impl<Set: EqSignature, SetB: BorrowedStructure<Set>> EqSignature
    for FiniteSubsetStructure<Set, SetB>
where
    Set::Elem: MetaType + Eq + Hash,
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().equal(a, b)
    }
}

impl<Set: PartialOrdSignature, SetB: BorrowedStructure<Set>> PartialOrdSignature
    for FiniteSubsetStructure<Set, SetB>
where
    Set::Elem: MetaType + Eq + Hash,
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().partial_cmp(a, b)
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> OrdSignature
    for FiniteSubsetStructure<Set, SetB>
where
    Set::Elem: MetaType + Eq + Hash,
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().cmp(a, b)
    }
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> CountableSetSignature
    for FiniteSubsetStructure<Set, SetB>
where
    Set::Elem: MetaType + Eq + Hash,
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> + Clone {
        self.elems.iter().cloned()
    }
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> FiniteSetSignature
    for FiniteSubsetStructure<Set, SetB>
where
    Set::Elem: MetaType + Eq + Hash,
{
}
