use algebraeon_structures::*;
use std::cmp::Ordering;
use std::collections::HashSet;
use std::hash::Hash;
use std::marker::PhantomData;

// A finite subset of a set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSubsetByHashStructure<Set: SetSignature, SetB: BorrowedStructure<Set>>
where
    Set::Elem: MetaType + Eq + Hash,
{
    _set: PhantomData<Set>,
    set: SetB,
    elems: HashSet<Set::Elem>,
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> FiniteSubsetByHashStructure<Set, SetB>
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

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> Signature
    for FiniteSubsetByHashStructure<Set, SetB>
where
    Set::Elem: MetaType + Eq + Hash,
{
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> SetSignature
    for FiniteSubsetByHashStructure<Set, SetB>
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
    for FiniteSubsetByHashStructure<Set, SetB>
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
    for FiniteSubsetByHashStructure<Set, SetB>
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
    for FiniteSubsetByHashStructure<Set, SetB>
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
    for FiniteSubsetByHashStructure<Set, SetB>
where
    Set::Elem: MetaType + Eq + Hash,
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.elems.into_iter()
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<Set: SetSignature, SetB: BorrowedStructure<Set>> FiniteSetSignature
    for FiniteSubsetByHashStructure<Set, SetB>
where
    Set::Elem: MetaType + Eq + Hash,
{
}
