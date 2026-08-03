use algebraeon_structures::*;
use std::cmp::Ordering;
use std::collections::HashSet;
use std::hash::Hash;
use std::rc::Rc;

// A finite subset of a set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSubsetByHashStructure<Set: SetSignature>
where
    Set::Elem: MetaType + Eq + Hash,
{
    set: Rc<Set>,
    elems: HashSet<Set::Elem>,
}

impl<Set: SetSignature> FiniteSubsetByHashStructure<Set>
where
    Set::Elem: MetaType + Eq + Hash,
{
    pub fn new(set: Rc<Set>, elems: HashSet<Set::Elem>) -> Rc<Self> {
        Self { set, elems }.into()
    }

    pub fn set(&self) -> &Rc<Set> {
        &self.set
    }
}

impl<Set: SetSignature> Signature for FiniteSubsetByHashStructure<Set> where
    Set::Elem: MetaType + Eq + Hash
{
}

impl<Set: SetSignature> SetSignature for FiniteSubsetByHashStructure<Set>
where
    Set::Elem: MetaType + Eq + Hash,
{
    type Elem = Set::Elem;

    fn validate_element(self: &Rc<Self>, x: &Self::Elem) -> Result<(), String> {
        if !self.elems.contains(x) {
            return Err("element not in finite subset".to_string());
        }
        Ok(())
    }
}

impl<Set: EqSignature> EqSignature for FiniteSubsetByHashStructure<Set>
where
    Set::Elem: MetaType + Eq + Hash,
{
    fn equal(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().equal(a, b)
    }
}

impl<Set: PartialOrdSignature> PartialOrdSignature for FiniteSubsetByHashStructure<Set>
where
    Set::Elem: MetaType + Eq + Hash,
{
    fn partial_cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().partial_cmp(a, b)
    }
}

impl<Set: OrdSignature> OrdSignature for FiniteSubsetByHashStructure<Set>
where
    Set::Elem: MetaType + Eq + Hash,
{
    fn cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().cmp(a, b)
    }
}

impl<Set: SetSignature> CountableSetSignature for FiniteSubsetByHashStructure<Set>
where
    Set::Elem: MetaType + Eq + Hash,
{
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        self.elems.clone().into_iter()
    }
}

impl<Set: SetSignature> FiniteSetSignature for FiniteSubsetByHashStructure<Set> where
    Set::Elem: MetaType + Eq + Hash
{
}
