use crate::*;
use std::cmp::Ordering;
use std::rc::Rc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeFiniteSubsetByOrd<const N: usize, Elem> {
    // ordered
    pub elems: [Elem; N],
}

impl<const N: usize, Elem> From<ConstSizeFiniteSubsetByOrd<N, Elem>> for FiniteSubsetByOrd<Elem> {
    fn from(value: ConstSizeFiniteSubsetByOrd<N, Elem>) -> Self {
        Self {
            elems: value.elems.into(),
        }
    }
}

impl<const N: usize, Set: OrdSignature> ConstSizeFiniteSubsetByOrd<N, Set> {
    pub fn size(&self) -> usize {
        N
    }
}

// A finite subset of a set
#[derive(Debug, Clone)]
pub struct ConstSizeFiniteSubsetByOrdStructure<const N: usize, Set: OrdSignature> {
    set: Rc<Set>,
    // ordered
    subset: ConstSizeFiniteSubsetByOrd<N, Set::Elem>,
}

pub trait SetToConstSizeFiniteSubsetByOrdSignature: OrdSignature {
    fn const_size_finite_subset<const N: usize>(
        self: Rc<Self>,
        elems: [Self::Elem; N],
    ) -> Rc<ConstSizeFiniteSubsetByOrdStructure<N, Self>> {
        ConstSizeFiniteSubsetByOrdStructure::new(self, elems).into()
    }
}
impl<Set: OrdSignature> SetToConstSizeFiniteSubsetByOrdSignature for Set {}

impl<const N: usize, Set: OrdSignature> PartialEq for ConstSizeFiniteSubsetByOrdStructure<N, Set> {
    fn eq(&self, other: &Self) -> bool {
        let set = self.set();
        if set != other.set() {
            return false;
        }
        for item in self
            .set()
            .merge_sorted_and_unique(
                self.subset.elems.iter().collect(),
                other.subset.elems.iter().collect(),
            )
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
        {
            match item {
                MergedUniqueSource::First(_) | MergedUniqueSource::Second(_) => {
                    return false;
                }
                MergedUniqueSource::Both(_, _) => {}
            }
        }
        true
    }
}

impl<const N: usize, Set: OrdSignature> Eq for ConstSizeFiniteSubsetByOrdStructure<N, Set> {}

impl<const N: usize, Set: OrdSignature> ConstSizeFiniteSubsetByOrdStructure<N, Set> {
    pub fn new(set: Rc<Set>, elems: [Set::Elem; N]) -> Self {
        debug_assert!(set.is_sorted_and_unique(&elems));
        Self {
            set,
            subset: ConstSizeFiniteSubsetByOrd { elems },
        }
    }

    pub fn set(&self) -> &Rc<Set> {
        &self.set
    }
}

impl<const N: usize, Set: OrdSignature> Signature for ConstSizeFiniteSubsetByOrdStructure<N, Set> {}

impl<const N: usize, Set: OrdSignature> SetSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set>
{
    type Elem = Set::Elem;

    fn validate_element(self: &Rc<Self>, x: &Self::Elem) -> Result<(), String> {
        if !self.set().binary_search(&self.subset.elems, x) {
            return Err("element not in finite subset".to_string());
        }
        Ok(())
    }
}

impl<const N: usize, Set: OrdSignature> EqSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set>
{
    fn equal(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().equal(a, b)
    }
}

impl<const N: usize, Set: OrdSignature> PartialOrdSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set>
{
    fn partial_cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().partial_cmp(a, b)
    }
}

impl<const N: usize, Set: OrdSignature> OrdSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set>
{
    fn cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().cmp(a, b)
    }
}

impl<const N: usize, Set: OrdSignature> CountableSetSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set>
{
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        self.subset.elems.clone().into_iter()
    }
}

impl<const N: usize, Set: OrdSignature> FiniteSetSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set>
{
    fn size(self: &Rc<Self>) -> Natural {
        Natural::from(N)
    }
}

impl<const N: usize, Set: OrderedFiniteSetSignature> OrderedFiniteSetSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set>
{
    fn list_all_elements_ordered(self: &Rc<Self>) -> Vec<Self::Elem> {
        self.list_all_elements()
    }

    fn element_to_enumeration(self: &Rc<Self>, elem: &Self::Elem) -> Natural {
        #[cfg(debug_assertions)]
        self.validate_element(elem).unwrap();
        Natural::from(
            self.set()
                .binary_search_index(&self.subset.elems, elem)
                .unwrap(),
        )
    }

    fn enumeration_to_element(self: &Rc<Self>, num: &Natural) -> Option<Self::Elem> {
        let num: Result<usize, ()> = num.try_into();
        if let Ok(num) = num {
            self.subset.elems.get(num).cloned()
        } else {
            None
        }
    }
}

impl<const N: usize, Set: OrderedFiniteSetSignature> ConstSizeFiniteSetSignature<N>
    for ConstSizeFiniteSubsetByOrdStructure<N, Set>
{
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sized() {
        let set = i32::structure().const_size_finite_subset([1, 2, 3, 4, 5, 6]);
        assert_eq!(set.size(), Natural::from(6usize));
    }

    #[test]
    fn enumeration() {
        assert_enumerated_ord_finite_set!(
            i32::structure().const_size_finite_subset([1, 2, 3, 4, 5]),
            5
        );
    }
}
