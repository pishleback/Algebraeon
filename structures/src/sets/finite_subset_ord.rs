use crate::*;
use std::cmp::Ordering;
use std::rc::Rc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSubsetByOrd<Elem> {
    // ordered
    pub elems: Vec<Elem>,
}

impl<Elem> FiniteSubsetByOrd<Elem> {
    pub fn size(&self) -> usize {
        self.elems.len()
    }
}

// A finite subset of a set
#[derive(Debug, Clone)]
pub struct FiniteSubsetByOrdStructure<Set: OrdSignature> {
    set: Rc<Set>,
    // ordered
    subset: FiniteSubsetByOrd<Set::Elem>,
}

pub trait SetToFiniteSubsetByOrdSignature: OrdSignature {
    fn finite_subset(
        self: Rc<Self>,
        elems: Vec<Self::Elem>,
    ) -> Rc<FiniteSubsetByOrdStructure<Self>> {
        FiniteSubsetByOrdStructure::new(self, elems).into()
    }
}
impl<Set: OrdSignature> SetToFiniteSubsetByOrdSignature for Set {}

impl<Set: OrdSignature> PartialEq for FiniteSubsetByOrdStructure<Set> {
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

impl<Set: OrdSignature> Eq for FiniteSubsetByOrdStructure<Set> {}

impl<Set: OrdSignature> FiniteSubsetByOrdStructure<Set> {
    pub fn new(set: Rc<Set>, elems: Vec<Set::Elem>) -> Self {
        debug_assert!(set.is_sorted_and_unique(&elems));
        Self {
            set,
            subset: FiniteSubsetByOrd { elems },
        }
    }

    pub fn set(&self) -> &Rc<Set> {
        &self.set
    }
}

impl<Set: OrdSignature> Signature for FiniteSubsetByOrdStructure<Set> {}

impl<Set: OrdSignature> SetSignature for FiniteSubsetByOrdStructure<Set> {
    type Elem = Set::Elem;

    fn validate_element(self: &Rc<Self>, x: &Self::Elem) -> Result<(), String> {
        if !self.set().binary_search(&self.subset.elems, x) {
            return Err("element not in finite subset".to_string());
        }
        Ok(())
    }
}

impl<Set: OrdSignature> EqSignature for FiniteSubsetByOrdStructure<Set> {
    fn equal(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().equal(a, b)
    }
}

impl<Set: OrdSignature> PartialOrdSignature for FiniteSubsetByOrdStructure<Set> {
    fn partial_cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().partial_cmp(a, b)
    }
}

impl<Set: OrdSignature> OrdSignature for FiniteSubsetByOrdStructure<Set> {
    fn cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().cmp(a, b)
    }
}

impl<Set: OrdSignature> CountableSetSignature for FiniteSubsetByOrdStructure<Set> {
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        self.subset.elems.clone().into_iter()
    }
}

impl<Set: OrdSignature> FiniteSetSignature for FiniteSubsetByOrdStructure<Set> {
    fn size(self: &Rc<Self>) -> Natural {
        Natural::from(self.subset.size())
    }
}

impl<Set: OrderedFiniteSetSignature> OrderedFiniteSetSignature for FiniteSubsetByOrdStructure<Set> {
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unsized() {
        let set = i32::structure().finite_subset(vec![1, 2, 3, 4, 5, 6]);
        assert_eq!(set.size(), Natural::from(6usize));
    }

    #[test]
    fn test_sized() {
        let set = i32::structure()
            .finite_subset(vec![1, 2, 3, 4, 5, 6])
            .try_to_const_sized::<6>()
            .unwrap();
        assert_eq!(set.size(), Natural::from(6usize));
    }

    #[test]
    fn enumeration() {
        assert_enumerated_ord_finite_set!(
            i32::structure()
                .finite_subset(vec![1, 2, 3, 4, 5])
                .try_to_const_sized::<5>()
                .unwrap(),
            5
        );
    }
}
