use crate::*;
use std::cmp::Ordering;
use std::marker::PhantomData;

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
pub struct ConstSizeFiniteSubsetByOrdStructure<
    const N: usize,
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
> {
    _set: PhantomData<Set>,
    set: SetB,
    // ordered
    subset: ConstSizeFiniteSubsetByOrd<N, Set::Elem>,
}

pub trait SetToConstSizeFiniteSubsetByOrdSignature: OrdSignature {
    fn const_size_finite_subset<const N: usize>(
        &self,
        elems: [Self::Elem; N],
    ) -> ConstSizeFiniteSubsetByOrdStructure<N, Self, &Self> {
        ConstSizeFiniteSubsetByOrdStructure::new(self, elems)
    }

    fn into_const_size_finite_subset<const N: usize>(
        self,
        elems: [Self::Elem; N],
    ) -> ConstSizeFiniteSubsetByOrdStructure<N, Self, Self> {
        ConstSizeFiniteSubsetByOrdStructure::new(self, elems)
    }
}
impl<Set: OrdSignature> SetToConstSizeFiniteSubsetByOrdSignature for Set {}

impl<const N: usize, Set: OrdSignature, SetB: BorrowedStructure<Set>> PartialEq
    for ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
{
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

impl<const N: usize, Set: OrdSignature, SetB: BorrowedStructure<Set>> Eq
    for ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
{
}

impl<const N: usize, Set: OrdSignature, SetB: BorrowedStructure<Set>>
    ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
{
    pub fn new(set: SetB, elems: [Set::Elem; N]) -> Self {
        debug_assert!(set.borrow().is_sorted_and_unique(&elems));
        Self {
            _set: PhantomData,
            set,
            subset: ConstSizeFiniteSubsetByOrd { elems },
        }
    }

    pub fn set(&self) -> &Set {
        self.set.borrow()
    }
}

impl<const N: usize, Set: OrdSignature, SetB: BorrowedStructure<Set>> Signature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
{
}

impl<const N: usize, Set: OrdSignature, SetB: BorrowedStructure<Set>> SetSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
{
    type Elem = Set::Elem;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if !self.set().binary_search(&self.subset.elems, x) {
            return Err("element not in finite subset".to_string());
        }
        Ok(())
    }
}

impl<const N: usize, Set: OrdSignature, SetB: BorrowedStructure<Set>> EqSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().equal(a, b)
    }
}

impl<const N: usize, Set: OrdSignature, SetB: BorrowedStructure<Set>> PartialOrdSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().partial_cmp(a, b)
    }
}

impl<const N: usize, Set: OrdSignature, SetB: BorrowedStructure<Set>> OrdSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        self.set().cmp(a, b)
    }
}

impl<const N: usize, Set: OrdSignature, SetB: BorrowedStructure<Set>> CountableSetSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.subset.elems.into_iter()
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<const N: usize, Set: OrdSignature, SetB: BorrowedStructure<Set>> FiniteSetSignature
    for ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
{
    fn size(&self) -> Natural {
        Natural::from(N)
    }
}

impl<const N: usize, Set: OrderedFiniteSetSignature, SetB: BorrowedStructure<Set>>
    OrderedFiniteSetSignature for ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
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

impl<const N: usize, Set: OrderedFiniteSetSignature, SetB: BorrowedStructure<Set>>
    ConstSizeFiniteSetSignature<N> for ConstSizeFiniteSubsetByOrdStructure<N, Set, SetB>
{
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sized() {
        let set = i32::structure().into_const_size_finite_subset([1, 2, 3, 4, 5, 6]);
        assert_eq!(set.size(), Natural::from(6usize));
    }

    #[test]
    fn enumeration() {
        assert_enumerated_ord_finite_set!(
            i32::structure().into_const_size_finite_subset([1, 2, 3, 4, 5]),
            5
        );
    }
}
