use crate::{
    combinatorics::subsets_colex,
    sets::{FiniteSubsetByOrd, SetToFiniteSubsetsByOrdSignature},
};
use algebraeon_structures::*;
use std::{cmp::Ordering, marker::PhantomData};

/// The set of all k-element subsets of a set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FixedSizeFiniteSubsetsByOrdStructure<Set: OrdSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
    k: usize,
}

pub trait SetToFixedSizeFiniteSubsetsByOrdSignature: OrdSignature {
    fn fixed_size_finite_subsets(
        &self,
        k: usize,
    ) -> FixedSizeFiniteSubsetsByOrdStructure<Self, &Self> {
        FixedSizeFiniteSubsetsByOrdStructure::new(self, k)
    }

    fn into_fixed_size_finite_subsets(
        self,
        k: usize,
    ) -> FixedSizeFiniteSubsetsByOrdStructure<Self, Self> {
        FixedSizeFiniteSubsetsByOrdStructure::new(self, k)
    }
}
impl<Set: OrdSignature> SetToFixedSizeFiniteSubsetsByOrdSignature for Set {}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>>
    FixedSizeFiniteSubsetsByOrdStructure<Set, SetB>
{
    pub fn new(set: SetB, k: usize) -> Self {
        Self {
            _set: PhantomData,
            set,
            k,
        }
    }

    pub fn set(&self) -> &Set {
        self.set.borrow()
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> Signature
    for FixedSizeFiniteSubsetsByOrdStructure<Set, SetB>
{
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> SetSignature
    for FixedSizeFiniteSubsetsByOrdStructure<Set, SetB>
{
    type Elem = FiniteSubsetByOrd<Set>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if !self.set().is_sorted_and_unique(&x.elems) {
            return Err("elems is not sorted and unique".to_string());
        }
        if x.size() != self.k {
            return Err("subset has the wrong size".to_string());
        }
        Ok(())
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> EqSignature
    for FixedSizeFiniteSubsetsByOrdStructure<Set, SetB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.set().finite_subsets().equal(a, b)
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> PartialOrdSignature
    for FixedSizeFiniteSubsetsByOrdStructure<Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        self.set().finite_subsets().partial_cmp(a, b)
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> OrdSignature
    for FixedSizeFiniteSubsetsByOrdStructure<Set, SetB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        self.set().finite_subsets().cmp(a, b)
    }
}

impl<Set: OrdSignature + CountableSetSignature, SetB: BorrowedStructure<Set>> CountableSetSignature
    for FixedSizeFiniteSubsetsByOrdStructure<Set, SetB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        // if the set has more than 64 elements then we'll never generate subsets including anything beyond the 64th element, so this is fine
        let elems = self
            .set()
            .generate_all_elements()
            .take(64)
            .collect::<Vec<_>>();
        subsets_colex(elems.len(), self.k).map(move |idx_subset| FiniteSubsetByOrd {
            elems: idx_subset
                .into_iter()
                .map(|idx| elems[idx].clone())
                .collect(),
        })
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<Set: OrdSignature + FiniteSetSignature, SetB: BorrowedStructure<Set>> FiniteSetSignature
    for FixedSizeFiniteSubsetsByOrdStructure<Set, SetB>
{
    fn size(&self) -> Natural {
        choose(self.set().size(), Natural::from(self.k))
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>
    EnumeratedOrdFiniteSetSignature for FixedSizeFiniteSubsetsByOrdStructure<Set, SetB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.list_all_elements()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        // colex rank
        let mut t = Natural::ZERO;
        for i in 0..self.k {
            t += choose(
                self.set().element_to_enumeration(&elem.elems[i]),
                Natural::from(i + 1),
            );
        }
        t
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        let n = self.set().size();
        if *num >= choose(&n, Natural::from(self.k)) {
            return None;
        }

        let mut elems = self
            .set()
            .list_all_elements_ordered()
            .into_iter()
            .map(Some)
            .collect::<Vec<_>>();
        debug_assert_eq!(Natural::from(elems.len()), n);

        let mut r = num.clone();
        let mut result = vec![None; self.k];

        for i in (0..self.k).rev() {
            let mut x = i;
            while Natural::from(x + 1) < n
                && choose(Natural::from(x + 1), Natural::from(i + 1)) <= r
            {
                x += 1;
            }
            r -= choose(Natural::from(x), Natural::from(i + 1));
            result[i] = Some(elems[x].take().unwrap());
        }

        Some(FiniteSubsetByOrd {
            elems: result.into_iter().map(|x| x.unwrap()).collect(),
        })
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>>
    FixedSizeFiniteSubsetsByOrdStructure<Set, SetB>
{
    pub fn subset(&self, elems: Vec<Set::Elem>) -> <Self as SetSignature>::Elem {
        for elem in &elems {
            #[cfg(debug_assertions)]
            self.set().validate_element(elem).unwrap();
        }
        FiniteSubsetByOrd {
            elems: self.set().unique(self.set().sort(elems)),
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::sets::{SetToFiniteSubsetByOrdSignature, SetToFixedSizeFiniteSubsetsByOrdSignature};
    use algebraeon_structures::*;

    #[test]
    fn test_enumerate_0_0() {
        let set = i32::structure().into_finite_subset(vec![]);
        let set_subsets = set.fixed_size_finite_subsets(0);

        let subsets = set_subsets.generate_all_elements().collect::<Vec<_>>();
        assert_eq!(subsets.len(), 1);
        assert_eq!(set_subsets.size(), Natural::from(1usize));

        assert_eq!(
            set_subsets.element_to_enumeration(&subsets[0]),
            Natural::from(0usize)
        );
        assert!(
            set_subsets.equal(
                &set_subsets
                    .enumeration_to_element(&Natural::from(0usize))
                    .unwrap(),
                &subsets[0]
            )
        );
        assert!(
            set_subsets
                .enumeration_to_element(&Natural::from(1usize))
                .is_none()
        );
    }

    #[test]
    fn test_enumerate_0_1() {
        let set = i32::structure().into_finite_subset(vec![]);
        let set_subsets = set.fixed_size_finite_subsets(1);

        let subsets = set_subsets.generate_all_elements().collect::<Vec<_>>();
        assert_eq!(subsets.len(), 0);
        assert_eq!(set_subsets.size(), Natural::from(0usize));

        assert!(
            set_subsets
                .enumeration_to_element(&Natural::from(0usize))
                .is_none()
        );
    }

    #[test]
    fn test_enumerate_5_3() {
        let set = i32::structure().into_finite_subset(vec![1, 2, 3, 4, 5]);
        let set_subsets = set.fixed_size_finite_subsets(3);

        let subsets = set_subsets.generate_all_elements().collect::<Vec<_>>();
        assert_eq!(subsets.len(), 10);
        assert_eq!(Natural::from(subsets.len()), set_subsets.size());
        for i in 0..9 {
            assert!(set_subsets.cmp(&subsets[i], &subsets[i + 1]).is_lt());
        }

        for (idx, subset) in set_subsets.generate_all_elements().enumerate() {
            println!(
                "{:?} {:?} {:?} {:?}",
                idx,
                subset,
                set_subsets.element_to_enumeration(&subset),
                set_subsets.enumeration_to_element(&Natural::from(idx)),
            );
            assert_eq!(
                Natural::from(idx),
                set_subsets.element_to_enumeration(&subset)
            );
            assert_eq!(
                set_subsets.enumeration_to_element(&Natural::from(idx)),
                Some(subset)
            );
        }
        assert_eq!(
            set_subsets.enumeration_to_element(&Natural::from(10usize)),
            None
        );
    }

    // #[test]
    // fn test_cmp() {
    //     let set = i32::structure().into_finite_subset(vec![1, 2, 3, 4]);
    //     let subsets = set.finite_subsets();

    //     // eq
    //     assert!(subsets.equal(
    //         &subsets.subset(vec![3, 3, 3, 2, 2, 1]),
    //         &subsets.subset(vec![1, 2, 3])
    //     ));
    //     assert!(subsets.equal(&subsets.subset(vec![4]), &subsets.subset(vec![4])));
    //     assert!(subsets.equal(&subsets.subset(vec![]), &subsets.subset(vec![])));

    //     // lt
    //     assert!(
    //         subsets
    //             .cmp(&subsets.subset(vec![]), &subsets.subset(vec![1]))
    //             .is_lt()
    //     );
    //     assert!(
    //         subsets
    //             .cmp(&subsets.subset(vec![]), &subsets.subset(vec![1, 2, 3, 4]))
    //             .is_lt()
    //     );
    //     assert!(
    //         subsets
    //             .cmp(&subsets.subset(vec![2, 3]), &subsets.subset(vec![3, 4]))
    //             .is_lt()
    //     );
    //     assert!(
    //         subsets
    //             .cmp(&subsets.subset(vec![2, 3]), &subsets.subset(vec![1, 2]))
    //             .is_gt()
    //     );
    // }
}
