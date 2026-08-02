use crate::{combinatorics::subsets_colex, sets::SetToFiniteSubsetsByOrdSignature};
use algebraeon_structures::*;
use std::{cmp::Ordering, sync::Arc};

/// The set of all k-element subsets of a set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FixedSizeFiniteSubsetsByOrdStructure<Set: OrdSignature> {
    set: Arc<Set>,
    k: usize,
}

pub trait SetToFixedSizeFiniteSubsetsByOrdSignature: OrdSignature {
    fn fixed_size_finite_subsets(
        self: &Arc<Self>,
        k: usize,
    ) -> Arc<FixedSizeFiniteSubsetsByOrdStructure<Self>> {
        FixedSizeFiniteSubsetsByOrdStructure::new(self.clone(), k)
    }
}
impl<Set: OrdSignature> SetToFixedSizeFiniteSubsetsByOrdSignature for Set {}

impl<Set: OrdSignature> FixedSizeFiniteSubsetsByOrdStructure<Set> {
    pub fn new(set: Arc<Set>, k: usize) -> Arc<Self> {
        Self { set, k }.into()
    }

    pub fn set(&self) -> &Arc<Set> {
        &self.set
    }
}

impl<Set: OrdSignature> Signature for FixedSizeFiniteSubsetsByOrdStructure<Set> {}

impl<Set: OrdSignature> SetSignature for FixedSizeFiniteSubsetsByOrdStructure<Set> {
    type Elem = FiniteSubsetByOrd<Set::Elem>;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        if !self.set().is_sorted_and_unique(&x.elems) {
            return Err("elems is not sorted and unique".to_string());
        }
        if x.size() != self.k {
            return Err("subset has the wrong size".to_string());
        }
        Ok(())
    }
}

impl<Set: OrdSignature> EqSignature for FixedSizeFiniteSubsetsByOrdStructure<Set> {
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.set().finite_subsets().equal(a, b)
    }
}

impl<Set: OrdSignature> PartialOrdSignature for FixedSizeFiniteSubsetsByOrdStructure<Set> {
    fn partial_cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        self.set().finite_subsets().partial_cmp(a, b)
    }
}

impl<Set: OrdSignature> OrdSignature for FixedSizeFiniteSubsetsByOrdStructure<Set> {
    fn cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        self.set().finite_subsets().cmp(a, b)
    }
}

impl<Set: OrdSignature + CountableSetSignature> CountableSetSignature
    for FixedSizeFiniteSubsetsByOrdStructure<Set>
{
    fn generate_all_elements(self: Arc<Self>) -> impl Iterator<Item = Self::Elem> {
        // if the set has more than 64 elements then we'll never generate subsets including anything beyond the 64th element, so this is fine
        let elems = self
            .set()
            .clone()
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
}

impl<Set: OrdSignature + FiniteSetSignature> FiniteSetSignature
    for FixedSizeFiniteSubsetsByOrdStructure<Set>
{
    fn size(self: &Arc<Self>) -> Natural {
        choose(self.set().size(), Natural::from(self.k))
    }
}

impl<Set: OrderedFiniteSetSignature> OrderedFiniteSetSignature
    for FixedSizeFiniteSubsetsByOrdStructure<Set>
{
    fn list_all_elements_ordered(self: &Arc<Self>) -> Vec<Self::Elem> {
        self.list_all_elements()
    }

    fn element_to_enumeration(self: &Arc<Self>, elem: &Self::Elem) -> Natural {
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

    fn enumeration_to_element(self: &Arc<Self>, num: &Natural) -> Option<Self::Elem> {
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

impl<Set: OrdSignature> FixedSizeFiniteSubsetsByOrdStructure<Set> {
    pub fn subset(&self, elems: Vec<Set::Elem>) -> <Self as SetSignature>::Elem {
        #[cfg(debug_assertions)]
        for elem in &elems {
            self.set().validate_element(elem).unwrap();
        }
        FiniteSubsetByOrd {
            elems: self.set().unique(self.set().sort(elems)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enumerate() {
        algebraeon_structures::assert_enumerated_ord_finite_set!(
            i32::structure()
                .finite_subset(vec![])
                .fixed_size_finite_subsets(0),
            1
        );

        algebraeon_structures::assert_enumerated_ord_finite_set!(
            i32::structure()
                .finite_subset(vec![])
                .fixed_size_finite_subsets(1),
            0
        );

        algebraeon_structures::assert_enumerated_ord_finite_set!(
            i32::structure()
                .finite_subset(vec![1, 2, 3, 4, 5])
                .fixed_size_finite_subsets(3),
            10
        );
    }
}
