use algebraeon_structures::*;
use std::{cmp::Ordering, sync::Arc};

use crate::combinatorics::all_subsets;

/// The set of all finite subsets of a set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSubsetsByOrdStructure<Set: OrdSignature> {
    set: Arc<Set>,
}

pub trait SetToFiniteSubsetsByOrdSignature: OrdSignature {
    fn finite_subsets(self: &Arc<Self>) -> Arc<FiniteSubsetsByOrdStructure<Self>> {
        FiniteSubsetsByOrdStructure::new(self.clone())
    }
}
impl<Set: OrdSignature> SetToFiniteSubsetsByOrdSignature for Set {}

impl<Set: OrdSignature> FiniteSubsetsByOrdStructure<Set> {
    pub fn new(set: Arc<Set>) -> Arc<Self> {
        Self { set }.into()
    }

    pub fn set(&self) -> &Arc<Set> {
        &self.set
    }
}

impl<Set: OrdSignature> Signature for FiniteSubsetsByOrdStructure<Set> {}

impl<Set: OrdSignature> SetSignature for FiniteSubsetsByOrdStructure<Set> {
    type Elem = FiniteSubsetByOrd<Set::Elem>;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        if !self.set().is_sorted_and_unique(&x.elems) {
            return Err("elems is not sorted and unique".to_string());
        }
        Ok(())
    }
}

impl<Set: OrdSignature> EqSignature for FiniteSubsetsByOrdStructure<Set> {
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.cmp(a, b).is_eq()
    }
}

impl<Set: OrdSignature> PartialOrdSignature for FiniteSubsetsByOrdStructure<Set> {
    fn partial_cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<Set: OrdSignature> OrdSignature for FiniteSubsetsByOrdStructure<Set> {
    fn cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        // needs to be such that a < b iff self.element_to_enumeration(a) < self.element_to_enumeration(b)
        // element_to_enumeration assigns a binary number with bits set based on the enumeration on the underlying set
        // so here we need to implement a generalized binary number comparission
        for item in self
            .set()
            .merge_sorted_and_unique(a.elems.iter().collect(), b.elems.iter().collect())
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
        {
            match item {
                MergedUniqueSource::First(_) => {
                    return Ordering::Greater;
                }
                MergedUniqueSource::Second(_) => {
                    return Ordering::Less;
                }
                MergedUniqueSource::Both(_, _) => {}
            }
        }
        Ordering::Equal
    }
}

impl<Set: OrdSignature + CountableSetSignature> CountableSetSignature
    for FiniteSubsetsByOrdStructure<Set>
{
    fn generate_all_elements(self: Arc<Self>) -> impl Iterator<Item = Self::Elem> {
        // if the set has more than 64 elements then we'll never generate subsets including anything beyond the 64th element, so this is fine
        let elems = self
            .set()
            .clone()
            .generate_all_elements()
            .take(64)
            .collect::<Vec<_>>();
        all_subsets(elems.len()).map(move |idx_subset| FiniteSubsetByOrd {
            elems: idx_subset
                .into_iter()
                .map(|idx| elems[idx].clone())
                .collect(),
        })
    }
}

impl<Set: OrdSignature + FiniteSetSignature> FiniteSetSignature
    for FiniteSubsetsByOrdStructure<Set>
{
    fn size(self: &Arc<Self>) -> Natural {
        Natural::TWO.pow(&self.set().size())
    }
}

impl<Set: OrderedFiniteSetSignature> OrderedFiniteSetSignature
    for FiniteSubsetsByOrdStructure<Set>
{
    fn list_all_elements_ordered(self: &Arc<Self>) -> Vec<Self::Elem> {
        self.list_all_elements()
    }

    fn element_to_enumeration(self: &Arc<Self>, elem: &Self::Elem) -> Natural {
        let mut t = Natural::ZERO;
        for x in &elem.elems {
            let i: usize = self
                .set()
                .element_to_enumeration(x)
                .try_into()
                .expect("too large");
            t |= Natural::ONE << i;
        }
        t
    }

    fn enumeration_to_element(self: &Arc<Self>, num: &Natural) -> Option<Self::Elem> {
        let len = num.bitcount();
        if Natural::from(len) > self.set().size() {
            return None;
        }
        Some(FiniteSubsetByOrd {
            elems: (0..len)
                .filter_map(move |i| {
                    if (num >> i) & Natural::ONE == Natural::ZERO {
                        None
                    } else {
                        Some(
                            self.set()
                                .enumeration_to_element(&Natural::from(i))
                                .unwrap(),
                        )
                    }
                })
                .collect(),
        })
    }
}

impl<Set: OrdSignature> FiniteSubsetsByOrdStructure<Set> {
    pub fn subset(self: &Arc<Self>, elems: Vec<Set::Elem>) -> <Self as SetSignature>::Elem {
        #[cfg(debug_assertions)]
        for elem in &elems {
            self.set().validate_element(elem).unwrap();
        }
        FiniteSubsetByOrd {
            elems: self.set().unique(self.set().sort(elems)),
        }
    }

    pub fn is_disjoint(
        self: &Arc<Self>,
        subset_1: &<Self as SetSignature>::Elem,
        subset_2: &<Self as SetSignature>::Elem,
    ) -> bool {
        debug_assert!(self.validate_element(subset_1).is_ok());
        debug_assert!(self.validate_element(subset_2).is_ok());
        for item in self
            .set()
            .merge_sorted_and_unique(subset_1.elems.clone(), subset_2.elems.clone())
        {
            match item {
                MergedUniqueSource::First(_) | MergedUniqueSource::Second(_) => {}
                MergedUniqueSource::Both(_, _) => {
                    return false;
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use crate::sets::SetToFiniteSubsetsByOrdSignature;
    use algebraeon_structures::*;

    #[test]
    fn test_enumerate() {
        algebraeon_structures::assert_enumerated_ord_finite_set!(
            i32::structure().finite_subset(vec![]).finite_subsets(),
            1
        );

        algebraeon_structures::assert_enumerated_ord_finite_set!(
            i32::structure().finite_subset(vec![1]).finite_subsets(),
            2
        );

        algebraeon_structures::assert_enumerated_ord_finite_set!(
            i32::structure()
                .finite_subset(vec![1, 2, 3, 4, 5, 6])
                .finite_subsets(),
            64
        );
    }

    #[test]
    fn test_cmp() {
        let set = i32::structure().finite_subset(vec![1, 2, 3, 4]);
        let subsets = set.finite_subsets();

        // eq
        assert!(subsets.equal(
            &subsets.subset(vec![3, 3, 3, 2, 2, 1]),
            &subsets.subset(vec![1, 2, 3])
        ));
        assert!(subsets.equal(&subsets.subset(vec![4]), &subsets.subset(vec![4])));
        assert!(subsets.equal(&subsets.subset(vec![]), &subsets.subset(vec![])));

        // lt
        assert!(
            subsets
                .cmp(&subsets.subset(vec![]), &subsets.subset(vec![1]))
                .is_lt()
        );
        assert!(
            subsets
                .cmp(&subsets.subset(vec![]), &subsets.subset(vec![1, 2, 3, 4]))
                .is_lt()
        );
        assert!(
            subsets
                .cmp(&subsets.subset(vec![2, 3]), &subsets.subset(vec![3, 4]))
                .is_lt()
        );
        assert!(
            subsets
                .cmp(&subsets.subset(vec![2, 3]), &subsets.subset(vec![1, 2]))
                .is_gt()
        );
    }
}
