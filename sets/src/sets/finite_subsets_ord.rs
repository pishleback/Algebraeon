use crate::combinatorics::all_subsets;
use algebraeon_structures::*;
use std::{cmp::Ordering, marker::PhantomData};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSubsetByOrd<Set: SetSignature> {
    // ordered
    pub elems: Vec<Set::Elem>,
}

// The set of finite subsets of a set
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSubsetsByOrdStructure<Set: OrdSignature, SetB: BorrowedStructure<Set>> {
    _set: PhantomData<Set>,
    set: SetB,
}

pub trait SetToFiniteSubsetsByOrdSignature: OrdSignature {
    fn finite_subsets(&self) -> FiniteSubsetsByOrdStructure<Self, &Self> {
        FiniteSubsetsByOrdStructure::new(self)
    }

    fn into_finite_subsets(self) -> FiniteSubsetsByOrdStructure<Self, Self> {
        FiniteSubsetsByOrdStructure::new(self)
    }
}
impl<Set: OrdSignature> SetToFiniteSubsetsByOrdSignature for Set {}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> FiniteSubsetsByOrdStructure<Set, SetB> {
    pub fn new(set: SetB) -> Self {
        Self {
            _set: PhantomData,
            set,
        }
    }

    pub fn set(&self) -> &Set {
        self.set.borrow()
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> Signature
    for FiniteSubsetsByOrdStructure<Set, SetB>
{
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> SetSignature
    for FiniteSubsetsByOrdStructure<Set, SetB>
{
    type Elem = FiniteSubsetByOrd<Set>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if !self.set().is_sorted_and_unique(&x.elems) {
            return Err("elems is not sorted and unique".to_string());
        }
        Ok(())
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> EqSignature
    for FiniteSubsetsByOrdStructure<Set, SetB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.cmp(a, b).is_eq()
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> PartialOrdSignature
    for FiniteSubsetsByOrdStructure<Set, SetB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> OrdSignature
    for FiniteSubsetsByOrdStructure<Set, SetB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
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

impl<Set: OrdSignature + CountableSetSignature, SetB: BorrowedStructure<Set>> CountableSetSignature
    for FiniteSubsetsByOrdStructure<Set, SetB>
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        // if the set has more than 64 elements then we'll never generate subsets including anything beyond the 64th element, so this is fine
        let elems = self
            .set()
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

impl<Set: OrdSignature + FiniteSetSignature, SetB: BorrowedStructure<Set>> FiniteSetSignature
    for FiniteSubsetsByOrdStructure<Set, SetB>
{
    fn size(&self) -> Natural {
        Natural::TWO.pow(&self.set().size())
    }
}

impl<Set: EnumeratedOrdFiniteSetSignature, SetB: BorrowedStructure<Set>>
    EnumeratedOrdFiniteSetSignature for FiniteSubsetsByOrdStructure<Set, SetB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.list_all_elements()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
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

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
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

impl<Set: OrdSignature, SetB: BorrowedStructure<Set>> FiniteSubsetsByOrdStructure<Set, SetB> {
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
    use crate::sets::{SetToFiniteSubsetByOrdSignature, SetToFiniteSubsetsByOrdSignature};
    use algebraeon_structures::*;

    #[test]
    fn test_enumerate_0() {
        let set = i32::structure().into_finite_subset(vec![]);
        let subsets = set.finite_subsets();
        for (idx, subset) in subsets.generate_all_elements().enumerate() {
            assert_eq!(Natural::from(idx), subsets.element_to_enumeration(&subset));
            assert_eq!(
                subsets.enumeration_to_element(&Natural::from(idx)),
                Some(subset)
            );
        }
        assert_eq!(subsets.enumeration_to_element(&Natural::from(1usize)), None);
    }

    #[test]
    fn test_enumerate_1() {
        let set = i32::structure().into_finite_subset(vec![1]);
        let subsets = set.finite_subsets();
        for (idx, subset) in subsets.generate_all_elements().enumerate() {
            assert_eq!(Natural::from(idx), subsets.element_to_enumeration(&subset));
            assert_eq!(
                subsets.enumeration_to_element(&Natural::from(idx)),
                Some(subset)
            );
        }
        assert_eq!(subsets.enumeration_to_element(&Natural::from(2usize)), None);
    }

    #[test]
    fn test_enumerate_6() {
        let set = i32::structure().into_finite_subset(vec![1, 2, 3, 4, 5, 6]);
        let subsets = set.finite_subsets();
        for (idx, subset) in subsets.generate_all_elements().enumerate() {
            assert_eq!(Natural::from(idx), subsets.element_to_enumeration(&subset));
            assert_eq!(
                subsets.enumeration_to_element(&Natural::from(idx)),
                Some(subset)
            );
        }
        assert_eq!(
            subsets.enumeration_to_element(&Natural::from(64usize)),
            None
        );
    }

    #[test]
    fn test_cmp() {
        let set = i32::structure().into_finite_subset(vec![1, 2, 3, 4]);
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
