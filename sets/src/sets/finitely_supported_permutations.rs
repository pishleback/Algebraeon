use algebraeon_structures::*;
use itertools::Itertools;
use std::borrow::Borrow;
use std::collections::HashSet;
use std::hash::Hash;
use std::sync::Arc;

#[derive(Debug, Clone)]
struct FromAndTo {
    from: usize, // index of the thing mapping to an item
    to: usize,   // index of the thing an item maps to
}

#[derive(Debug, Clone)]
pub struct FinitelySupportedPermutation<T> {
    perm: Vec<(T, FromAndTo)>,
}

impl<Elem: MetaType> PartialEq for FinitelySupportedPermutation<Elem>
where
    Elem::Signature: OrdSignature,
{
    fn eq(&self, other: &Self) -> bool {
        self.equal(other)
    }
}

impl<Elem: MetaType> Eq for FinitelySupportedPermutation<Elem> where Elem::Signature: OrdSignature {}

impl<Elem: MetaType + Hash> Hash for FinitelySupportedPermutation<Elem>
where
    Elem::Signature: OrdSignature,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        for (elem, mapping) in &self.perm {
            elem.hash(state);
            mapping.to.hash(state);
        }
    }
}

impl<T> FinitelySupportedPermutation<T> {
    fn identity() -> Self {
        Self { perm: vec![] }
    }

    fn inverse(self) -> Self {
        Self {
            perm: self
                .perm
                .into_iter()
                .map(|(elem, mapping)| {
                    (
                        elem,
                        FromAndTo {
                            from: mapping.to,
                            to: mapping.from,
                        },
                    )
                })
                .collect(),
        }
    }
}

impl<Elem: MetaType> MetaType for FinitelySupportedPermutation<Elem>
where
    Elem::Signature: OrdSignature,
{
    type Signature = FinitelySupportedPermutationsStructure<Elem::Signature>;

    fn structure() -> Arc<Self::Signature> {
        FinitelySupportedPermutationsStructure::new(Elem::structure())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelySupportedPermutationsStructure<Set: SetSignature> {
    set: Arc<Set>,
}

impl<Set: SetSignature> FinitelySupportedPermutationsStructure<Set> {
    pub fn new(set: Arc<Set>) -> Arc<Self> {
        Self { set }.into()
    }
}

pub trait SetToFinitelySupportedPermutationsStructure: SetSignature {
    fn finitely_supported_permutations(
        self: &Arc<Self>,
    ) -> Arc<FinitelySupportedPermutationsStructure<Self>> {
        FinitelySupportedPermutationsStructure::new(self.clone())
    }
}
impl<Set: SetSignature> SetToFinitelySupportedPermutationsStructure for Set {}

pub trait FiniteSetToFinitelySupportedPermutationsStructure: FiniteSetSignature {
    fn permutations(self: &Arc<Self>) -> Arc<FinitelySupportedPermutationsStructure<Self>> {
        FinitelySupportedPermutationsStructure::new(self.clone())
    }
}
impl<Set: FiniteSetSignature> FiniteSetToFinitelySupportedPermutationsStructure for Set {}

impl<Set: SetSignature> Signature for FinitelySupportedPermutationsStructure<Set> {}

impl<Set: OrdSignature> SetSignature for FinitelySupportedPermutationsStructure<Set> {
    type Elem = FinitelySupportedPermutation<Set::Elem>;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        if !self.set().is_sorted_and_unique_by_key(&x.perm, |(a, _)| a) {
            return Err("Permutation first items are not sorted and unique".to_string());
        }
        if (0..x.perm.len()).collect::<HashSet<_>>()
            != x.perm.iter().map(|item| item.1.to).collect::<HashSet<_>>()
        {
            return Err("Permutation images are not valid".to_string());
        }
        if (0..x.perm.len()).collect::<HashSet<_>>()
            != x.perm
                .iter()
                .map(|item| item.1.from)
                .collect::<HashSet<_>>()
        {
            return Err("Permutation preimages are not valid".to_string());
        }
        for item_1 in &x.perm {
            let item_2 = &x.perm[item_1.1.to];
            if self.set().equal(&item_1.0, &item_2.0) {
                return Err("Permutation includes a redundant fixed point".to_string());
            }
            let item_3 = &x.perm[item_2.1.from];
            if !self.set().equal(&item_1.0, &item_3.0) {
                return Err("Permutation image followed by preimage is not identity".to_string());
            }
        }
        for item_1 in &x.perm {
            let item_2 = &x.perm[item_1.1.from];
            if self.set().equal(&item_1.0, &item_2.0) {
                return Err("Permutation includes a redundant fixed point".to_string());
            }
            let item_3 = &x.perm[item_2.1.to];
            if !self.set().equal(&item_1.0, &item_3.0) {
                return Err("Permutation preimage followed by image is not identity".to_string());
            }
        }
        Ok(())
    }
}

impl<Set: OrdSignature> EqSignature for FinitelySupportedPermutationsStructure<Set> {
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        let mut a_tos = vec![];
        let mut b_tos = vec![];
        for merged in self.set().merge_sorted_and_unique_by_key(
            a.perm.iter().collect(),
            b.perm.iter().collect(),
            |item| &item.0,
        ) {
            match merged {
                MergedUniqueSource::First(_) | MergedUniqueSource::Second(_) => {
                    return false;
                }
                MergedUniqueSource::Both(a_item, b_item) => {
                    a_tos.push(a_item.1.to);
                    b_tos.push(b_item.1.to);
                }
            }
        }
        a_tos == b_tos
    }
}

impl<Set: OrdSignature> FinitelySupportedPermutationsStructure<Set> {
    pub fn map<T>(
        self: &Arc<Self>,
        perm: <Self as SetSignature>::Elem,
        f: impl Fn(Set::Elem) -> T,
    ) -> FinitelySupportedPermutation<T> {
        FinitelySupportedPermutation {
            perm: perm
                .perm
                .into_iter()
                .map(|(elem, from_and_to)| (f(elem), from_and_to))
                .collect(),
        }
    }

    fn image_ref<'a>(
        self: &Arc<Self>,
        perm: &'a <Self as SetSignature>::Elem,
        elem: &'a Set::Elem,
    ) -> &'a Set::Elem {
        debug_assert!(self.is_element(perm));
        if let Some(to_and_from) = self
            .set()
            .binary_search_by_key(&perm.perm, elem, |(elem, _)| elem)
        {
            &perm.perm[to_and_from.1.to].0
        } else {
            elem
        }
    }

    fn preimage_ref<'a>(
        self: &Arc<Self>,
        perm: &'a <Self as SetSignature>::Elem,
        elem: &'a Set::Elem,
    ) -> &'a Set::Elem {
        debug_assert!(self.is_element(perm));
        if let Some(to_and_from) = self
            .set()
            .binary_search_by_key(&perm.perm, elem, |(elem, _)| elem)
        {
            &perm.perm[to_and_from.1.from].0
        } else {
            elem
        }
    }
}

impl<Set: OrdSignature> PermutationsSignature<Set> for FinitelySupportedPermutationsStructure<Set> {
    fn set(self: &Arc<Self>) -> &Arc<Set> {
        &self.set
    }

    fn new_cycle(self: &Arc<Self>, cycle: Vec<Set::Elem>) -> Result<Self::Elem, ()> {
        let sorted_enumerated_cycle = self
            .set()
            .sort_by_key(cycle.iter().enumerate().collect(), &|item| item.1);
        if !self
            .set()
            .is_sorted_and_unique_by_key(&sorted_enumerated_cycle, |item| item.1)
        {
            return Err(());
        }
        let n: usize = sorted_enumerated_cycle.len();
        let enumerated_sorted_enumerated_cycle = sorted_enumerated_cycle
            .into_iter()
            .enumerate()
            .collect::<Vec<_>>();
        let perm = enumerated_sorted_enumerated_cycle
            .iter()
            .map(|(_, (elem_idx, elem))| {
                ((*elem).clone(), {
                    let elem_idx_from = elem_idx.checked_sub(1).unwrap_or(n - 1); // (i-1) % n
                    let elem_idx_to = if elem_idx + 1 == n { 0 } else { elem_idx + 1 }; // (i+1) % n
                    FromAndTo {
                        from: self
                            .set()
                            .binary_search_by_key(
                                &enumerated_sorted_enumerated_cycle,
                                &cycle[elem_idx_from],
                                |item| item.1.1,
                            )
                            .unwrap()
                            .0,
                        to: self
                            .set()
                            .binary_search_by_key(
                                &enumerated_sorted_enumerated_cycle,
                                &cycle[elem_idx_to],
                                |item| item.1.1,
                            )
                            .unwrap()
                            .0,
                    }
                })
            })
            .collect();
        let s = FinitelySupportedPermutation { perm };
        debug_assert!(self.is_element(&s));
        Ok(s)
    }

    fn new_perm(
        self: &Arc<Self>,
        perm: Vec<(impl Borrow<Set::Elem>, impl Borrow<Set::Elem>)>,
    ) -> Result<Self::Elem, ()> {
        let perm_sorted_froms = self
            .set()
            .sort_by_key(perm.iter().collect(), &|(from, _)| from.borrow());
        let perm_sorted_tos = self
            .set()
            .sort_by_key(perm.iter().collect(), &|(_, to)| to.borrow());

        if !self
            .set()
            .is_sorted_and_unique_by_key(&perm_sorted_froms, |(from, _)| from.borrow())
        {
            return Err(());
        }

        if !self
            .set()
            .is_sorted_and_unique_by_key(&perm_sorted_tos, |(_, to)| to.borrow())
        {
            return Err(());
        }

        for merged in self.set().merge_sorted_and_unique_by_key(
            perm_sorted_froms
                .iter()
                .map(|(from, _)| from.borrow())
                .collect(),
            perm_sorted_tos.iter().map(|(_, to)| to.borrow()).collect(),
            |item| item,
        ) {
            match merged {
                MergedUniqueSource::First(_) | MergedUniqueSource::Second(_) => {
                    return Err(());
                }
                MergedUniqueSource::Both(from_elem, to_elem) => {
                    debug_assert!(self.set().equal(from_elem, to_elem));
                }
            }
        }

        // at this point we know the input is a valid permutation
        let elems_sorted = perm_sorted_tos
            .iter()
            .filter(|(from, to)| !self.set().equal(from.borrow(), to.borrow()))
            .map(|(_, to)| to.borrow())
            .collect::<Vec<_>>(); // could just as well use froms_sorted here

        let s = FinitelySupportedPermutation {
            perm: elems_sorted
                .iter()
                .map(|elem| {
                    ((*elem).borrow().clone(), {
                        FromAndTo {
                            from: self
                                .set()
                                .binary_search_index(
                                    &elems_sorted,
                                    self.set()
                                        .binary_search_by_key(
                                            &perm_sorted_tos,
                                            elem.borrow(),
                                            |item| item.1.borrow(),
                                        )
                                        .unwrap()
                                        .0
                                        .borrow(),
                                )
                                .unwrap(),
                            to: self
                                .set()
                                .binary_search_index(
                                    &elems_sorted,
                                    self.set()
                                        .binary_search_by_key(
                                            &perm_sorted_froms,
                                            elem.borrow(),
                                            |item| item.0.borrow(),
                                        )
                                        .unwrap()
                                        .1
                                        .borrow(),
                                )
                                .unwrap(),
                        }
                    })
                })
                .collect(),
        };
        debug_assert!(self.is_element(&s));
        Ok(s)
    }

    fn support(self: &Arc<Self>, perm: &Self::Elem) -> Vec<Set::Elem> {
        debug_assert!(self.is_element(perm));
        perm.perm.iter().map(|(elem, _)| elem.clone()).collect()
    }

    fn support_size(self: &Arc<Self>, perm: &Self::Elem) -> usize {
        debug_assert!(self.is_element(perm));
        perm.perm.len()
    }

    fn image(self: &Arc<Self>, perm: &Self::Elem, elem: &Set::Elem) -> Set::Elem {
        self.image_ref(perm, elem).clone()
    }

    fn preimage(self: &Arc<Self>, perm: &Self::Elem, elem: &Set::Elem) -> Set::Elem {
        self.preimage_ref(perm, elem).clone()
    }

    fn disjoint_cycles(self: &Arc<Self>, perm: &Self::Elem) -> Vec<Vec<Set::Elem>> {
        debug_assert!(self.is_element(perm));
        // vector of pairs of moved elements and whether they have been accounted for
        let mut elems_todo = perm
            .perm
            .iter()
            .map(|(elem, _)| (elem, false))
            .collect::<Vec<_>>();
        debug_assert!(
            self.set()
                .is_sorted_and_unique_by_key(&elems_todo, |item| item.0)
        );

        let mut cycles = vec![];
        for i in 0..elems_todo.len() {
            if elems_todo[i].1 {
                continue;
            }
            let mut cycle = vec![];
            let mut j = i;
            while !elems_todo[j].1 {
                cycle.push(elems_todo[j].0.clone());
                elems_todo[j].1 = true;
                j = self
                    .set()
                    .binary_search_by_key(&perm.perm, elems_todo[j].0, |item| &item.0)
                    .unwrap()
                    .1
                    .to;
            }
            cycles.push(cycle);
        }
        cycles
    }
}

impl<Set: OrdSignature> CompositionSignature for FinitelySupportedPermutationsStructure<Set> {
    fn compose(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));

        // filter out elements which end up fixed after the composition
        let unfixed_elems = self
            .set()
            .merge_sorted_and_unique_by_key(a.perm.clone(), b.perm.clone(), |item| &item.0)
            .filter(|merged_item| match merged_item {
                MergedUniqueSource::Second(b_item) | MergedUniqueSource::Both(_, b_item) => {
                    if let Some(a_item) =
                        self.set()
                            .binary_search_by_key(&a.perm, &b.perm[b_item.1.to].0, |item| &item.0)
                        && let Some(b_item_2) = self.set().binary_search_by_key(
                            &b.perm,
                            &a.perm[a_item.1.to].0,
                            |item| &item.0,
                        )
                    {
                        !self.set().equal(&b_item.0, &b_item_2.0)
                    } else {
                        true
                    }
                }
                MergedUniqueSource::First(_) => true,
            })
            .map(|merged_item| match merged_item {
                MergedUniqueSource::First(a_item) => a_item.0,
                MergedUniqueSource::Second(b_item) => b_item.0,
                MergedUniqueSource::Both(a_item, b_item) => {
                    debug_assert!(self.set().equal(&a_item.0, &b_item.0));
                    a_item.0
                }
            })
            .collect::<Vec<_>>();

        // compute the composition
        let s = FinitelySupportedPermutation {
            perm: unfixed_elems
                .iter()
                .map(|elem| {
                    let elem_image = self.image_ref(a, self.image_ref(b, elem));
                    let elem_preimage = self.preimage_ref(b, self.preimage_ref(a, elem));
                    (
                        elem.clone(),
                        FromAndTo {
                            from: self
                                .set()
                                .binary_search_index(&unfixed_elems, elem_preimage)
                                .unwrap(),
                            to: self
                                .set()
                                .binary_search_index(&unfixed_elems, elem_image)
                                .unwrap(),
                        },
                    )
                })
                .collect(),
        };
        debug_assert!(self.is_element(&s));
        s
    }
}

impl<Set: OrdSignature> AssociativeCompositionSignature
    for FinitelySupportedPermutationsStructure<Set>
{
}

impl<Set: OrdSignature> LeftCancellativeCompositionSignature
    for FinitelySupportedPermutationsStructure<Set>
{
    fn try_left_difference(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.compose(&self.inverse(b), a))
    }
}

impl<Set: OrdSignature> RightCancellativeCompositionSignature
    for FinitelySupportedPermutationsStructure<Set>
{
    fn try_right_difference(
        self: &Arc<Self>,
        a: &Self::Elem,
        b: &Self::Elem,
    ) -> Option<Self::Elem> {
        Some(self.compose(a, &self.inverse(b)))
    }
}

impl<Set: OrdSignature> IdentitySignature for FinitelySupportedPermutationsStructure<Set> {
    fn identity(self: &Arc<Self>) -> Self::Elem {
        FinitelySupportedPermutation::identity()
    }
}

impl<Set: OrdSignature> MonoidSignature for FinitelySupportedPermutationsStructure<Set> {}

impl<Set: OrdSignature> TryLeftInverseSignature for FinitelySupportedPermutationsStructure<Set> {
    fn try_left_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<Set: OrdSignature> TryRightInverseSignature for FinitelySupportedPermutationsStructure<Set> {
    fn try_right_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<Set: OrdSignature> TryInverseSignature for FinitelySupportedPermutationsStructure<Set> {
    fn try_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<Set: OrdSignature> GroupSignature for FinitelySupportedPermutationsStructure<Set> {
    fn inverse(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        a.clone().inverse()
    }
}

impl<Set: OrdSignature + FiniteSetSignature> CountableSetSignature
    for FinitelySupportedPermutationsStructure<Set>
{
    fn generate_all_elements(self: Arc<Self>) -> impl Iterator<Item = Self::Elem> {
        let all_elems = self.set().list_all_elements();
        let n = all_elems.len();
        (0..n)
            .permutations(n)
            .map(|perm| {
                self.new_perm(
                    perm.into_iter()
                        .enumerate()
                        .map(|(from, to)| (&all_elems[from], &all_elems[to]))
                        .collect(),
                )
                .unwrap()
            })
            .collect::<Vec<_>>()
            .into_iter()
    }
}

impl<Set: OrdSignature + FiniteSetSignature> FiniteSetSignature
    for FinitelySupportedPermutationsStructure<Set>
{
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PermutationActionStructure<Set: OrdSignature> {
    set: Arc<Set>,
}

impl<Set: OrdSignature> PermutationActionStructure<Set> {
    fn new(set: Arc<Set>) -> Arc<Self> {
        Self { set }.into()
    }
}

pub trait SetToFinitelySupportedPermutationAction: OrdSignature {
    fn finitely_supported_permutation_action(
        self: Arc<Self>,
    ) -> Arc<impl LeftGroupActionSignature<FinitelySupportedPermutationsStructure<Self>, Self>>
    {
        PermutationActionStructure::new(self)
    }
}
impl<Set: OrdSignature> SetToFinitelySupportedPermutationAction for Set {}

impl<Set: OrdSignature> Signature for PermutationActionStructure<Set> {}

impl<Set: OrdSignature> LeftGroupActionSignature<FinitelySupportedPermutationsStructure<Set>, Set>
    for PermutationActionStructure<Set>
{
    fn group(self: &Arc<Self>) -> Arc<FinitelySupportedPermutationsStructure<Set>> {
        self.set.finitely_supported_permutations()
    }

    fn set(self: &Arc<Self>) -> Arc<Set> {
        self.set.clone()
    }

    fn apply(
        self: &Arc<Self>,
        g: &FinitelySupportedPermutation<Set::Elem>,
        x: &Set::Elem,
    ) -> Set::Elem {
        self.group().image(g, x)
    }
}

#[cfg(test)]
mod partition_tests {
    use super::*;

    #[test]
    fn test_image_and_preimage() {
        let x = FinitelySupportedPermutation::<i32>::new_cycle(vec![1, 2, 3, 4]).unwrap();
        assert_eq!(x.image(&3), 4);
        assert_eq!(x.preimage(&3), 2);
    }

    #[test]
    fn test_new_perm() {
        let x =
            FinitelySupportedPermutation::<i32>::new_perm(vec![(1, 2), (2, 3), (3, 1)]).unwrap();
        let y = FinitelySupportedPermutation::<i32>::new_cycle(vec![1, 2, 3]).unwrap();
        assert!(x.equal(&y));
        println!("{:?}", x.disjoint_cycles());
    }

    #[test]
    fn test_composition_and_equal() {
        let a = FinitelySupportedPermutation::<i32>::new_cycle(vec![1, 2, 3, 4]).unwrap();
        let b = FinitelySupportedPermutation::<i32>::new_cycle(vec![3, 4, 5, 6]).unwrap();

        let c = FinitelySupportedPermutation::<i32>::new_cycle(vec![1, 2, 3]).unwrap();
        let d = FinitelySupportedPermutation::<i32>::new_cycle(vec![4, 5, 6]).unwrap();

        assert!(a.compose(&b).equal(&c.compose(&d)));
    }

    #[test]
    fn test_composition_and_equal_2() {
        let a = FinitelySupportedPermutation::<i32>::new_perm(vec![(0, 0), (1, 2), (2, 1), (3, 3)])
            .unwrap();
        let b = FinitelySupportedPermutation::<i32>::new_perm(vec![(0, 1), (1, 2), (2, 0), (3, 3)])
            .unwrap();
        let c =
            FinitelySupportedPermutation::<i32>::new_perm(vec![(0, 2), (1, 1), (2, 0)]).unwrap();

        println!("a = {:?}", a.disjoint_cycles());
        println!("b = {:?}", b.disjoint_cycles());
        println!("ab = {:?}", a.compose(&b).disjoint_cycles());
        println!("c = {:?}", c.disjoint_cycles());

        assert!(a.compose(&b).equal(&c));
    }
}
