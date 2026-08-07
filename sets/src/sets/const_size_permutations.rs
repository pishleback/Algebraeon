use algebraeon_structures::*;
use itertools::Itertools;
use std::borrow::Borrow;
use std::collections::HashSet;
use std::marker::PhantomData;
use std::sync::Arc;

#[derive(Debug, Clone)]
pub struct ConstSizePermutation<const N: usize, Elem> {
    _elem: PhantomData<Elem>,
    forward: [usize; N],
    backward: [usize; N],
}

impl<const N: usize, Elem> PartialEq for ConstSizePermutation<N, Elem> {
    fn eq(&self, other: &Self) -> bool {
        self.forward == other.forward
    }
}

impl<const N: usize, Elem> Eq for ConstSizePermutation<N, Elem> {}

impl<const N: usize, Elem> ConstSizePermutation<N, Elem> {
    pub fn validate(&self) -> Result<(), String> {
        let expected = (0..N).collect::<HashSet<_>>();
        let forward = self.forward.iter().cloned().collect();
        if expected != forward {
            return Err("`forward` does not contain all the values 0..N exactly one".to_string());
        }

        let backward = self.backward.iter().cloned().collect();
        if expected != backward {
            return Err("`backward` does not contain all the values 0..N exactly one".to_string());
        }

        for i in 0..N {
            if i != self.forward[self.backward[i]] || i != self.backward[self.forward[i]] {
                return Err("`forward` and `backward` are not mutually inverse".to_string());
            }
        }

        Ok(())
    }

    fn identity() -> Self {
        Self {
            _elem: PhantomData,
            forward: std::array::from_fn(|i| i),
            backward: std::array::from_fn(|i| i),
        }
    }

    fn inverse(self) -> Self {
        Self {
            _elem: PhantomData,
            forward: self.backward,
            backward: self.forward,
        }
    }

    fn disjoint_cycles(&self) -> Vec<Vec<usize>> {
        // vector of pairs of moved elements and whether they have been accounted for
        let mut elems_todo = (0..N).map(|elem| (elem, false)).collect::<Vec<_>>();

        let mut cycles = vec![];
        for i in 0..elems_todo.len() {
            if elems_todo[i].1 {
                continue;
            }
            let mut cycle = vec![];
            let mut j = i;
            while !elems_todo[j].1 {
                cycle.push(elems_todo[j].0);
                elems_todo[j].1 = true;
                j = self.forward[j];
            }
            cycles.push(cycle);
        }
        cycles
    }
}

impl<const N: usize, Elem: MetaType> MetaType for ConstSizePermutation<N, Elem>
where
    Elem::Signature: ConstSizeFiniteSetSignature<N>,
{
    type Signature = ConstSizePermutationsStructure<N, Elem::Signature>;

    fn structure() -> Arc<Self::Signature> {
        ConstSizePermutationsStructure::new(Elem::structure())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstSizePermutationsStructure<const N: usize, Set: ConstSizeFiniteSetSignature<N>> {
    set: Arc<Set>,
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N>> ConstSizePermutationsStructure<N, Set> {
    pub fn new(set: Arc<Set>) -> Arc<Self> {
        Self { set }.into()
    }
}

pub trait SetToConstSizePermutationsStructure<const N: usize>:
    ConstSizeFiniteSetSignature<N>
{
    fn const_size_permutations(self: &Arc<Self>) -> Arc<ConstSizePermutationsStructure<N, Self>> {
        ConstSizePermutationsStructure::new(self.clone())
    }
}
impl<const N: usize, Set: ConstSizeFiniteSetSignature<N>> SetToConstSizePermutationsStructure<N>
    for Set
{
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N>> Signature
    for ConstSizePermutationsStructure<N, Set>
{
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N>> SetSignature
    for ConstSizePermutationsStructure<N, Set>
{
    type Elem = ConstSizePermutation<N, Set::Elem>;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        x.validate()
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N>> EqSignature
    for ConstSizePermutationsStructure<N, Set>
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        a == b
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    PermutationsSignature<Set> for ConstSizePermutationsStructure<N, Set>
{
    fn set(self: &Arc<Self>) -> &Arc<Set> {
        &self.set
    }

    fn new_cycle(self: &Arc<Self>, cycle: Vec<Set::Elem>) -> Result<Self::Elem, ()> {
        let k = cycle.len();
        if k == 0 {
            return Ok(ConstSizePermutation::identity());
        }
        let cycle = cycle
            .into_iter()
            .map(|elem| self.set().element_to_enumeration(&elem).try_into().unwrap());
        let cycle_unique = cycle.unique().collect::<Vec<usize>>();
        if k != cycle_unique.len() {
            return Err(());
        }
        let cycle = cycle_unique;

        let mut forward = std::array::from_fn(|i| i);
        let mut backward = std::array::from_fn(|i| i);
        for i in 0..k {
            let a = cycle[i];
            let b = cycle[(i + 1) % k];
            forward[a] = b;
            backward[b] = a;
        }

        let perm = ConstSizePermutation {
            _elem: PhantomData,
            forward,
            backward,
        };
        #[cfg(debug_assertions)]
        perm.validate().unwrap();
        Ok(perm)
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
        let perm = perm.into_iter().map(|(from, to)| {
            let from_idx: usize = self
                .set()
                .element_to_enumeration(from.borrow())
                .try_into()
                .unwrap();
            let to_idx: usize = self
                .set()
                .element_to_enumeration(to.borrow())
                .try_into()
                .unwrap();
            (from_idx, to_idx)
        });

        let mut forward = std::array::from_fn(|i| i);
        let mut backward = std::array::from_fn(|i| i);
        for (a, b) in perm {
            forward[a] = b;
            backward[b] = a;
        }

        let perm = ConstSizePermutation {
            _elem: PhantomData,
            forward,
            backward,
        };

        debug_assert!(perm.validate().is_ok());
        Ok(perm)
    }

    fn support(self: &Arc<Self>, perm: &Self::Elem) -> Vec<Set::Elem> {
        debug_assert!(self.is_element(perm));
        (0..N)
            .filter(|i| *i != perm.forward[*i])
            .map(|i| {
                self.set()
                    .enumeration_to_element(&Natural::from(i))
                    .unwrap()
            })
            .collect()
    }

    fn support_size(self: &Arc<Self>, perm: &Self::Elem) -> usize {
        debug_assert!(self.is_element(perm));
        self.support(perm).len()
    }

    fn image(self: &Arc<Self>, perm: &Self::Elem, elem: &Set::Elem) -> Set::Elem {
        let num: usize = self.set().element_to_enumeration(elem).try_into().unwrap();
        self.set()
            .enumeration_to_element(&Natural::from(perm.forward[num]))
            .unwrap()
    }

    fn preimage(self: &Arc<Self>, perm: &Self::Elem, elem: &Set::Elem) -> Set::Elem {
        let num: usize = self.set().element_to_enumeration(elem).try_into().unwrap();
        self.set()
            .enumeration_to_element(&Natural::from(perm.backward[num]))
            .unwrap()
    }

    fn disjoint_cycles(self: &Arc<Self>, perm: &Self::Elem) -> Vec<Vec<Set::Elem>> {
        debug_assert!(self.is_element(perm));
        perm.disjoint_cycles()
            .into_iter()
            .map(|cycle| {
                cycle
                    .into_iter()
                    .map(|i| self.set().enumeration_to_element(&i.into()).unwrap())
                    .collect()
            })
            .collect()
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    CompositionSignature for ConstSizePermutationsStructure<N, Set>
{
    fn compose(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        let s = ConstSizePermutation {
            _elem: PhantomData,
            forward: std::array::from_fn(|i| a.forward[b.forward[i]]),
            backward: std::array::from_fn(|i| b.backward[a.backward[i]]),
        };
        debug_assert!(self.is_element(&s));
        s
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    AssociativeCompositionSignature for ConstSizePermutationsStructure<N, Set>
{
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    LeftCancellativeCompositionSignature for ConstSizePermutationsStructure<N, Set>
{
    fn try_left_difference(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.compose(&self.inverse(b), a))
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    RightCancellativeCompositionSignature for ConstSizePermutationsStructure<N, Set>
{
    fn try_right_difference(
        self: &Arc<Self>,
        a: &Self::Elem,
        b: &Self::Elem,
    ) -> Option<Self::Elem> {
        Some(self.compose(a, &self.inverse(b)))
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    IdentitySignature for ConstSizePermutationsStructure<N, Set>
{
    fn identity(self: &Arc<Self>) -> Self::Elem {
        ConstSizePermutation::identity()
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    MonoidSignature for ConstSizePermutationsStructure<N, Set>
{
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    TryLeftInverseSignature for ConstSizePermutationsStructure<N, Set>
{
    fn try_left_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    TryRightInverseSignature for ConstSizePermutationsStructure<N, Set>
{
    fn try_right_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    TryInverseSignature for ConstSizePermutationsStructure<N, Set>
{
    fn try_inverse(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature> GroupSignature
    for ConstSizePermutationsStructure<N, Set>
{
    fn inverse(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        a.clone().inverse()
    }
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    CountableSetSignature for ConstSizePermutationsStructure<N, Set>
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

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    FiniteSetSignature for ConstSizePermutationsStructure<N, Set>
{
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PermutationActionStructure<
    const N: usize,
    Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
> {
    set: Arc<Set>,
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    PermutationActionStructure<N, Set>
{
    fn new(set: Arc<Set>) -> Arc<Self> {
        Self { set }.into()
    }
}

pub trait SetToConstSizePermutationAction<const N: usize>:
    ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature
{
    fn const_size_permutation_action(
        self: Arc<Self>,
    ) -> Arc<impl LeftGroupActionSignature<ConstSizePermutationsStructure<N, Self>, Self>> {
        PermutationActionStructure::new(self)
    }
}
impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    SetToConstSizePermutationAction<N> for Set
{
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature> Signature
    for PermutationActionStructure<N, Set>
{
}

impl<const N: usize, Set: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature>
    LeftGroupActionSignature<ConstSizePermutationsStructure<N, Set>, Set>
    for PermutationActionStructure<N, Set>
{
    fn group(self: &Arc<Self>) -> Arc<ConstSizePermutationsStructure<N, Set>> {
        self.set.const_size_permutations()
    }

    fn set(self: &Arc<Self>) -> Arc<Set> {
        self.set.clone()
    }

    fn apply(self: &Arc<Self>, g: &ConstSizePermutation<N, Set::Elem>, x: &Set::Elem) -> Set::Elem {
        self.group().image(g, x)
    }
}

#[cfg(test)]
mod partition_tests {
    use super::*;
    use crate::sets::ConstSizeEnumeratedFiniteSet;

    #[test]
    fn test_image_and_preimage() {
        let perms = ConstSizeEnumeratedFiniteSet::<4>::structure().const_size_permutations();

        let x = perms
            .new_cycle(vec![
                0.try_into().unwrap(),
                1.try_into().unwrap(),
                2.try_into().unwrap(),
                3.try_into().unwrap(),
            ])
            .unwrap();

        assert_eq!(x.image(&2.try_into().unwrap()), 3.try_into().unwrap());
        assert_eq!(x.preimage(&2.try_into().unwrap()), 1.try_into().unwrap());
    }

    #[test]
    fn test_new_perm() {
        type S = ConstSizeEnumeratedFiniteSet<4>;
        let perms = S::structure().const_size_permutations();

        let x = perms
            .new_perm(vec![
                (S::try_from(1).unwrap(), S::try_from(2).unwrap()),
                (S::try_from(2).unwrap(), S::try_from(3).unwrap()),
                (S::try_from(3).unwrap(), S::try_from(1).unwrap()),
            ])
            .unwrap();
        let y = perms
            .new_cycle(vec![
                S::try_from(1).unwrap(),
                S::try_from(2).unwrap(),
                S::try_from(3).unwrap(),
            ])
            .unwrap();

        assert!(x.equal(&y));
        println!("{:?}", x.disjoint_cycles());
    }

    #[test]
    fn test_composition_and_equal() {
        type S = ConstSizeEnumeratedFiniteSet<7>;
        let perms = S::structure().const_size_permutations();

        let a = perms
            .new_cycle(
                [1, 2, 3, 4]
                    .into_iter()
                    .map(|x| S::try_from(x).unwrap())
                    .collect(),
            )
            .unwrap();
        let b = perms
            .new_cycle(
                [3, 4, 5, 6]
                    .into_iter()
                    .map(|x| S::try_from(x).unwrap())
                    .collect(),
            )
            .unwrap();

        let c = perms
            .new_cycle(
                [1, 2, 3]
                    .into_iter()
                    .map(|x| S::try_from(x).unwrap())
                    .collect(),
            )
            .unwrap();
        let d = perms
            .new_cycle(
                [4, 5, 6]
                    .into_iter()
                    .map(|x| S::try_from(x).unwrap())
                    .collect(),
            )
            .unwrap();

        assert!(a.compose(&b).equal(&c.compose(&d)));
    }

    #[test]
    fn test_composition_and_equal_2() {
        type S = ConstSizeEnumeratedFiniteSet<7>;
        let perms = S::structure().const_size_permutations();

        let a = perms
            .new_perm(
                [(0, 0), (1, 2), (2, 1), (3, 3)]
                    .into_iter()
                    .map(|(x, y)| (S::try_from(x).unwrap(), S::try_from(y).unwrap()))
                    .collect(),
            )
            .unwrap();
        let b = perms
            .new_perm(
                [(0, 1), (1, 2), (2, 0), (3, 3)]
                    .into_iter()
                    .map(|(x, y)| (S::try_from(x).unwrap(), S::try_from(y).unwrap()))
                    .collect(),
            )
            .unwrap();
        let c = perms
            .new_perm(
                [(0, 2), (1, 1), (2, 0)]
                    .into_iter()
                    .map(|(x, y)| (S::try_from(x).unwrap(), S::try_from(y).unwrap()))
                    .collect(),
            )
            .unwrap();

        println!("a = {:?}", a.disjoint_cycles());
        println!("b = {:?}", b.disjoint_cycles());
        println!("ab = {:?}", a.compose(&b).disjoint_cycles());
        println!("c = {:?}", c.disjoint_cycles());

        assert!(a.compose(&b).equal(&c));
    }
}
