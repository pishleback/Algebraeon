use crate::sets::{
    ConstSizePermutationsStructure, FiniteSetToFinitelySupportedPermutationsStructure,
    FinitelySupportedPermutationsStructure, FunctionsStructure,
    SetToConstSizePermutationsStructure,
};
use algebraeon_structures::*;
use itertools::Itertools;
use std::{
    cmp::Ordering,
    fmt::Debug,
    marker::PhantomData,
    ops::{Index, IndexMut},
    sync::Arc,
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Function<const N: usize, DomainElem, RangeElem> {
    _domain: PhantomData<DomainElem>,
    images: [RangeElem; N],
}

impl<const N: usize, DomainElem: MetaType, RangeElem> Function<N, DomainElem, RangeElem> {
    pub fn new(f: impl FnMut(DomainElem) -> RangeElem) -> Self
    where
        DomainElem::Signature: OrderedFiniteSetSignature + ConstSizeFiniteSetSignature<N>,
    {
        Self {
            _domain: PhantomData,
            images: DomainElem::structure()
                .list_all_elements_ordered()
                .into_iter()
                .map(f)
                .collect::<Vec<_>>()
                .try_into()
                .map_err(|_| ())
                .unwrap(),
        }
    }
}

impl<const N: usize, DomainElem, RangeElem> Function<N, DomainElem, RangeElem> {
    pub fn iter(&self) -> std::slice::Iter<'_, RangeElem> {
        self.images.iter()
    }

    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, RangeElem> {
        self.images.iter_mut()
    }

    pub fn map<NewRangeElem>(
        self,
        f: impl FnMut(RangeElem) -> NewRangeElem,
    ) -> Function<N, DomainElem, NewRangeElem> {
        Function {
            _domain: PhantomData,
            images: self.images.map(f),
        }
    }
}

impl<const N: usize, DomainElem, RangeElem> IntoIterator for Function<N, DomainElem, RangeElem> {
    type Item = RangeElem;
    type IntoIter = std::array::IntoIter<RangeElem, N>;

    fn into_iter(self) -> Self::IntoIter {
        self.images.into_iter()
    }
}

impl<const N: usize, DomainElem, RangeElem> Index<usize> for Function<N, DomainElem, RangeElem> {
    type Output = RangeElem;

    fn index(&self, index: usize) -> &Self::Output {
        &self.images[index]
    }
}

impl<const N: usize, DomainElem, RangeElem> IndexMut<usize> for Function<N, DomainElem, RangeElem> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.images[index]
    }
}

impl<const N: usize, DomainElem, RangeElem> From<[RangeElem; N]>
    for Function<N, DomainElem, RangeElem>
{
    fn from(images: [RangeElem; N]) -> Self {
        Self {
            _domain: PhantomData,
            images,
        }
    }
}

impl<const N: usize, DomainElem, RangeElem> From<Function<N, DomainElem, RangeElem>>
    for [RangeElem; N]
{
    fn from(func: Function<N, DomainElem, RangeElem>) -> Self {
        func.images
    }
}

impl<const N: usize, DomainElem, RangeElem: Clone> From<&Function<N, DomainElem, RangeElem>>
    for [RangeElem; N]
{
    fn from(func: &Function<N, DomainElem, RangeElem>) -> Self {
        func.images.clone()
    }
}

impl<const N: usize, DomainElem, RangeElem> TryFrom<Vec<RangeElem>>
    for Function<N, DomainElem, RangeElem>
{
    type Error = ();

    fn try_from(images: Vec<RangeElem>) -> Result<Self, Self::Error> {
        if let Ok(images) = TryInto::<[RangeElem; N]>::try_into(images) {
            Ok(images.into())
        } else {
            Err(())
        }
    }
}

impl<const N: usize, DomainElem, RangeElem> From<Function<N, DomainElem, RangeElem>>
    for Vec<RangeElem>
{
    fn from(func: Function<N, DomainElem, RangeElem>) -> Self {
        func.images.into()
    }
}

impl<const N: usize, DomainElem, RangeElem: Clone> From<&Function<N, DomainElem, RangeElem>>
    for Vec<RangeElem>
{
    fn from(func: &Function<N, DomainElem, RangeElem>) -> Self {
        func.images.clone().into()
    }
}

impl<const N: usize, DomainElem: MetaType, RangeElem: MetaType> MetaType
    for Function<N, DomainElem, RangeElem>
where
    DomainElem::Signature: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
{
    type Signature = ConstSizeFunctionsStructure<N, DomainElem::Signature, RangeElem::Signature>;

    fn structure() -> Arc<Self::Signature> {
        ConstSizeFunctionsStructure::new(DomainElem::structure(), RangeElem::structure())
    }
}

/// Represent all functions from `domain` to `range`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeFunctionsStructure<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N>,
    Range: SetSignature,
> {
    domain: Arc<Domain>,
    range: Arc<Range>,
}

impl<const N: usize, Domain: ConstSizeFiniteSetSignature<N>, Range: SetSignature>
    ConstSizeFunctionsStructure<N, Domain, Range>
{
    pub fn new(domain: Arc<Domain>, range: Arc<Range>) -> Arc<Self> {
        Self { domain, range }.into()
    }

    pub fn forget_const(self: &Arc<Self>) -> Arc<FunctionsStructure<Domain, Range>> {
        FunctionsStructure::new(self.domain.clone(), self.range.clone())
    }
}

pub trait SetToConstSizeFunctionsToSignature: SetSignature {
    fn const_size_functions_to<const N: usize, Range: SetSignature>(
        self: &Arc<Self>,
        range: &Arc<Range>,
    ) -> Arc<ConstSizeFunctionsStructure<N, Self, Range>>
    where
        Self: ConstSizeFiniteSetSignature<N>,
    {
        ConstSizeFunctionsStructure::new(self.clone(), range.clone())
    }
}
impl<Set: SetSignature> SetToConstSizeFunctionsToSignature for Set {}

pub trait SetToConstSizeFunctionsFromSignature: SetSignature {
    fn functions_from<const N: usize, Domain: ConstSizeFiniteSetSignature<N>>(
        self: &Arc<Self>,
        domain: &Arc<Domain>,
    ) -> Arc<ConstSizeFunctionsStructure<N, Domain, Self>> {
        ConstSizeFunctionsStructure::new(domain.clone(), self.clone())
    }
}
impl<Set: SetSignature> SetToConstSizeFunctionsFromSignature for Set {}

impl<const N: usize, Domain: ConstSizeFiniteSetSignature<N>, Range: SetSignature>
    ConstSizeFunctionsStructure<N, Domain, Range>
{
    pub fn domain(&self) -> &Arc<Domain> {
        &self.domain
    }

    pub fn range(&self) -> &Arc<Range> {
        &self.range
    }
}

impl<const N: usize, Domain: ConstSizeFiniteSetSignature<N>, Range: SetSignature> Signature
    for ConstSizeFunctionsStructure<N, Domain, Range>
{
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: SetSignature,
> FunctionsSignature<Domain, Range> for ConstSizeFunctionsStructure<N, Domain, Range>
{
    fn function(
        self: &Arc<Self>,
        f: impl Fn(&Domain::Elem) -> Range::Elem,
    ) -> Option<Function<N, Domain::Elem, Range::Elem>> {
        let s: [_; N] = self
            .domain()
            .list_all_elements_ordered()
            .iter()
            .map(f)
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        for y in &s {
            if !self.range().is_element(y) {
                return None;
            }
        }
        Some(s.into())
    }

    fn image<'a>(
        self: &Arc<Self>,
        f: &'a <Self as SetSignature>::Elem,
        x: &Domain::Elem,
    ) -> &'a Range::Elem {
        debug_assert!(self.is_element(f));
        debug_assert!(self.domain().is_element(x));
        &f[TryInto::<usize>::try_into(self.domain().element_to_enumeration(x)).unwrap()]
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: SetSignature,
> SetSignature for ConstSizeFunctionsStructure<N, Domain, Range>
{
    type Elem = Function<N, Domain::Elem, Range::Elem>;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        for y in &x.images {
            self.range().validate_element(y)?;
        }
        Ok(())
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: EqSignature,
> EqSignature for ConstSizeFunctionsStructure<N, Domain, Range>
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        (0..N).all(|i| self.range().equal(&a[i], &b[i]))
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: OrdSignature,
> PartialOrdSignature for ConstSizeFunctionsStructure<N, Domain, Range>
{
    fn partial_cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: OrdSignature,
> OrdSignature for ConstSizeFunctionsStructure<N, Domain, Range>
{
    fn cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        for i in (0..N).rev() {
            match self.range().cmp(&a[i], &b[i]) {
                Ordering::Less => {
                    return Ordering::Less;
                }
                Ordering::Greater => {
                    return Ordering::Greater;
                }
                Ordering::Equal => {}
            }
        }
        Ordering::Equal
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: OrderedFiniteSetSignature,
> CountableSetSignature for ConstSizeFunctionsStructure<N, Domain, Range>
{
    fn generate_all_elements(self: Arc<Self>) -> impl Iterator<Item = Self::Elem> {
        let n: usize = self.domain().size().try_into().unwrap_or(usize::MAX);
        (0..n)
            .map(|_| self.range().list_all_elements())
            .multi_cartesian_product()
            .map(|f| f.try_into().unwrap())
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: OrderedFiniteSetSignature,
> FiniteSetSignature for ConstSizeFunctionsStructure<N, Domain, Range>
{
    fn size(self: &Arc<Self>) -> Natural {
        self.range().size().pow(&self.domain().size())
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: OrderedFiniteSetSignature,
> OrderedFiniteSetSignature for ConstSizeFunctionsStructure<N, Domain, Range>
{
    fn list_all_elements_ordered(self: &Arc<Self>) -> Vec<Self::Elem> {
        let n: usize = self.domain().size().try_into().unwrap_or(usize::MAX);
        (0..n)
            .map(|_| {
                self.range()
                    .list_all_elements()
                    .into_iter()
                    .rev()
                    .collect::<Vec<_>>()
            })
            .multi_cartesian_product()
            .collect::<Vec<_>>()
            .into_iter()
            .rev()
            .map(|images| {
                images
                    .into_iter()
                    .rev()
                    .collect::<Vec<_>>()
                    .try_into()
                    .unwrap()
            })
            .collect()
    }

    fn element_to_enumeration(self: &Arc<Self>, elem: &Self::Elem) -> Natural {
        debug_assert!(self.is_element(elem));
        let r = self.range().size();
        debug_assert_eq!(self.range().size(), (&r).into());
        let mut num = Natural::ZERO;
        for i in (0..N).rev() {
            num = num * &r + self.range().element_to_enumeration(&elem[i]);
        }
        num
    }

    fn enumeration_to_element(self: &Arc<Self>, num: &Natural) -> Option<Self::Elem> {
        if *num >= self.size() {
            return None;
        }
        let d: usize = self.domain().size().try_into().unwrap();
        let mut f = vec![];
        let r = self.range().size();
        let mut n = num.clone();
        let mut k;
        for _ in 0..d {
            (n, k) = n.div_mod(&r);
            f.push(self.range().enumeration_to_element(&k).unwrap())
        }
        let f = f.try_into().unwrap();
        debug_assert!(self.is_element(&f));
        Some(f)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RightPermutationActionOnConstSizeFunctionsStructure<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: SetSignature,
    DomainPerms: PermutationsSignature<Domain>,
> {
    functions: Arc<ConstSizeFunctionsStructure<N, Domain, Range>>,
    domain_perms: Arc<DomainPerms>,
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: SetSignature,
    DomainPerms: PermutationsSignature<Domain>,
> RightPermutationActionOnConstSizeFunctionsStructure<N, Domain, Range, DomainPerms>
{
    fn new(
        functions: Arc<ConstSizeFunctionsStructure<N, Domain, Range>>,
        domain_perms: Arc<DomainPerms>,
    ) -> Arc<Self> {
        Self {
            functions,
            domain_perms,
        }
        .into()
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: SetSignature,
    DomainPerms: PermutationsSignature<Domain>,
> Signature for RightPermutationActionOnConstSizeFunctionsStructure<N, Domain, Range, DomainPerms>
{
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: SetSignature,
> ConstSizeFunctionsStructure<N, Domain, Range>
{
    pub fn domain_precomposition_finitely_supported_permutation_action(
        self: &Arc<Self>,
    ) -> Arc<impl RightGroupActionSignature<Self, FinitelySupportedPermutationsStructure<Domain>>>
    {
        RightPermutationActionOnConstSizeFunctionsStructure::new(
            self.clone(),
            self.domain().permutations(),
        )
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: SetSignature,
> ConstSizeFunctionsStructure<N, Domain, Range>
{
    pub fn output_finitely_supported_permutation_action(
        self: &Arc<Self>,
    ) -> Arc<impl LeftGroupActionSignature<FinitelySupportedPermutationsStructure<Domain>, Self>>
    {
        self.domain_precomposition_finitely_supported_permutation_action()
            .opposite()
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: SetSignature,
> ConstSizeFunctionsStructure<N, Domain, Range>
{
    pub fn domain_precomposition_const_size_permutation_action(
        self: &Arc<Self>,
    ) -> Arc<impl RightGroupActionSignature<Self, ConstSizePermutationsStructure<N, Domain>>> {
        RightPermutationActionOnConstSizeFunctionsStructure::new(
            self.clone(),
            self.domain().const_size_permutations(),
        )
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: SetSignature,
> ConstSizeFunctionsStructure<N, Domain, Range>
{
    pub fn output_const_size_permutation_action(
        self: &Arc<Self>,
    ) -> Arc<impl LeftGroupActionSignature<ConstSizePermutationsStructure<N, Domain>, Self>> {
        self.domain_precomposition_const_size_permutation_action()
            .opposite()
    }
}

/// Sym(D) has a right action on Fun(D -> R) by composition on the right
impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: SetSignature,
    DomainPerms: PermutationsSignature<Domain>,
> RightGroupActionSignature<ConstSizeFunctionsStructure<N, Domain, Range>, DomainPerms>
    for RightPermutationActionOnConstSizeFunctionsStructure<N, Domain, Range, DomainPerms>
{
    fn group(self: &Arc<Self>) -> &Arc<DomainPerms> {
        &self.domain_perms
    }

    fn set(self: &Arc<Self>) -> &Arc<ConstSizeFunctionsStructure<N, Domain, Range>> {
        &self.functions
    }

    fn apply(
        self: &Arc<Self>,
        g: &<DomainPerms>::Elem,
        f: &<ConstSizeFunctionsStructure<N, Domain, Range> as SetSignature>::Elem,
    ) -> <ConstSizeFunctionsStructure<N, Domain, Range> as SetSignature>::Elem {
        self.functions
            .function(|x| {
                self.functions
                    .image(f, &self.domain_perms.image(g, x))
                    .clone()
            })
            .unwrap()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LeftPermutationActionOnConstSizeFunctionsStructure<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N>,
    Range: OrderedFiniteSetSignature,
    RangePerms: PermutationsSignature<Range>,
> {
    functions: Arc<ConstSizeFunctionsStructure<N, Domain, Range>>,
    range_perms: Arc<RangePerms>,
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N>,
    Range: OrderedFiniteSetSignature,
    RangePerms: PermutationsSignature<Range>,
> LeftPermutationActionOnConstSizeFunctionsStructure<N, Domain, Range, RangePerms>
{
    fn new(
        functions: Arc<ConstSizeFunctionsStructure<N, Domain, Range>>,
        range_perms: Arc<RangePerms>,
    ) -> Arc<Self> {
        Self {
            functions,
            range_perms,
        }
        .into()
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N>,
    Range: OrderedFiniteSetSignature,
    RangePerms: PermutationsSignature<Range>,
> Signature for LeftPermutationActionOnConstSizeFunctionsStructure<N, Domain, Range, RangePerms>
{
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: OrderedFiniteSetSignature,
> ConstSizeFunctionsStructure<N, Domain, Range>
{
    pub fn range_finitely_supported_permutation_action(
        self: &Arc<Self>,
    ) -> Arc<impl LeftGroupActionSignature<FinitelySupportedPermutationsStructure<Range>, Self>>
    {
        LeftPermutationActionOnConstSizeFunctionsStructure::new(
            self.clone(),
            self.range().permutations(),
        )
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: OrderedFiniteSetSignature,
> ConstSizeFunctionsStructure<N, Domain, Range>
{
    pub fn range_const_size_permutation_action<const M: usize>(
        self: &Arc<Self>,
    ) -> Arc<impl LeftGroupActionSignature<ConstSizePermutationsStructure<M, Range>, Self>>
    where
        Range: ConstSizeFiniteSetSignature<M>,
    {
        LeftPermutationActionOnConstSizeFunctionsStructure::new(
            self.clone(),
            self.range().const_size_permutations(),
        )
    }
}

/// Sym(R) has a left action on Fun(D -> R) by composition on the left
impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + OrderedFiniteSetSignature,
    Range: OrderedFiniteSetSignature,
    RangePerms: PermutationsSignature<Range>,
> LeftGroupActionSignature<RangePerms, ConstSizeFunctionsStructure<N, Domain, Range>>
    for LeftPermutationActionOnConstSizeFunctionsStructure<N, Domain, Range, RangePerms>
{
    fn group(self: &Arc<Self>) -> &Arc<RangePerms> {
        &self.range_perms
    }

    fn set(self: &Arc<Self>) -> &Arc<ConstSizeFunctionsStructure<N, Domain, Range>> {
        &self.functions
    }

    fn apply(
        self: &Arc<Self>,
        g: &<RangePerms>::Elem,
        f: &<ConstSizeFunctionsStructure<N, Domain, Range> as SetSignature>::Elem,
    ) -> <ConstSizeFunctionsStructure<N, Domain, Range> as SetSignature>::Elem {
        self.functions
            .function(|x| self.range_perms.image(g, self.functions.image(f, x)))
            .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sets::FiniteSetToFinitelySupportedPermutationsStructure;

    #[test]
    fn enumerate() {
        let set_a = i32::structure().const_size_finite_subset([1, 2, 3, 4, 5]);
        let set_b = i32::structure().const_size_finite_subset([1, 2, 3]);
        let fns = set_a.const_size_functions_to(&set_b);
        algebraeon_structures::assert_enumerated_ord_finite_set!(fns, 243);
    }

    #[test]
    fn test_image() {
        let set_a = i32::structure().const_size_finite_subset([1, 2, 3, 4, 5]);
        let set_b = i32::structure().const_size_finite_subset([1, 2, 3]);
        let fns = set_a.const_size_functions_to(&set_b);
        let f = fns
            .function(|i| match i {
                1 => 1,
                2 => 2,
                3 => 3,
                4 => 2,
                5 => 1,
                _ => unreachable!(),
            })
            .unwrap();
        debug_assert!(set_b.equal(fns.image(&f, &1), &1));
        debug_assert!(set_b.equal(fns.image(&f, &2), &2));
        debug_assert!(set_b.equal(fns.image(&f, &3), &3));
        debug_assert!(set_b.equal(fns.image(&f, &4), &2));
        debug_assert!(set_b.equal(fns.image(&f, &5), &1));
    }

    #[test]
    fn test_permutation_actions() {
        let set_a = i32::structure().const_size_finite_subset([1, 2, 3, 4, 5]);
        let set_a_perms = set_a.permutations();
        let set_b = i32::structure().const_size_finite_subset([1, 2, 3]);
        let set_b_perms = set_b.permutations();
        let fns = set_a.const_size_functions_to(&set_b);

        let x = fns
            .function(|i| match i {
                1 => 1,
                2 => 2,
                3 => 3,
                4 => 2,
                5 => 1,
                _ => unreachable!(),
            })
            .unwrap();

        assert!(
            fns.equal(
                &fns.domain_precomposition_finitely_supported_permutation_action()
                    .apply(&set_a_perms.new_cycle(vec![1, 2, 3, 4, 5]).unwrap(), &x),
                &fns.function(|i| match i {
                    1 => 2,
                    2 => 3,
                    3 => 2,
                    4 => 1,
                    5 => 1,
                    _ => unreachable!(),
                })
                .unwrap()
            )
        );

        assert!(
            fns.equal(
                &fns.range_finitely_supported_permutation_action()
                    .apply(&set_b_perms.new_cycle(vec![1, 2, 3]).unwrap(), &x),
                &fns.function(|i| match i {
                    1 => 2,
                    2 => 3,
                    3 => 1,
                    4 => 3,
                    5 => 2,
                    _ => unreachable!(),
                })
                .unwrap()
            )
        );
    }
}
