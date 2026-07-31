use crate::sets::{
    FiniteSetToFinitelySupportedPermutationsStructure, FinitelySupportedPermutationsStructure,
    FunctionsStructure,
};
use algebraeon_structures::*;
use itertools::Itertools;
use std::{
    cmp::Ordering,
    fmt::Debug,
    marker::PhantomData,
    ops::{Index, IndexMut},
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Function<const N: usize, DomainElem, RangeElem> {
    _domain: PhantomData<DomainElem>,
    images: [RangeElem; N],
}

impl<const N: usize, DomainElem: MetaType, RangeElem> Function<N, DomainElem, RangeElem> {
    pub fn new(f: impl FnMut(DomainElem) -> RangeElem) -> Self
    where
        DomainElem::Signature: EnumeratedOrdFiniteSetSignature + ConstSizeFiniteSetSignature<N>,
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
    DomainElem::Signature: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
{
    type Signature = ConstSizeFunctionsStructure<
        N,
        DomainElem::Signature,
        DomainElem::Signature,
        RangeElem::Signature,
        RangeElem::Signature,
    >;

    fn structure() -> Self::Signature {
        ConstSizeFunctionsStructure::new(DomainElem::structure(), RangeElem::structure())
    }
}

/// Represent all functions from `domain` to `range`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConstSizeFunctionsStructure<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N>,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
> {
    _domain: PhantomData<Domain>,
    domain: DomainB,
    _range: PhantomData<Range>,
    range: RangeB,
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N>,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
> ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    pub fn new(domain: DomainB, range: RangeB) -> Self {
        Self {
            _domain: PhantomData,
            domain,
            _range: PhantomData,
            range,
        }
    }

    pub fn into_domain_and_range(self) -> (DomainB, RangeB) {
        (self.domain, self.range)
    }

    pub fn forget_const(&self) -> FunctionsStructure<Domain, &Domain, Range, &Range> {
        FunctionsStructure::new(self.domain.borrow(), self.range.borrow())
    }

    pub fn into_forget_const(self) -> FunctionsStructure<Domain, DomainB, Range, RangeB> {
        FunctionsStructure::new(self.domain, self.range)
    }
}

pub trait SetToConstSizeFunctionsToSignature: SetSignature {
    fn into_const_size_functions_to<const N: usize, Range: SetSignature>(
        self,
        range: Range,
    ) -> ConstSizeFunctionsStructure<N, Self, Self, Range, Range>
    where
        Self: ConstSizeFiniteSetSignature<N>,
    {
        ConstSizeFunctionsStructure::new(self, range)
    }

    fn const_size_functions_to<'a, const N: usize, Range: SetSignature>(
        &self,
        range: &'a Range,
    ) -> ConstSizeFunctionsStructure<N, Self, &Self, Range, &'a Range>
    where
        Self: ConstSizeFiniteSetSignature<N>,
    {
        ConstSizeFunctionsStructure::new(self, range)
    }
}
impl<Set: SetSignature> SetToConstSizeFunctionsToSignature for Set {}

pub trait SetToConstSizeFunctionsFromSignature: SetSignature {
    fn into_functions_from<const N: usize, Domain: ConstSizeFiniteSetSignature<N>>(
        self,
        domain: Domain,
    ) -> ConstSizeFunctionsStructure<N, Domain, Domain, Self, Self> {
        ConstSizeFunctionsStructure::new(domain, self)
    }

    fn functions_from<'a, const N: usize, Domain: ConstSizeFiniteSetSignature<N>>(
        &self,
        domain: &'a Domain,
    ) -> ConstSizeFunctionsStructure<N, Domain, &'a Domain, Self, &Self> {
        ConstSizeFunctionsStructure::new(domain, self)
    }
}
impl<Set: SetSignature> SetToConstSizeFunctionsFromSignature for Set {}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N>,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
> ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    pub fn domain(&self) -> &Domain {
        self.domain.borrow()
    }

    pub fn range(&self) -> &Range {
        self.range.borrow()
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N>,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
> Signature for ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
> FunctionsSignature<Domain, Range>
    for ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    fn function(
        &self,
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

    fn image<'a>(&self, f: &'a <Self as SetSignature>::Elem, x: &Domain::Elem) -> &'a Range::Elem {
        debug_assert!(self.is_element(f));
        debug_assert!(self.domain().is_element(x));
        &f[TryInto::<usize>::try_into(self.domain().element_to_enumeration(x)).unwrap()]
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
> SetSignature for ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    type Elem = Function<N, Domain::Elem, Range::Elem>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        for y in &x.images {
            self.range().validate_element(y)?;
        }
        Ok(())
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EqSignature,
    RangeB: BorrowedStructure<Range>,
> EqSignature for ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        (0..N).all(|i| self.range().equal(&a[i], &b[i]))
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: OrdSignature,
    RangeB: BorrowedStructure<Range>,
> PartialOrdSignature for ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: OrdSignature,
    RangeB: BorrowedStructure<Range>,
> OrdSignature for ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
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
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> CountableSetSignature for ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        let n: usize = self.domain().size().try_into().unwrap_or(usize::MAX);
        (0..n)
            .map(|_| self.range().list_all_elements())
            .multi_cartesian_product()
            .map(|f| f.try_into().unwrap())
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> FiniteSetSignature for ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    fn size(&self) -> Natural {
        self.range().size().pow(&self.domain().size())
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> EnumeratedOrdFiniteSetSignature
    for ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
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

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        debug_assert!(self.is_element(elem));
        let r = self.range().size();
        debug_assert_eq!(self.range().size(), (&r).into());
        let mut num = Natural::ZERO;
        for i in (0..N).rev() {
            num = num * &r + self.range().element_to_enumeration(&elem[i]);
        }
        num
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
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
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>>,
    DomainPerms: PermutationsSignature<Domain>,
    DomainPermsB: BorrowedStructure<DomainPerms>,
> {
    _functions: PhantomData<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>>,
    functions: FunctionsB,
    _domain_perms: PhantomData<DomainPerms>,
    domain_perms: DomainPermsB,
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>>,
    DomainPerms: PermutationsSignature<Domain>,
    DomainPermsB: BorrowedStructure<DomainPerms>,
>
    RightPermutationActionOnConstSizeFunctionsStructure<
        N,
        Domain,
        DomainB,
        Range,
        RangeB,
        FunctionsB,
        DomainPerms,
        DomainPermsB,
    >
{
    fn new(functions: FunctionsB, domain_perms: DomainPermsB) -> Self {
        Self {
            _functions: PhantomData,
            functions,
            _domain_perms: PhantomData,
            domain_perms,
        }
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>>,
    DomainPerms: PermutationsSignature<Domain>,
    DomainPermsB: BorrowedStructure<DomainPerms>,
> Signature
    for RightPermutationActionOnConstSizeFunctionsStructure<
        N,
        Domain,
        DomainB,
        Range,
        RangeB,
        FunctionsB,
        DomainPerms,
        DomainPermsB,
    >
{
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
> FunctionsDomainPermutationActionSignature<Domain, Range>
    for ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    type DomainPerms = FinitelySupportedPermutationsStructure<Domain, Domain>;
    type DomainPermsRef<'a>
        = FinitelySupportedPermutationsStructure<Domain, &'a Domain>
    where
        Self: 'a;

    fn into_domain_permutation_action(
        self,
    ) -> impl RightGroupActionSignature<Self, Self::DomainPerms> {
        let domain_perms = self.domain().clone().into_permutations();
        RightPermutationActionOnConstSizeFunctionsStructure::new(self, domain_perms)
    }

    fn domain_permutation_action<'a>(
        &'a self,
    ) -> impl RightGroupActionSignature<Self, Self::DomainPermsRef<'a>> {
        RightPermutationActionOnConstSizeFunctionsStructure::new(self, self.domain().permutations())
    }
}

/// Sym(D) has a right action on Fun(D -> R) by composition on the right
impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>>,
    DomainPerms: PermutationsSignature<Domain>,
    DomainPermsB: BorrowedStructure<DomainPerms>,
>
    RightGroupActionSignature<
        ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>,
        DomainPerms,
    >
    for RightPermutationActionOnConstSizeFunctionsStructure<
        N,
        Domain,
        DomainB,
        Range,
        RangeB,
        FunctionsB,
        DomainPerms,
        DomainPermsB,
    >
{
    fn group(&self) -> &DomainPerms {
        self.domain_perms.borrow()
    }

    fn set(&self) -> &ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB> {
        self.functions.borrow()
    }

    fn apply(
        &self,
        g: &<DomainPerms>::Elem,
        f: &<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB> as SetSignature>::Elem,
    ) -> <ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB> as SetSignature>::Elem
    {
        let functions = self.functions.borrow();
        let domain_perms = self.domain_perms.borrow();
        functions
            .function(|x| functions.image(f, &domain_perms.image(g, x)).clone())
            .unwrap()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LeftPermutationActionOnConstSizeFunctionsStructure<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N>,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>>,
    RangePerms: PermutationsSignature<Range>,
    RangePermsB: BorrowedStructure<RangePerms>,
> {
    _functions: PhantomData<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>>,
    functions: FunctionsB,
    _range_perms: PhantomData<RangePerms>,
    range_perms: RangePermsB,
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N>,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>>,
    RangePerms: PermutationsSignature<Range>,
    RangePermsB: BorrowedStructure<RangePerms>,
>
    LeftPermutationActionOnConstSizeFunctionsStructure<
        N,
        Domain,
        DomainB,
        Range,
        RangeB,
        FunctionsB,
        RangePerms,
        RangePermsB,
    >
{
    fn new(functions: FunctionsB, range_perms: RangePermsB) -> Self {
        Self {
            _functions: PhantomData,
            functions,
            _range_perms: PhantomData,
            range_perms,
        }
    }
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N>,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>>,
    RangePerms: PermutationsSignature<Range>,
    RangePermsB: BorrowedStructure<RangePerms>,
> Signature
    for LeftPermutationActionOnConstSizeFunctionsStructure<
        N,
        Domain,
        DomainB,
        Range,
        RangeB,
        FunctionsB,
        RangePerms,
        RangePermsB,
    >
{
}

impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> FunctionsRangePermutationActionSignature<Domain, Range>
    for ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>
{
    type RangePerms = FinitelySupportedPermutationsStructure<Range, Range>;
    type RangePermsRef<'a>
        = FinitelySupportedPermutationsStructure<Range, &'a Range>
    where
        Self: 'a;

    fn into_range_permutation_action(
        self,
    ) -> impl LeftGroupActionSignature<Self::RangePerms, Self> {
        let range_perms = self.range().clone().into_permutations();
        LeftPermutationActionOnConstSizeFunctionsStructure::new(self, range_perms)
    }

    fn range_permutation_action<'a>(
        &'a self,
    ) -> impl LeftGroupActionSignature<Self::RangePermsRef<'a>, Self> {
        LeftPermutationActionOnConstSizeFunctionsStructure::new(self, self.range().permutations())
    }
}

/// Sym(R) has a left action on Fun(D -> R) by composition on the left
impl<
    const N: usize,
    Domain: ConstSizeFiniteSetSignature<N> + EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>>,
    RangePerms: PermutationsSignature<Range>,
    RangePermsB: BorrowedStructure<RangePerms>,
>
    LeftGroupActionSignature<
        RangePerms,
        ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB>,
    >
    for LeftPermutationActionOnConstSizeFunctionsStructure<
        N,
        Domain,
        DomainB,
        Range,
        RangeB,
        FunctionsB,
        RangePerms,
        RangePermsB,
    >
{
    fn group(&self) -> &RangePerms {
        self.range_perms.borrow()
    }

    fn set(&self) -> &ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB> {
        self.functions.borrow()
    }

    fn apply(
        &self,
        g: &<RangePerms>::Elem,
        f: &<ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB> as SetSignature>::Elem,
    ) -> <ConstSizeFunctionsStructure<N, Domain, DomainB, Range, RangeB> as SetSignature>::Elem
    {
        let functions = self.functions.borrow();
        let range_perms = self.range_perms.borrow();
        functions
            .function(|x| range_perms.image(g, functions.image(f, x)))
            .unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sets::FiniteSetToFinitelySupportedPermutationsStructure;

    #[test]
    fn enumerate() {
        let set_a = i32::structure().into_const_size_finite_subset([1, 2, 3, 4, 5]);
        let set_b = i32::structure().into_const_size_finite_subset([1, 2, 3]);
        let fns = set_a.const_size_functions_to(&set_b);
        algebraeon_structures::assert_enumerated_ord_finite_set!(fns, 243);
    }

    #[test]
    fn test_image() {
        let set_a = i32::structure().into_const_size_finite_subset([1, 2, 3, 4, 5]);
        let set_b = i32::structure().into_const_size_finite_subset([1, 2, 3]);
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
        let set_a = i32::structure().into_const_size_finite_subset([1, 2, 3, 4, 5]);
        let set_a_perms = set_a.permutations();
        let set_b = i32::structure().into_const_size_finite_subset([1, 2, 3]);
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
                &fns.domain_permutation_action()
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
                &fns.range_permutation_action()
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
