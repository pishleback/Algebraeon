use algebraeon_structures::*;
use itertools::Itertools;
use std::{cmp::Ordering, fmt::Debug, marker::PhantomData};

use crate::sets::{
    FiniteSetToFinitelySupportedPermutationsStructure, FinitelySupportedPermutationsStructure,
};

#[derive(Debug, Clone)]
pub struct Function<ElemTo> {
    pub(crate) images: Vec<ElemTo>,
}

impl<ElemTo> Function<ElemTo> {
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.images.len()
    }
}

/// Represent all functions from `domain` to `range`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FunctionsStructure<
    Domain: SetSignature,
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
    Domain: SetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
> FunctionsStructure<Domain, DomainB, Range, RangeB>
{
    pub fn new(domain: DomainB, range: RangeB) -> Self {
        Self {
            _domain: PhantomData,
            domain,
            _range: PhantomData,
            range,
        }
    }
}

pub trait SetToFunctionsToSignature: SetSignature {
    fn into_functions_to<Range: SetSignature>(
        self,
        range: Range,
    ) -> FunctionsStructure<Self, Self, Range, Range> {
        FunctionsStructure::new(self, range)
    }

    fn functions_to<'a, Range: SetSignature>(
        &self,
        range: &'a Range,
    ) -> FunctionsStructure<Self, &Self, Range, &'a Range> {
        FunctionsStructure::new(self, range)
    }
}
impl<Set: SetSignature> SetToFunctionsToSignature for Set {}

pub trait SetToFunctionsFromSignature: SetSignature {
    fn into_functions_from<Domain: SetSignature>(
        self,
        domain: Domain,
    ) -> FunctionsStructure<Domain, Domain, Self, Self> {
        FunctionsStructure::new(domain, self)
    }

    fn functions_from<'a, Domain: SetSignature>(
        &self,
        domain: &'a Domain,
    ) -> FunctionsStructure<Domain, &'a Domain, Self, &Self> {
        FunctionsStructure::new(domain, self)
    }
}
impl<Set: SetSignature> SetToFunctionsFromSignature for Set {}

impl<
    Domain: SetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
> FunctionsStructure<Domain, DomainB, Range, RangeB>
{
    pub fn domain(&self) -> &Domain {
        self.domain.borrow()
    }

    pub fn range(&self) -> &Range {
        self.range.borrow()
    }
}

impl<
    Domain: SetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: SetSignature,
    RangeB: BorrowedStructure<Range>,
> Signature for FunctionsStructure<Domain, DomainB, Range, RangeB>
{
}

impl<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> FunctionsStructure<Domain, DomainB, Range, RangeB>
{
    pub fn function(
        &self,
        f: impl Fn(&Domain::Elem) -> Range::Elem,
    ) -> Option<Function<Range::Elem>> {
        let s = Function {
            images: self
                .domain()
                .list_all_elements_ordered()
                .iter()
                .map(f)
                .collect(),
        };
        for y in &s.images {
            if !self.range().is_element(y) {
                return None;
            }
        }
        Some(s)
    }

    pub fn image<'a>(&self, f: &'a Function<Range::Elem>, x: &Domain::Elem) -> &'a Range::Elem {
        debug_assert!(self.is_element(f));
        debug_assert!(self.domain().is_element(x));
        &f.images[TryInto::<usize>::try_into(self.domain().element_to_enumeration(x)).unwrap()]
    }
}

impl<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> SetSignature for FunctionsStructure<Domain, DomainB, Range, RangeB>
{
    type Elem = Function<Range::Elem>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if Natural::from(x.len()) != self.domain().size() {
            return Err("Incorrect vector length".to_string());
        }
        for y in &x.images {
            self.range().validate_element(y)?;
        }
        Ok(())
    }
}

impl<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> EqSignature for FunctionsStructure<Domain, DomainB, Range, RangeB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        let n = a.len();
        debug_assert_eq!(n, b.len());
        (0..n).all(|i| self.range().equal(&a.images[i], &b.images[i]))
    }
}

impl<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> PartialOrdSignature for FunctionsStructure<Domain, DomainB, Range, RangeB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> OrdSignature for FunctionsStructure<Domain, DomainB, Range, RangeB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        let n = a.len();
        debug_assert_eq!(n, b.len());
        for i in (0..n).rev() {
            match self.range().cmp(&a.images[i], &b.images[i]) {
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
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> CountableSetSignature for FunctionsStructure<Domain, DomainB, Range, RangeB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        let n: usize = self.domain().size().try_into().unwrap_or(usize::MAX);
        (0..n)
            .map(|_| self.range().list_all_elements())
            .multi_cartesian_product()
            .map(|images| Function { images })
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.clone().into_generate_all_elements()
    }
}

impl<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> FiniteSetSignature for FunctionsStructure<Domain, DomainB, Range, RangeB>
{
    fn size(&self) -> Natural {
        self.range().size().pow(&self.domain().size())
    }
}

impl<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> EnumeratedOrdFiniteSetSignature for FunctionsStructure<Domain, DomainB, Range, RangeB>
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
            .map(|images| Function {
                images: images.into_iter().rev().collect(),
            })
            .collect()
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        debug_assert!(self.is_element(elem));
        let d = elem.images.len();
        let r = self.range().size();
        debug_assert_eq!(self.range().size(), (&r).into());
        let mut num = Natural::ZERO;
        for i in (0..d).rev() {
            num = num * &r + self.range().element_to_enumeration(&elem.images[i]);
        }
        num
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        if *num >= self.size() {
            return None;
        }
        let d: usize = self.domain().size().try_into().unwrap();
        let mut images = vec![];
        let r = self.range().size();
        let mut n = num.clone();
        let mut k;
        for _ in 0..d {
            (n, k) = n.div_mod(&r);
            images.push(self.range().enumeration_to_element(&k).unwrap())
        }
        let f = Function { images };
        debug_assert!(self.is_element(&f));
        Some(f)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RightPermutationActionOnFunctionsStructure<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<FunctionsStructure<Domain, DomainB, Range, RangeB>>,
    DomainPerms: PermutationsSignature<Domain>,
    DomainPermsB: BorrowedStructure<DomainPerms>,
> {
    _functions: PhantomData<FunctionsStructure<Domain, DomainB, Range, RangeB>>,
    functions: FunctionsB,
    _domain_perms: PhantomData<DomainPerms>,
    domain_perms: DomainPermsB,
}

impl<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<FunctionsStructure<Domain, DomainB, Range, RangeB>>,
    DomainPerms: PermutationsSignature<Domain>,
    DomainPermsB: BorrowedStructure<DomainPerms>,
>
    RightPermutationActionOnFunctionsStructure<
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
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<FunctionsStructure<Domain, DomainB, Range, RangeB>>,
    DomainPerms: PermutationsSignature<Domain>,
    DomainPermsB: BorrowedStructure<DomainPerms>,
> Signature
    for RightPermutationActionOnFunctionsStructure<
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
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> FunctionsStructure<Domain, DomainB, Range, RangeB>
{
    pub fn domain_permutation_action(
        &self,
    ) -> impl RightGroupActionSignature<Self, FinitelySupportedPermutationsStructure<Domain, &Domain>>
    {
        RightPermutationActionOnFunctionsStructure::new(self, self.domain().permutations())
    }
}

/// Sym(D) has a right action on Fun(D -> R) by composition on the right
impl<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<FunctionsStructure<Domain, DomainB, Range, RangeB>>,
    DomainPerms: PermutationsSignature<Domain>,
    DomainPermsB: BorrowedStructure<DomainPerms>,
> RightGroupActionSignature<FunctionsStructure<Domain, DomainB, Range, RangeB>, DomainPerms>
    for RightPermutationActionOnFunctionsStructure<
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

    fn set(&self) -> &FunctionsStructure<Domain, DomainB, Range, RangeB> {
        self.functions.borrow()
    }

    fn apply(
        &self,
        g: &<DomainPerms>::Elem,
        f: &<FunctionsStructure<Domain, DomainB, Range, RangeB> as SetSignature>::Elem,
    ) -> <FunctionsStructure<Domain, DomainB, Range, RangeB> as SetSignature>::Elem {
        let functions = self.functions.borrow();
        let domain_perms = self.domain_perms.borrow();
        functions
            .function(|x| functions.image(f, &domain_perms.image(g, x)).clone())
            .unwrap()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LeftPermutationActionOnFunctionsStructure<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<FunctionsStructure<Domain, DomainB, Range, RangeB>>,
    RangePerms: PermutationsSignature<Range>,
    RangePermsB: BorrowedStructure<RangePerms>,
> {
    _functions: PhantomData<FunctionsStructure<Domain, DomainB, Range, RangeB>>,
    functions: FunctionsB,
    _range_perms: PhantomData<RangePerms>,
    range_perms: RangePermsB,
}

impl<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<FunctionsStructure<Domain, DomainB, Range, RangeB>>,
    RangePerms: PermutationsSignature<Range>,
    RangePermsB: BorrowedStructure<RangePerms>,
>
    LeftPermutationActionOnFunctionsStructure<
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
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<FunctionsStructure<Domain, DomainB, Range, RangeB>>,
    RangePerms: PermutationsSignature<Range>,
    RangePermsB: BorrowedStructure<RangePerms>,
> Signature
    for LeftPermutationActionOnFunctionsStructure<
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
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
> FunctionsStructure<Domain, DomainB, Range, RangeB>
{
    pub fn range_permutation_action(
        &self,
    ) -> impl LeftGroupActionSignature<FinitelySupportedPermutationsStructure<Range, &Range>, Self>
    {
        LeftPermutationActionOnFunctionsStructure::new(self, self.range().permutations())
    }
}

/// Sym(R) has a left action on Fun(D -> R) by composition on the left
impl<
    Domain: EnumeratedOrdFiniteSetSignature,
    DomainB: BorrowedStructure<Domain>,
    Range: EnumeratedOrdFiniteSetSignature,
    RangeB: BorrowedStructure<Range>,
    FunctionsB: BorrowedStructure<FunctionsStructure<Domain, DomainB, Range, RangeB>>,
    RangePerms: PermutationsSignature<Range>,
    RangePermsB: BorrowedStructure<RangePerms>,
> LeftGroupActionSignature<RangePerms, FunctionsStructure<Domain, DomainB, Range, RangeB>>
    for LeftPermutationActionOnFunctionsStructure<
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

    fn set(&self) -> &FunctionsStructure<Domain, DomainB, Range, RangeB> {
        self.functions.borrow()
    }

    fn apply(
        &self,
        g: &<RangePerms>::Elem,
        f: &<FunctionsStructure<Domain, DomainB, Range, RangeB> as SetSignature>::Elem,
    ) -> <FunctionsStructure<Domain, DomainB, Range, RangeB> as SetSignature>::Elem {
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
    use crate::sets::{
        FiniteSetToFinitelySupportedPermutationsStructure, SetToFiniteSubsetByOrdSignature,
    };

    #[test]
    fn test_enumeration() {
        let set_a = i32::structure().into_finite_subset(vec![1, 2, 3, 4, 5]);
        let set_b = i32::structure().into_finite_subset(vec![1, 2, 3]);
        let fns = set_a.functions_to(&set_b);

        assert_eq!(fns.size(), fns.list_all_elements().len().into());
        assert_eq!(fns.size(), fns.list_all_elements_ordered().len().into());

        let all_fns = fns.list_all_elements_ordered();
        for i in 0..(all_fns.len() - 1) {
            assert!(fns.cmp(&all_fns[i], &all_fns[i + 1]).is_lt());
        }
        for (i, s) in all_fns.iter().enumerate() {
            assert_eq!(Natural::from(i), fns.element_to_enumeration(s));
            assert!(fns.equal(&fns.enumeration_to_element(&Natural::from(i)).unwrap(), s));
        }
        assert!(fns.enumeration_to_element(&fns.size()).is_none());
    }

    #[test]
    fn test_image() {
        let set_a = i32::structure().into_finite_subset(vec![1, 2, 3, 4, 5]);
        let set_b = i32::structure().into_finite_subset(vec![1, 2, 3]);
        let fns = set_a.functions_to(&set_b);
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
        let set_a = i32::structure().into_finite_subset(vec![1, 2, 3, 4, 5]);
        let set_a_perms = set_a.permutations();
        let set_b = i32::structure().into_finite_subset(vec![1, 2, 3]);
        let set_b_perms = set_b.permutations();
        let fns = set_a.functions_to(&set_b);

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
