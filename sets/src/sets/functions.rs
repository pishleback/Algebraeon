use crate::sets::{
    FiniteSetToFinitelySupportedPermutationsStructure, FinitelySupportedPermutationsStructure,
};
use algebraeon_structures::*;
use itertools::Itertools;
use std::{cmp::Ordering, fmt::Debug, sync::Arc};

/// Represent all functions from `domain` to `range`
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FunctionsStructure<Domain: SetSignature, Range: SetSignature> {
    domain: Arc<Domain>,
    range: Arc<Range>,
}

impl<Domain: SetSignature, Range: SetSignature> FunctionsStructure<Domain, Range> {
    pub fn new(domain: Arc<Domain>, range: Arc<Range>) -> Arc<Self> {
        Self { domain, range }.into()
    }

    pub fn into_domain_and_range(self) -> (Arc<Domain>, Arc<Range>) {
        (self.domain, self.range)
    }
}

pub trait SetToFunctionsToSignature: SetSignature {
    fn functions_to<Range: SetSignature>(
        self: &Arc<Self>,
        range: &Arc<Range>,
    ) -> Arc<FunctionsStructure<Self, Range>> {
        FunctionsStructure::new(self.clone(), range.clone())
    }
}
impl<Set: SetSignature> SetToFunctionsToSignature for Set {}

pub trait SetToFunctionsFromSignature: SetSignature {
    fn functions_from<'a, Domain: SetSignature>(
        self: &Arc<Self>,
        domain: &Arc<Domain>,
    ) -> Arc<FunctionsStructure<Domain, Self>> {
        FunctionsStructure::new(domain.clone(), self.clone())
    }
}
impl<Set: SetSignature> SetToFunctionsFromSignature for Set {}

impl<Domain: SetSignature, Range: SetSignature> FunctionsStructure<Domain, Range> {
    pub fn domain(&self) -> &Domain {
        self.domain.borrow()
    }

    pub fn range(&self) -> &Range {
        self.range.borrow()
    }
}

impl<Domain: SetSignature, Range: SetSignature> Signature for FunctionsStructure<Domain, Range> {}

impl<Domain: OrderedFiniteSetSignature, Range: SetSignature> FunctionsSignature<Domain, Range>
    for FunctionsStructure<Domain, Range>
{
    fn function(
        self: &Arc<Self>,
        f: impl Fn(&Domain::Elem) -> Range::Elem,
    ) -> Option<Vec<Range::Elem>> {
        let s = self
            .domain()
            .list_all_elements_ordered()
            .iter()
            .map(f)
            .collect();
        for y in &s {
            if !self.range().is_element(y) {
                return None;
            }
        }
        Some(s)
    }

    fn image<'a>(self: &Arc<Self>, f: &'a Vec<Range::Elem>, x: &Domain::Elem) -> &'a Range::Elem {
        debug_assert!(self.is_element(f));
        debug_assert!(self.domain().is_element(x));
        &f[TryInto::<usize>::try_into(self.domain().element_to_enumeration(x)).unwrap()]
    }
}

impl<Domain: OrderedFiniteSetSignature, Range: SetSignature> SetSignature
    for FunctionsStructure<Domain, Range>
{
    type Elem = Vec<Range::Elem>;

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        if Natural::from(x.len()) != self.domain().size() {
            return Err("Incorrect vector length".to_string());
        }
        for y in x {
            self.range().validate_element(y)?;
        }
        Ok(())
    }
}

impl<Domain: OrderedFiniteSetSignature, Range: EqSignature> EqSignature
    for FunctionsStructure<Domain, Range>
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        let n = a.len();
        debug_assert_eq!(n, b.len());
        (0..n).all(|i| self.range().equal(&a[i], &b[i]))
    }
}

impl<Domain: OrderedFiniteSetSignature, Range: OrdSignature> PartialOrdSignature
    for FunctionsStructure<Domain, Range>
{
    fn partial_cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<std::cmp::Ordering> {
        Some(self.cmp(a, b))
    }
}

impl<Domain: OrderedFiniteSetSignature, Range: OrdSignature> OrdSignature
    for FunctionsStructure<Domain, Range>
{
    fn cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        debug_assert!(self.is_element(a));
        debug_assert!(self.is_element(b));
        let n = a.len();
        debug_assert_eq!(n, b.len());
        for i in (0..n).rev() {
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

impl<Domain: OrderedFiniteSetSignature, Range: OrderedFiniteSetSignature> CountableSetSignature
    for FunctionsStructure<Domain, Range>
{
    fn generate_all_elements(self: Arc<Self>) -> impl Iterator<Item = Self::Elem> {
        let n: usize = self.domain().size().try_into().unwrap_or(usize::MAX);
        (0..n)
            .map(|_| self.range().list_all_elements())
            .multi_cartesian_product()
    }
}

impl<Domain: OrderedFiniteSetSignature, Range: OrderedFiniteSetSignature> FiniteSetSignature
    for FunctionsStructure<Domain, Range>
{
    fn size(self: &Arc<Self>) -> Natural {
        self.range().size().pow(&self.domain().size())
    }
}

impl<Domain: OrderedFiniteSetSignature, Range: OrderedFiniteSetSignature> OrderedFiniteSetSignature
    for FunctionsStructure<Domain, Range>
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
            .map(|images| images.into_iter().rev().collect())
            .collect()
    }

    fn element_to_enumeration(self: &Arc<Self>, elem: &Self::Elem) -> Natural {
        debug_assert!(self.is_element(elem));
        let d = elem.len();
        let r = self.range().size();
        debug_assert_eq!(self.range().size(), (&r).into());
        let mut num = Natural::ZERO;
        for i in (0..d).rev() {
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
        debug_assert!(self.is_element(&f));
        Some(f)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RightPermutationActionOnFunctionsStructure<
    Domain: OrderedFiniteSetSignature,
    Range: SetSignature,
    DomainPerms: PermutationsSignature<Domain>,
> {
    functions: Arc<FunctionsStructure<Domain, Range>>,
    domain_perms: Arc<DomainPerms>,
}

impl<
    Domain: OrderedFiniteSetSignature,
    Range: SetSignature,
    DomainPerms: PermutationsSignature<Domain>,
> RightPermutationActionOnFunctionsStructure<Domain, Range, DomainPerms>
{
    fn new(
        functions: Arc<FunctionsStructure<Domain, Range>>,
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
    Domain: OrderedFiniteSetSignature,
    Range: SetSignature,
    DomainPerms: PermutationsSignature<Domain>,
> Signature for RightPermutationActionOnFunctionsStructure<Domain, Range, DomainPerms>
{
}

impl<Domain: OrderedFiniteSetSignature, Range: SetSignature> FunctionsStructure<Domain, Range> {
    pub fn domain_finitely_supported_permutation_action(
        self: &Arc<Self>,
    ) -> impl RightGroupActionSignature<Self, FinitelySupportedPermutationsStructure<Domain, &Domain>>
    {
        RightPermutationActionOnFunctionsStructure::new(self.clone(), self.domain().permutations())
    }
}

/// Sym(D) has a right action on Fun(D -> R) by composition on the right
impl<
    Domain: OrderedFiniteSetSignature,
    Range: SetSignature,
    DomainPerms: PermutationsSignature<Domain>,
> RightGroupActionSignature<FunctionsStructure<Domain, Range>, DomainPerms>
    for RightPermutationActionOnFunctionsStructure<Domain, Range, DomainPerms>
{
    fn group(self: &Arc<Self>) -> &Arc<DomainPerms> {
        &self.domain_perms
    }

    fn set(self: &Arc<Self>) -> &Arc<FunctionsStructure<Domain, Range>> {
        &self.functions
    }

    fn apply(
        self: &Arc<Self>,
        g: &<DomainPerms>::Elem,
        f: &<FunctionsStructure<Domain, Range> as SetSignature>::Elem,
    ) -> <FunctionsStructure<Domain, Range> as SetSignature>::Elem {
        let functions = self.functions.borrow();
        let domain_perms = self.domain_perms.borrow();
        functions
            .function(|x| functions.image(f, &domain_perms.image(g, x)).clone())
            .unwrap()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LeftPermutationActionOnFunctionsStructure<
    Domain: SetSignature,
    Range: OrderedFiniteSetSignature,
    RangePerms: PermutationsSignature<Range>,
> {
    functions: Arc<FunctionsStructure<Domain, Range>>,
    range_perms: Arc<RangePerms>,
}

impl<
    Domain: SetSignature,
    Range: OrderedFiniteSetSignature,
    RangePerms: PermutationsSignature<Range>,
> LeftPermutationActionOnFunctionsStructure<Domain, Range, RangePerms>
{
    fn new(
        functions: Arc<FunctionsStructure<Domain, Range>>,
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
    Domain: SetSignature,
    Range: OrderedFiniteSetSignature,
    RangePerms: PermutationsSignature<Range>,
> Signature for LeftPermutationActionOnFunctionsStructure<Domain, Range, RangePerms>
{
}

impl<Domain: OrderedFiniteSetSignature, Range: OrderedFiniteSetSignature>
    FunctionsStructure<Domain, Range>
{
    pub fn range_finitely_supported_permutation_action(
        self: &Arc<Self>,
    ) -> impl LeftGroupActionSignature<FinitelySupportedPermutationsStructure<Range, &Range>, Self>
    {
        LeftPermutationActionOnFunctionsStructure::new(self.clone(), self.range().permutations())
    }
}

/// Sym(R) has a left action on Fun(D -> R) by composition on the left
impl<
    Domain: OrderedFiniteSetSignature,
    Range: OrderedFiniteSetSignature,
    RangePerms: PermutationsSignature<Range>,
> LeftGroupActionSignature<RangePerms, FunctionsStructure<Domain, Range>>
    for LeftPermutationActionOnFunctionsStructure<Domain, Range, RangePerms>
{
    fn group(self: &Arc<Self>) -> &Arc<RangePerms> {
        &self.range_perms
    }

    fn set(self: &Arc<Self>) -> &Arc<FunctionsStructure<Domain, Range>> {
        &self.functions
    }

    fn apply(
        self: &Arc<Self>,
        g: &<RangePerms>::Elem,
        f: &<FunctionsStructure<Domain, Range> as SetSignature>::Elem,
    ) -> <FunctionsStructure<Domain, Range> as SetSignature>::Elem {
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
        let set_a = i32::structure().into_finite_subset(vec![1, 2, 3, 4, 5]);
        let set_b = i32::structure().into_finite_subset(vec![1, 2, 3]);
        let fns = set_a.functions_to(&set_b);
        algebraeon_structures::assert_enumerated_ord_finite_set!(fns, 243);
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
                &fns.domain_finitely_supported_permutation_action()
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
