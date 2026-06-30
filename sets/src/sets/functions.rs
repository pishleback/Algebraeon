use algebraeon_structures::*;
use itertools::Itertools;
use std::{cmp::Ordering, fmt::Debug, marker::PhantomData};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sets::SetToFiniteSubsetByOrdSignature;

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
}
