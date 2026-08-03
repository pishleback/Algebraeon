use crate::structure::*;
use algebraeon_structures::*;
use std::{borrow::Cow, rc::Rc};

#[derive(Debug, Clone)]
pub struct FreeModuleOverOrderedSetStructure<Set: OrdSignature, Ring: SemiRingSignature> {
    set: Rc<Set>,
    ring: Rc<Ring>,
}

impl<Set: OrdSignature, Ring: SemiRingSignature> PartialEq
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
    fn eq(&self, other: &Self) -> bool {
        self.set == other.set && self.ring == other.ring
    }
}

impl<Set: OrdSignature, Ring: SemiRingSignature> Eq
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
}

impl<Set: OrdSignature, Ring: SemiRingSignature + EqSignature>
    FreeModuleOverOrderedSetStructure<Set, Ring>
{
    pub fn new(set: Rc<Set>, ring: Rc<Ring>) -> Rc<Self> {
        Self { set, ring }.into()
    }

    pub fn set(&self) -> &Rc<Set> {
        &self.set
    }

    /// Input: vector of (Set, Ring)
    /// Output: vector of (Set, Ring) which is
    ///  - ordered
    ///  - has no duplicate set elements
    ///  - has no ring elements equal to 0
    ///
    /// and is equal to the sum of the input vectors terms. In other words, the returned vector will pass self.is_element(..)
    pub fn collapse_terms(
        self: &Rc<Self>,
        v: <Self as SetSignature>::Elem,
    ) -> <Self as SetSignature>::Elem {
        let mut v = self.set().sort_by_key(v, &|(x, _)| x).into_iter();

        let mut current_x = None;
        let mut current_a = self.ring().zero();

        let mut w = vec![];
        loop {
            enum ItemResult<S, R> {
                Same(R),
                Change(S, R),
                End,
            }

            let result = if let Some((x, a)) = v.next() {
                if let Some(current_x) = current_x.as_ref() {
                    if self.set().equal(current_x, &x) {
                        ItemResult::Same(a)
                    } else {
                        ItemResult::Change(x, a)
                    }
                } else {
                    current_x = Some(x);
                    ItemResult::Same(a)
                }
            } else {
                ItemResult::End
            };

            match result {
                ItemResult::Same(a) => {
                    self.ring().add_mut(&mut current_a, &a);
                }
                ItemResult::Change(x, a) => {
                    w.push((current_x.unwrap(), current_a));
                    current_x = Some(x);
                    current_a = a;
                }
                ItemResult::End => {
                    w.push((current_x.unwrap(), current_a));
                    break;
                }
            }
        }
        let w = w
            .into_iter()
            .filter(|(_, a)| !self.ring().is_zero(a))
            .collect();
        debug_assert!(self.validate_element(&w).is_ok());
        w
    }
}

impl<Set: OrdSignature, Ring: SemiRingSignature> Signature
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
}

impl<Set: OrdSignature, Ring: SemiRingSignature + EqSignature> SetSignature
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
    // must be ordered and contain no duplicates wrt the first argument
    // all ring elements in the second argument must be non-zero
    type Elem = Vec<(Set::Elem, Ring::Elem)>;

    fn validate_element(self: &Rc<Self>, v: &Self::Elem) -> Result<(), String> {
        if !self.set().is_sorted_and_unique_by_key(v, |(x, _)| x) {
            return Err("not sorted or has duplicate".to_string());
        }
        for (_, a) in v {
            if self.ring().is_zero(a) {
                return Err("multiplicity zero".to_string());
            }
        }
        Ok(())
    }
}

impl<Set: OrdSignature, Ring: SemiRingSignature + EqSignature> EqSignature
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
    fn equal(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> bool {
        debug_assert!(self.validate_element(v).is_ok());
        debug_assert!(self.validate_element(w).is_ok());
        // since elements are sorted and exclude entries with zero coefficients, we just need to check if they are identically equal
        let n = v.len();
        if n != w.len() {
            false
        } else {
            (0..n).all(|i| {
                let (vx, va) = &v[i];
                let (wx, wa) = &w[i];
                self.set().equal(vx, wx) && self.ring().equal(va, wa)
            })
        }
    }
}

impl<Set: OrdSignature, Ring: SemiRingSignature + EqSignature> RinglikeSpecializationSignature
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
}

impl<Set: OrdSignature, Ring: SemiRingSignature + EqSignature> ZeroSignature
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
    fn zero(self: &Rc<Self>) -> Self::Elem {
        vec![]
    }
}

impl<Set: OrdSignature, Ring: SemiRingSignature + EqSignature> AdditionSignature
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
    fn add(self: &Rc<Self>, v: &Self::Elem, w: &Self::Elem) -> Self::Elem {
        self.collapse_terms(v.iter().chain(w.iter()).cloned().collect())
    }
}

impl<Set: OrdSignature, Ring: SemiRingSignature + EqSignature> CancellativeAdditionSignature
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
    fn try_sub(self: &Rc<Self>, _a: &Self::Elem, _b: &Self::Elem) -> Option<Self::Elem> {
        todo!()
    }
}

impl<Set: OrdSignature, Ring: SemiRingSignature + EqSignature> TryNegateSignature
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
    fn try_neg(self: &Rc<Self>, v: &Self::Elem) -> Option<Self::Elem> {
        v.iter()
            .map(|(x, a)| Some((x.clone(), self.ring().try_neg(a)?)))
            .collect::<Option<_>>()
    }
}

impl<Set: OrdSignature, Ring: SemiRingSignature + EqSignature> AdditiveMonoidSignature
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
}

impl<Set: OrdSignature, Ring: RingSignature + EqSignature> AdditiveGroupSignature
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
    fn neg(self: &Rc<Self>, v: &Self::Elem) -> Self::Elem {
        v.iter()
            .map(|(x, a)| (x.clone(), self.ring().neg(a)))
            .collect()
    }
}

impl<Set: OrdSignature, Ring: SemiRingSignature + EqSignature> SemiModuleSignature<Ring>
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
    fn ring(self: &Rc<Self>) -> Rc<Ring> {
        self.ring.clone()
    }

    fn scalar_mul(self: &Rc<Self>, v: &Self::Elem, b: &Ring::Elem) -> Self::Elem {
        v.iter()
            .map(|(x, a)| (x.clone(), self.ring().mul(a, b)))
            .filter(|(_, a)| !self.ring().is_zero(a))
            .collect()
    }
}

impl<Set: OrdSignature, Ring: RingSignature + EqSignature> FreeModuleSignature<Set, Ring>
    for FreeModuleOverOrderedSetStructure<Set, Ring>
{
    fn basis_set(self: &Rc<Self>) -> Rc<Set> {
        self.set().clone()
    }

    fn to_component<'a>(self: &Rc<Self>, x: &Set::Elem, v: &'a Self::Elem) -> Cow<'a, Ring::Elem> {
        if let Some((_, a)) = self.set().binary_search_by_key(v, x, |(x, _)| x) {
            Cow::Borrowed(a)
        } else {
            Cow::Owned(self.ring().zero())
        }
    }

    fn from_component(self: &Rc<Self>, x: &Set::Elem, a: &Ring::Elem) -> Self::Elem {
        if self.ring().is_zero(a) {
            vec![]
        } else {
            vec![(x.clone(), a.clone())]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::FreeModuleOverOrderedSetStructure;
    use algebraeon_structures::*;

    #[test]
    fn test_ordered_set_free_module() {
        let m = FreeModuleOverOrderedSetStructure::new(Natural::structure(), Integer::structure());

        let v = vec![
            (0u32.into(), 1.into()),
            (1u32.into(), 1.into()),
            (0u32.into(), (-1).into()),
            (1u32.into(), 1.into()),
        ];
        let w = m.collapse_terms(v);
        assert_eq!(w, vec![(1u32.into(), 2.into())]);
    }
}
