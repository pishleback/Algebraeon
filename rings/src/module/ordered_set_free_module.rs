use crate::structure::*;
use algebraeon_sets::structure::*;
use std::{borrow::Cow, marker::PhantomData};

#[derive(Debug, Clone)]
pub struct FreeModuleOverOrderedSetStructure<
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> {
    _set: PhantomData<Set>,
    set: SetB,
    _ring: PhantomData<Ring>,
    ring: RingB,
}

impl<
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> PartialEq for FreeModuleOverOrderedSetStructure<Set, SetB, Ring, RingB>
{
    fn eq(&self, other: &Self) -> bool {
        self.set == other.set && self.ring == other.ring
    }
}

impl<
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> Eq for FreeModuleOverOrderedSetStructure<Set, SetB, Ring, RingB>
{
}

impl<
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> FreeModuleOverOrderedSetStructure<Set, SetB, Ring, RingB>
{
    pub fn new(set: SetB, ring: RingB) -> Self {
        Self {
            _set: PhantomData,
            set,
            _ring: PhantomData,
            ring,
        }
    }

    pub fn set(&self) -> &Set {
        self.set.borrow()
    }

    /// Input: vector of (Set, Ring)
    /// Output: vector of (Set, Ring) which is
    ///  - ordered
    ///  - has no duplicate set elements
    ///  - has no ring elements equal to 0
    ///
    /// and is equal to the sum of the input vectors terms. In other words, the returned vector will pass self.is_element(..)
    pub fn collapse_terms(&self, v: <Self as SetSignature>::Set) -> <Self as SetSignature>::Set {
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
        debug_assert!(self.is_element(&w).is_ok());
        w
    }
}

impl<
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> Signature for FreeModuleOverOrderedSetStructure<Set, SetB, Ring, RingB>
{
}

impl<
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> SetSignature for FreeModuleOverOrderedSetStructure<Set, SetB, Ring, RingB>
{
    // must be ordered and contain no duplicates wrt the first argument
    // all ring elements in the second argument must be non-zero
    type Set = Vec<(Set::Set, Ring::Set)>;

    fn is_element(&self, v: &Self::Set) -> Result<(), String> {
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

impl<
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> EqSignature for FreeModuleOverOrderedSetStructure<Set, SetB, Ring, RingB>
{
    fn equal(&self, v: &Self::Set, w: &Self::Set) -> bool {
        debug_assert!(self.is_element(v).is_ok());
        debug_assert!(self.is_element(w).is_ok());
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

impl<
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditiveMonoidSignature for FreeModuleOverOrderedSetStructure<Set, SetB, Ring, RingB>
{
    fn zero(&self) -> Self::Set {
        vec![]
    }

    fn add(&self, v: &Self::Set, w: &Self::Set) -> Self::Set {
        self.collapse_terms(v.iter().chain(w.iter()).cloned().collect())
    }
}

impl<
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditiveGroupSignature for FreeModuleOverOrderedSetStructure<Set, SetB, Ring, RingB>
{
    fn neg(&self, v: &Self::Set) -> Self::Set {
        v.iter()
            .map(|(x, a)| (x.clone(), self.ring().neg(a)))
            .collect()
    }
}

impl<
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> SemiModuleSignature<Ring> for FreeModuleOverOrderedSetStructure<Set, SetB, Ring, RingB>
{
    fn ring(&self) -> Cow<Ring> {
        Cow::Borrowed(self.ring.borrow())
    }

    fn scalar_mul(&self, v: &Self::Set, b: &Ring::Set) -> Self::Set {
        v.iter()
            .map(|(x, a)| (x.clone(), self.ring().mul(a, b)))
            .filter(|(_, a)| !self.ring().is_zero(a))
            .collect()
    }
}

impl<
    Set: OrdSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> FreeModuleSignature<Ring> for FreeModuleOverOrderedSetStructure<Set, SetB, Ring, RingB>
{
    type Basis = Set;

    fn basis_set(&self) -> impl std::borrow::Borrow<Self::Basis> {
        self.set()
    }

    fn to_component<'a>(&self, x: &Set::Set, v: &'a Self::Set) -> Cow<'a, Ring::Set> {
        if let Some((_, a)) = self.set().binary_search_by_key(v, x, |(x, _)| x) {
            Cow::Borrowed(a)
        } else {
            Cow::Owned(self.ring().zero())
        }
    }

    fn from_component(&self, x: &Set::Set, a: &Ring::Set) -> Self::Set {
        if self.ring().is_zero(a) {
            vec![]
        } else {
            vec![(x.clone(), a.clone())]
        }
    }
}

impl<
    Set: OrdSignature + FiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> FinitelyFreeModuleSignature<Ring> for FreeModuleOverOrderedSetStructure<Set, SetB, Ring, RingB>
{
}

#[cfg(test)]
mod tests {
    use algebraeon_nzq::{Integer, Natural};
    use algebraeon_sets::structure::MetaType;

    use super::FreeModuleOverOrderedSetStructure;

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
