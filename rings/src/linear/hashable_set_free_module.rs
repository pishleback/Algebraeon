use crate::structure::*;
use algebraeon_sets::structure::*;
use itertools::Itertools;
use std::fmt::Debug;
use std::{collections::HashMap, hash::Hash, marker::PhantomData};

pub trait FreeModuleOverHashableSetElement: Debug + Clone + Eq + Hash {}
impl<Set: Debug + Clone + Eq + Hash> FreeModuleOverHashableSetElement for Set {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FreeModuleOverHashableSetStructure<
    Set: FreeModuleOverHashableSetElement,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> {
    _set: PhantomData<Set>,
    _ring: PhantomData<Ring>,
    ring: RingB,
}

impl<Set: FreeModuleOverHashableSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    FreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    pub fn new(ring: RingB) -> Self {
        Self {
            _set: PhantomData::default(),
            _ring: PhantomData::default(),
            ring,
        }
    }

    pub fn reduce(&self, v: <Self as SetSignature>::Set) -> <Self as SetSignature>::Set {
        v.into_iter()
            .filter(|(_, r)| !self.ring().is_zero(r))
            .collect()
    }

    pub fn basis_element(&self, b: Set) -> <Self as SetSignature>::Set {
        [(b, self.ring().one())].into()
    }
}

impl<Set: FreeModuleOverHashableSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    Signature for FreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
}

impl<Set: FreeModuleOverHashableSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    SetSignature for FreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    type Set = HashMap<Set, Ring::Set>;

    fn is_element(&self, v: &Self::Set) -> bool {
        v.iter().all(|(_, r)| self.ring().is_element(r))
    }
}

impl<Set: FreeModuleOverHashableSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    EqSignature for FreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    fn equal(&self, v: &Self::Set, w: &Self::Set) -> bool {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        let mut keys = v.keys().chain(w.keys()).unique();
        keys.all(|k| match (v.get(k), w.get(k)) {
            (None, None) => unreachable!(),
            (None, Some(s)) => self.ring().is_zero(s),
            (Some(r), None) => self.ring().is_zero(r),
            (Some(r), Some(s)) => self.ring().equal(r, s),
        })
    }
}

impl<Set: FreeModuleOverHashableSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    AdditiveMonoidSignature for FreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    fn zero(&self) -> Self::Set {
        [].into()
    }

    fn add(&self, v: &Self::Set, w: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        let keys = v.keys().chain(w.keys()).unique();
        self.reduce(
            keys.map(|k| {
                (
                    k.clone(),
                    match (v.get(k), w.get(k)) {
                        (None, None) => unreachable!(),
                        (None, Some(s)) => s.clone(),
                        (Some(r), None) => r.clone(),
                        (Some(r), Some(s)) => self.ring().add(r, s),
                    },
                )
            })
            .collect(),
        )
    }
}

impl<Set: FreeModuleOverHashableSetElement, Ring: RingSignature, RingB: BorrowedStructure<Ring>>
    AdditiveGroupSignature for FreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    fn neg(&self, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        v.iter()
            .map(|(b, r)| (b.clone(), self.ring().neg(r)))
            .collect()
    }
}

impl<Set: FreeModuleOverHashableSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    SemiModuleSignature<Ring> for FreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    fn ring(&self) -> &Ring {
        self.ring.borrow()
    }

    fn scalar_mul(&self, r: &<Ring>::Set, v: &Self::Set) -> Self::Set {
        debug_assert!(self.ring().is_element(r));
        debug_assert!(self.is_element(v));
        if self.ring().is_zero(r) {
            [].into()
        } else {
            v.iter()
                .map(|(b, s)| (b.clone(), self.ring().mul(r, s)))
                .collect()
        }
    }
}

impl<Set: FreeModuleOverHashableSetElement, Ring: RingSignature, RingB: BorrowedStructure<Ring>>
    FreeModuleSignature<Ring> for FreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
}

#[cfg(test)]
mod tests {
    use algebraeon_nzq::Natural;

    use super::*;

    #[test]
    fn test_semimodule() {
        #[derive(Debug, Clone, PartialEq, Eq, Hash)]
        struct T(usize);

        let m = FreeModuleOverHashableSetStructure::<T, _, _>::new(Natural::structure());

        let v = [(T(5), Natural::from(2u32)), (T(7), Natural::from(3u32))].into();
        let w = [(T(5), Natural::from(1u32)), (T(10), Natural::from(4u32))].into();

        debug_assert!(
            m.equal(
                &m.add(&v, &w),
                &[
                    (T(5), Natural::from(3u32)),
                    (T(7), Natural::from(3u32)),
                    (T(10), Natural::from(4u32))
                ]
                .into()
            )
        );

        debug_assert!(m.equal(
            &m.scalar_mul(&Natural::from(3u32), &v),
            &[(T(5), Natural::from(6u32)), (T(7), Natural::from(9u32)),].into()
        ));
    }
}
