use crate::structure::*;
use algebraeon_sets::structure::*;
use std::fmt::Debug;
use std::marker::PhantomData;

pub trait FreeModuleOverSetElement: Debug + Clone + Eq {}
impl<Set: Debug + Clone + Eq> FreeModuleOverSetElement for Set {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FreeModuleOverSetStructure<
    Set: FreeModuleOverSetElement,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> {
    _set: PhantomData<Set>,
    _ring: PhantomData<Ring>,
    ring: RingB,
}

pub trait RingToFreeModuleOverSetStructure: SemiRingSignature {
    fn free_module_on_set<'a, Set: FreeModuleOverSetElement>(
        &'a self,
    ) -> FreeModuleOverSetStructure<Set, Self, &'a Self> {
        FreeModuleOverSetStructure::new(self)
    }
    fn into_free_module_on_set<Set: FreeModuleOverSetElement>(
        self,
    ) -> FreeModuleOverSetStructure<Set, Self, Self> {
        FreeModuleOverSetStructure::new(self)
    }
}
impl<Ring: SemiRingSignature> RingToFreeModuleOverSetStructure for Ring {}

impl<Set: FreeModuleOverSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    FreeModuleOverSetStructure<Set, Ring, RingB>
{
    pub fn new(ring: RingB) -> Self {
        Self {
            _set: PhantomData::default(),
            _ring: PhantomData::default(),
            ring,
        }
    }

    pub fn reduce(&self, v: &<Self as SetSignature>::Set) -> <Self as SetSignature>::Set {
        let mut v_reduced = vec![];
        for (b, c) in v {
            'SEARCH: {
                for (existing_b, existing_c) in &mut v_reduced {
                    if b == existing_b {
                        self.ring().add_mut(existing_c, c);
                        break 'SEARCH;
                    }
                }
                v_reduced.push((b.clone(), c.clone()));
            }
        }
        v_reduced
            .into_iter()
            .filter(|(_, c)| !self.ring().is_zero(c))
            .collect()
    }

    pub fn basis_element(&self, b: Set) -> <Self as SetSignature>::Set {
        vec![(b, self.ring().one())]
    }
}

impl<Set: FreeModuleOverSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    Signature for FreeModuleOverSetStructure<Set, Ring, RingB>
{
}

impl<Set: FreeModuleOverSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    SetSignature for FreeModuleOverSetStructure<Set, Ring, RingB>
{
    type Set = Vec<(Set, Ring::Set)>;

    fn is_element(&self, v: &Self::Set) -> bool {
        // duplicates in Set field allowed
        v.iter().all(|(_, r)| self.ring().is_element(r))
    }
}

impl<Set: FreeModuleOverSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    EqSignature for FreeModuleOverSetStructure<Set, Ring, RingB>
{
    fn equal(&self, v: &Self::Set, w: &Self::Set) -> bool {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));

        let v = self.reduce(v);
        let w = self.reduce(w);

        //every non-zero component of v is in w
        for (vb, vc) in &v {
            if !self.ring().is_zero(vc) {
                if !w.iter().any(|(wb, _)| vb == wb) {
                    return false;
                }
            }
        }
        //every non-zero component of w is in v
        for (wb, wc) in &w {
            if !self.ring().is_zero(wc) {
                if !v.iter().any(|(vb, _)| vb == wb) {
                    return false;
                }
            }
        }
        //the components of all basis elements are equal
        for (vb, vc) in &v {
            for (wb, wc) in &w {
                if vb == wb {
                    if !self.ring().equal(vc, wc) {
                        return false;
                    }
                }
            }
        }
        true
    }
}

impl<Set: FreeModuleOverSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    AdditiveMonoidSignature for FreeModuleOverSetStructure<Set, Ring, RingB>
{
    fn zero(&self) -> Self::Set {
        vec![]
    }

    fn add(&self, v: &Self::Set, w: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        self.reduce(&v.clone().into_iter().chain(w.clone().into_iter()).collect())
    }
}

impl<Set: FreeModuleOverSetElement, Ring: RingSignature, RingB: BorrowedStructure<Ring>>
    AdditiveGroupSignature for FreeModuleOverSetStructure<Set, Ring, RingB>
{
    fn neg(&self, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        v.iter()
            .map(|(b, r)| (b.clone(), self.ring().neg(r)))
            .collect()
    }
}

impl<Set: FreeModuleOverSetElement, Ring: SemiRingSignature, RingB: BorrowedStructure<Ring>>
    SemiModuleSignature<Ring> for FreeModuleOverSetStructure<Set, Ring, RingB>
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

impl<Set: FreeModuleOverSetElement, Ring: RingSignature, RingB: BorrowedStructure<Ring>>
    FreeModuleSignature<Ring> for FreeModuleOverSetStructure<Set, Ring, RingB>
{
    type Basis = Set;

    fn to_component(&self, b: &Self::Basis, v: &Self::Set) -> Ring::Set {
        let v = self.reduce(v);
        for (vb, vc) in v {
            if *b == vb {
                return vc;
            }
        }
        self.ring().zero()
    }

    fn from_component(&self, b: &Self::Basis, r: &<Ring>::Set) -> Self::Set {
        vec![(b.clone(), r.clone())]
    }
}

#[cfg(test)]
mod tests {
    use algebraeon_nzq::Natural;

    use super::*;

    #[test]
    fn test_semimodule() {
        #[derive(Debug, Clone, PartialEq, Eq, Hash)]
        struct T(usize);

        let m = Natural::structure().into_free_module_on_set();

        let v = [(T(5), Natural::from(2u32)), (T(7), Natural::from(3u32))].into();
        let w = [(T(5), Natural::from(1u32)), (T(10), Natural::from(4u32))].into();

        debug_assert!(
            m.equal(
                &m.add(&v, &w),
                &[
                    (T(5), Natural::from(3u32)),
                    (T(7), Natural::from(3u32)),
                    (T(10), Natural::from(1u32)),
                    (T(10), Natural::from(3u32)),
                    (T(6), Natural::from(0u32)),
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
