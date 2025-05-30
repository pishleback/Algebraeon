use crate::structure::*;
use algebraeon_sets::structure::*;
use std::fmt::Debug;
use std::{collections::HashMap, hash::Hash, marker::PhantomData};

use super::hashable_set_free_module::FreeModuleOverHashableSetStructure;

pub trait FinitelyFreeModuleOverHashableSetElement: Debug + Clone + Eq + Hash {}
impl<Set: Debug + Clone + Eq + Hash> FinitelyFreeModuleOverHashableSetElement for Set {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeModuleOverHashableSetStructure<
    Set: FinitelyFreeModuleOverHashableSetElement,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> {
    _ring: PhantomData<Ring>,
    basis: Vec<Set>,
    ring: RingB,
}

pub trait RingToFinitelyFreeModuleOverHashableSetStructure: SemiRingSignature {
    fn finitely_free_module_on_hashable_set<'a, Set: FinitelyFreeModuleOverHashableSetElement>(
        &'a self,
        basis: Vec<Set>,
    ) -> FinitelyFreeModuleOverHashableSetStructure<Set, Self, &'a Self> {
        FinitelyFreeModuleOverHashableSetStructure::new(self, basis)
    }
    fn into_finitely_free_module_on_hashable_set<Set: FinitelyFreeModuleOverHashableSetElement>(
        self,
        basis: Vec<Set>,
    ) -> FinitelyFreeModuleOverHashableSetStructure<Set, Self, Self> {
        FinitelyFreeModuleOverHashableSetStructure::new(self, basis)
    }
}
impl<Ring: SemiRingSignature> RingToFinitelyFreeModuleOverHashableSetStructure for Ring {}

impl<
    Set: FinitelyFreeModuleOverHashableSetElement,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> FinitelyFreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    pub fn new(ring: RingB, basis: Vec<Set>) -> Self {
        Self {
            _ring: PhantomData::default(),
            basis,
            ring,
        }
    }

    pub fn into_free_module(self) -> FreeModuleOverHashableSetStructure<Set, Ring, RingB> {
        FreeModuleOverHashableSetStructure::new(self.ring)
    }

    pub fn free_module<'a>(&'a self) -> FreeModuleOverHashableSetStructure<Set, Ring, &'a Ring> {
        FreeModuleOverHashableSetStructure::new(self.ring.borrow())
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

impl<
    Set: FinitelyFreeModuleOverHashableSetElement,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> Signature for FinitelyFreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
}

impl<
    Set: FinitelyFreeModuleOverHashableSetElement,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> SetSignature for FinitelyFreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    type Set = HashMap<Set, Ring::Set>;

    fn is_element(&self, v: &Self::Set) -> bool {
        v.iter()
            .all(|(b, r)| self.basis.contains(b) && self.ring().is_element(r))
    }
}

impl<
    Set: FinitelyFreeModuleOverHashableSetElement,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> EqSignature for FinitelyFreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    fn equal(&self, v: &Self::Set, w: &Self::Set) -> bool {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        self.free_module().equal(v, w)
    }
}

impl<
    Set: FinitelyFreeModuleOverHashableSetElement,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditiveMonoidSignature for FinitelyFreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    fn zero(&self) -> Self::Set {
        self.free_module().zero()
    }

    fn add(&self, v: &Self::Set, w: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        self.free_module().add(v, w)
    }
}

impl<
    Set: FinitelyFreeModuleOverHashableSetElement,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> AdditiveGroupSignature for FinitelyFreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    fn neg(&self, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        self.free_module().neg(v)
    }
}

impl<
    Set: FinitelyFreeModuleOverHashableSetElement,
    Ring: SemiRingSignature,
    RingB: BorrowedStructure<Ring>,
> SemiModuleSignature<Ring> for FinitelyFreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    fn ring(&self) -> &Ring {
        self.ring.borrow()
    }

    fn scalar_mul(&self, r: &<Ring>::Set, v: &Self::Set) -> Self::Set {
        debug_assert!(self.ring().is_element(r));
        debug_assert!(self.is_element(v));
        self.free_module().scalar_mul(r, v)
    }
}

impl<
    Set: FinitelyFreeModuleOverHashableSetElement,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> FreeModuleSignature<Ring> for FinitelyFreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
}

impl<
    Set: FinitelyFreeModuleOverHashableSetElement,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
> FinitelyFreeModuleSignature<Ring>
    for FinitelyFreeModuleOverHashableSetStructure<Set, Ring, RingB>
{
    fn rank(&self) -> usize {
        self.basis.len()
    }

    fn to_vec(&self, v: &Self::Set) -> Vec<Ring::Set> {
        debug_assert!(self.is_element(v));
        self.basis
            .iter()
            .map(|b| match v.get(&b) {
                Some(c) => c.clone(),
                None => self.ring().zero(),
            })
            .collect()
    }

    fn from_vec(&self, v: &Vec<Ring::Set>) -> Self::Set {
        self.reduce(
            self.basis
                .iter()
                .enumerate()
                .map(|(i, b)| (b.clone(), v[i].clone()))
                .collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::{Integer, Natural};

    #[test]
    fn test_finite_rank_modules() {
        #[derive(Debug, Clone, PartialEq, Eq, Hash)]
        enum Basis {
            A,
            B,
            C,
        }

        let m = FinitelyFreeModuleOverHashableSetStructure::<Basis, _, _>::new(
            Integer::structure(),
            vec![Basis::A, Basis::B, Basis::C],
        );

        let a = m.basis_element(Basis::A);
        let b = m.basis_element(Basis::B);
        let c = m.basis_element(Basis::C);

        assert_eq!(
            m.add(&m.neg(&b), &m.add(&a, &b)),
            [(Basis::A, Integer::ONE)].into()
        );

        assert_eq!(
            m.add(&m.add(&a, &b), &m.add(&b, &c)),
            [
                (Basis::A, Integer::ONE),
                (Basis::B, Integer::TWO),
                (Basis::C, Integer::ONE)
            ]
            .into()
        );

        assert_eq!(m.scalar_mul(&5.into(), &a), [(Basis::A, 5.into())].into());

        assert_eq!(m.basis(), vec![a, b, c]);
    }

    #[test]
    fn test_semimodule() {
        #[derive(Debug, Clone, PartialEq, Eq, Hash)]
        struct T(usize);

        let m = Natural::structure().into_finitely_free_module_on_hashable_set(vec![
            T(5),
            T(7),
            T(10),
            T(4),
        ]);

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
