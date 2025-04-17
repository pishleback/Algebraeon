use crate::structure::*;
use algebraeon_sets::structure::*;
use itertools::Itertools;
use std::{collections::HashMap, hash::Hash};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FreeModuleOverSetStructure<Set: SetSignature, Ring: RingSignature>
where
    Set::Set: Eq + Hash,
{
    set: Set,
    ring: Ring,
}

impl<Set: SetSignature, Ring: RingSignature> Signature for FreeModuleOverSetStructure<Set, Ring> where
    Set::Set: Eq + Hash
{
}

impl<Set: SetSignature, Ring: RingSignature> SetSignature for FreeModuleOverSetStructure<Set, Ring>
where
    Set::Set: Eq + Hash,
{
    type Set = HashMap<Set::Set, Ring::Set>;

    fn is_element(&self, v: &Self::Set) -> bool {
        v.iter()
            .all(|(e, r)| self.set.is_element(e) && self.ring.is_element(r))
    }
}

impl<Set: SetSignature, Ring: RingSignature> EqSignature for FreeModuleOverSetStructure<Set, Ring>
where
    Set::Set: Eq + Hash,
{
    fn equal(&self, v: &Self::Set, w: &Self::Set) -> bool {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        let mut keys = v.keys().chain(w.keys()).unique();
        keys.all(|k| match (v.get(k), w.get(k)) {
            (None, None) => unreachable!(),
            (None, Some(s)) => self.ring.is_zero(s),
            (Some(r), None) => self.ring.is_zero(r),
            (Some(r), Some(s)) => self.ring.equal(r, s),
        })
    }
}

impl<Set: SetSignature, Ring: RingSignature> ModuleSignature<Ring>
    for FreeModuleOverSetStructure<Set, Ring>
where
    Set::Set: Eq + Hash,
{
    fn ring(&self) -> &Ring {
        &self.ring
    }

    fn add(&self, v: &Self::Set, w: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        let keys = v.keys().chain(w.keys()).unique();
        keys.map(|k| {
            (
                k.clone(),
                match (v.get(k), w.get(k)) {
                    (None, None) => unreachable!(),
                    (None, Some(s)) => s.clone(),
                    (Some(r), None) => r.clone(),
                    (Some(r), Some(s)) => self.ring.add(r, s),
                },
            )
        })
        .collect()
    }

    fn neg(&self, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        v.iter()
            .map(|(b, r)| (b.clone(), self.ring.neg(r)))
            .collect()
    }

    fn scalar_mul(&self, r: &<Ring>::Set, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        v.iter()
            .map(|(b, s)| (b.clone(), self.ring.mul(r, s)))
            .collect()
    }
}

// #[derive(Debug, Clone, PartialEq, Eq)]
// pub struct FreeModuleFiniteNumberedBasisStructure<Ring: RingSignature> {
//     ring: Ring,
//     rank: usize,
// }
