use crate::structure::*;
use algebraeon_sets::structure::*;
use itertools::Itertools;
// use itertools::Itertools;
use std::{collections::HashMap, hash::Hash};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FreeModuleOverSetStructure<Set: SetSignature, Ring: RingSignature>
where
    Set::Set: Eq + Hash,
{
    set: Set,
    ring: Ring,
}

impl<Set: SetSignature, Ring: RingSignature> FreeModuleOverSetStructure<Set, Ring>
where
    Set::Set: Eq + Hash,
{
    pub fn new(set: Set, ring: Ring) -> Self {
        Self { set, ring }
    }

    pub fn reduce(&self, v: <Self as SetSignature>::Set) -> <Self as SetSignature>::Set {
        v.into_iter()
            .filter(|(_, r)| !self.ring.is_zero(r))
            .collect()
    }

    pub fn basis_element(&self, b: Set::Set) -> <Self as SetSignature>::Set {
        debug_assert!(self.set.is_element(&b));
        [(b, self.ring.one())].into()
    }
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
                        (Some(r), Some(s)) => self.ring.add(r, s),
                    },
                )
            })
            .collect(),
        )
    }

    fn neg(&self, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        v.iter()
            .map(|(b, r)| (b.clone(), self.ring.neg(r)))
            .collect()
    }

    fn scalar_mul(&self, r: &<Ring>::Set, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        if self.ring.is_zero(r) {
            [].into()
        } else {
            v.iter()
                .map(|(b, s)| (b.clone(), self.ring.mul(r, s)))
                .collect()
        }
    }
}

impl<Set: SetSignature, Ring: RingSignature> FreeModuleSignature<Ring>
    for FreeModuleOverSetStructure<Set, Ring>
where
    Set::Set: Eq + Hash,
{
}

impl<Set: FiniteSetSignature, Ring: RingSignature> FinitelyFreeModuleSignature<Ring>
    for FreeModuleOverSetStructure<Set, Ring>
where
    Set::Set: Eq + Hash,
{
    fn rank(&self) -> usize {
        self.set.size()
    }

    fn to_vec(&self, v: &Self::Set) -> Vec<Ring::Set> {
        self.set
            .list_all_elements()
            .into_iter()
            .map(|b| match v.get(&b) {
                Some(c) => c.clone(),
                None => self.ring.zero(),
            })
            .collect()
    }

    fn from_vec(&self, v: &Vec<Ring::Set>) -> Self::Set {
        self.set
            .list_all_elements()
            .into_iter()
            .enumerate()
            .map(|(i, b)| (b, v[i].clone()))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;

    #[test]
    fn test_finite_rank_modules() {
        #[derive(Debug, Clone, PartialEq, Eq, Hash, CanonicalStructure)]
        enum Basis {
            A,
            B,
            C,
        }

        impl CountableSetSignature for BasisCanonicalStructure {
            fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> {
                vec![Basis::A, Basis::B, Basis::C].into_iter()
            }
        }

        impl FiniteSetSignature for BasisCanonicalStructure {}

        let m = FreeModuleOverSetStructure::new(Basis::structure(), Integer::structure());

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
}
