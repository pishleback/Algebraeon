use crate::structure::*;
use algebraeon_sets::structure::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FreeModuleFiniteNumberedBasisStructure<Ring: RingSignature> {
    ring: Ring,
    rank: usize,
}

impl<Ring: RingSignature> FreeModuleFiniteNumberedBasisStructure<Ring> {
    pub fn new(ring: Ring, rank: usize) -> Self {
        Self { ring, rank }
    }

    pub fn basis_element(&self, i: usize) -> <Self as SetSignature>::Set {
        debug_assert!(i < self.rank);
        (0..self.rank)
            .map(|j| {
                if i == j {
                    self.ring.one()
                } else {
                    self.ring.zero()
                }
            })
            .collect()
    }
}

impl<Ring: RingSignature> Signature for FreeModuleFiniteNumberedBasisStructure<Ring> {}

impl<Ring: RingSignature> SetSignature for FreeModuleFiniteNumberedBasisStructure<Ring> {
    type Set = Vec<Ring::Set>;

    fn is_element(&self, v: &Self::Set) -> bool {
        self.rank == v.len() && v.iter().all(|r| self.ring.is_element(r))
    }
}

impl<Ring: RingSignature> EqSignature for FreeModuleFiniteNumberedBasisStructure<Ring> {
    fn equal(&self, v: &Self::Set, w: &Self::Set) -> bool {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        (0..self.rank).all(|i| self.ring.equal(&v[i], &w[i]))
    }
}

impl<Ring: RingSignature> ModuleSignature<Ring> for FreeModuleFiniteNumberedBasisStructure<Ring> {
    fn ring(&self) -> &Ring {
        &self.ring
    }

    fn add(&self, v: &Self::Set, w: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        debug_assert!(self.is_element(w));
        (0..self.rank)
            .map(|i| self.ring.add(&v[i], &w[i]))
            .collect()
    }

    fn neg(&self, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        v.iter().map(|r| self.ring.neg(r)).collect()
    }

    fn scalar_mul(&self, r: &<Ring>::Set, v: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(v));
        v.iter().map(|s| self.ring.mul(r, s)).collect()
    }
}

impl<Ring: RingSignature> FreeModuleSignature<Ring>
    for FreeModuleFiniteNumberedBasisStructure<Ring>
{
}

impl<Ring: RingSignature> FiniteRankModuleSignature<Ring>
    for FreeModuleFiniteNumberedBasisStructure<Ring>
{
    fn basis(&self) -> Vec<Self::Set> {
        (0..self.rank)
            .into_iter()
            .map(|i| self.basis_element(i))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;

    #[test]
    fn test_finite_rank_modules() {
        let m = FreeModuleFiniteNumberedBasisStructure::new(Integer::structure(), 3);

        let a = m.basis_element(0);
        let b = m.basis_element(1);
        let c = m.basis_element(2);

        assert_eq!(
            m.add(&m.neg(&b), &m.add(&a, &b)),
            vec![Integer::from(1), Integer::from(0), Integer::from(0)]
        );

        assert_eq!(
            m.add(&m.add(&a, &b), &m.add(&b, &c)),
            vec![Integer::from(1), Integer::from(2), Integer::from(1)]
        );

        assert_eq!(
            m.scalar_mul(&5.into(), &a),
            vec![Integer::from(5), Integer::from(0), Integer::from(0)]
        );

        assert_eq!(m.basis(), vec![a, b, c]);
    }
}
