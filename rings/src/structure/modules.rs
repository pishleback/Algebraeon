use crate::structure::*;
use algebraeon_sets::structure::*;

pub trait ModuleSignature<Ring: RingSignature>: SetSignature {
    fn ring(&self) -> &Ring;
    fn zero(&self) -> Self::Set;
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set;
    fn neg(&self, a: &Self::Set) -> Self::Set;
    fn scalar_mul(&self, x: &Ring::Set, a: &Self::Set) -> Self::Set;
}

pub trait FreeModuleSignature<Ring: RingSignature>: ModuleSignature<Ring> {}

pub trait FinitelyFreeModuleSignature<Ring: RingSignature>: FreeModuleSignature<Ring> {
    fn basis(&self) -> Vec<Self::Set> {
        (0..self.rank())
            .map(|j| {
                self.from_vec(
                    &(0..self.rank())
                        .map(|i| {
                            if i == j {
                                self.ring().one()
                            } else {
                                self.ring().zero()
                            }
                        })
                        .collect(),
                )
            })
            .collect()
    }
    fn rank(&self) -> usize;
    fn to_vec(&self, v: &Self::Set) -> Vec<Ring::Set>;
    fn from_vec(&self, v: &Vec<Ring::Set>) -> Self::Set {
        debug_assert_eq!(v.len(), self.rank());
        let mut t = self.zero();
        for (i, b) in self.basis().into_iter().enumerate() {
            t = self.add(&t, &self.scalar_mul(&v[i], &b));
        }
        t
    }
}

pub trait LinearTransformation<
    Ring: RingSignature,
    Domain: ModuleSignature<Ring>,
    Range: ModuleSignature<Ring>,
>: Function<Domain, Range>
{
}
