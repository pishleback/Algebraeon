use crate::structure::*;
use algebraeon_sets::structure::*;

pub trait SemiModuleSignature<Ring: SemiRingSignature>: AdditiveMonoidSignature {
    fn ring(&self) -> &Ring;
    fn scalar_mul(&self, x: &Ring::Set, a: &Self::Set) -> Self::Set;
}

pub trait ModuleSignature<Ring: RingSignature>:
    SemiModuleSignature<Ring> + AdditiveGroupSignature
{
}
impl<Ring: RingSignature, Module: SemiModuleSignature<Ring> + AdditiveGroupSignature>
    ModuleSignature<Ring> for Module
{
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

// pub trait SubModuleSignature<Ring: RingSignature, Module: ModuleSignature<Ring>>:
//     EqSignature
// {
//     fn ring(&self) -> &Ring;
//     fn module(&self) -> &Module;
//     fn zero_submodule(&self) -> Self::Set {
//         self.generated(vec![])
//     }
//     fn improper_submodule(&self) -> Self::Set;
//     fn add(&self, x: &Self::Set, y: &Self::Set) -> Self::Set;
//     fn intersect(&self, x: &Self::Set, y: &Self::Set) -> Self::Set;
//     fn generated(&self, generators: Vec<&Module::Set>) -> Self::Set;
//     /// Does x contain p
//     fn contains_element(&self, x: &Self::Set, p: &Module::Set) -> bool;
//     /// Does x contain y
//     fn contains(&self, x: &Self::Set, y: &Self::Set) -> bool;
// }

// pub trait SubModuleCosetSignature<Ring: RingSignature, Module: ModuleSignature<Ring>>:
//     EqSignature
// {
//     fn ring(&self) -> &Ring;
//     fn module(&self) -> &Module;
//     fn add(&self, x: &Self::Set, y: &Self::Set) -> Self::Set;
//     fn intersect(&self, x: &Self::Set, y: &Self::Set) -> Self::Set;
// }
