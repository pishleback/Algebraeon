use crate::structure::*;
use algebraeon_sets::structure::*;

pub trait ModuleSignature<Ring: RingSignature>: SetSignature {
    fn ring(&self) -> &Ring;
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set;
    fn neg(&self, a: &Self::Set) -> Self::Set;
    fn scalar_mul(&self, x: &Ring::Set, a: &Self::Set) -> Self::Set;
}

pub trait FreeModuleSignature<Ring: RingSignature>: ModuleSignature<Ring> {}

pub trait FiniteRankModuleSignature<Ring: RingSignature>: FreeModuleSignature<Ring> {
    fn basis(&self) -> Vec<Self::Set>;
    fn rank(&self) -> usize {
        self.basis().len()
    }
}

pub trait LinearTransformation<
    Ring: RingSignature,
    Domain: ModuleSignature<Ring>,
    Range: ModuleSignature<Ring>,
>: Function<Domain, Range>
{
}