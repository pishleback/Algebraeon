use crate::structure::*;

/// Algebras over `Ring`
pub trait AlgebraSignature<Ring: RingSignature>: ModuleSignature<Ring> + RingSignature {}

pub trait FiniteDimensionalAlgebraSignature<K: FieldSignature>:
    AlgebraSignature<K> + FinitelyFreeModuleSignature<K>
where
    <Self as modules::FreeModuleSignature<K>>::Basis:
        algebraeon_sets::structure::FiniteSetSignature,
{
}

/// Order in a finite dimensional algebra as in https://en.wikipedia.org/wiki/Order_(ring_theory)
pub trait OrderSignature<
    Ring: IntegralDomainSignature,
    K: FieldSignature,
    KIsFOF: FieldOfFractionsInclusion<Ring, K>,
    A: FiniteDimensionalAlgebraSignature<K>,
>: FullLatticeSignature<Ring, K, KIsFOF, A> + RingSignature where
    <A as modules::FreeModuleSignature<K>>::Basis: algebraeon_sets::structure::FiniteSetSignature,
{
    fn is_maximal(&self) -> bool;
}
