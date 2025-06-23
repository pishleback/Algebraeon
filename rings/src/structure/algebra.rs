use crate::structure::*;
use algebraeon_sets::structure::FiniteSetSignature;

/// Algebras over `Ring`
pub trait AlgebraSignature<Ring: RingSignature>: ModuleSignature<Ring> + RingSignature {}

pub trait FiniteDimensionalAlgebraSignature<Basis: FiniteSetSignature, Field: FieldSignature>:
    AlgebraSignature<Field> + FinitelyFreeModuleSignature<Basis, Field>
{
}

/// Order in a finite dimensional algebra as in https://en.wikipedia.org/wiki/Order_(ring_theory)
pub trait OrderSignature<
    Basis: FiniteSetSignature,
    Ring: IntegralDomainSignature,
    K: FieldSignature,
    KIsFOF: FieldOfFractionsInclusion<Ring, K>,
    A: FiniteDimensionalAlgebraSignature<Basis, K>,
>: FullLatticeSignature<Ring, K, KIsFOF, A> + RingSignature
{
    fn is_maximal(&self) -> bool;
}
