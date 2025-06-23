use crate::structure::*;

/// A Lattice as in https://en.wikipedia.org/wiki/Lattice_(module)
pub trait LatticeSignature<
    Ring: IntegralDomainSignature,
    K: FieldSignature,
    KIsFOF: FieldOfFractionsInclusion<Ring, K>,
    A: AlgebraSignature<K>,
>: FinitelyGeneratedModuleSignature<Ring>
{
    fn is_full(&self) -> bool;
}

pub trait FullLatticeSignature<
    Ring: IntegralDomainSignature,
    K: FieldSignature,
    KIsFOF: FieldOfFractionsInclusion<Ring, K>,
    A: AlgebraSignature<K>,
>: LatticeSignature<Ring, K, KIsFOF, A>
{
}
