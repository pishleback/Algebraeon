use crate::structure::*;
use algebraeon_sets::structure::FiniteSetSignature;

/// Algebras over `Ring`
pub trait AlgebraSignature<Ring: RingSignature>: ModuleSignature<Ring> + RingSignature {}

pub trait FiniteDimensionalAlgebraSignature<Field: FieldSignature>:
    AlgebraSignature<Field> + FinitelyFreeModuleSignature<Field>
where
    Self::Basis: FiniteSetSignature,
{
}

// /// Order in a finite dimensional algebra as in <https://en.wikipedia.org/wiki/Order_(ring_theory)>
// pub trait OrderSignature<
//     Ring: IntegralDomainSignature,
//     K: FieldSignature,
//     KIsFOF: FieldOfFractionsInclusion<Ring, K>,
//     A: FiniteDimensionalAlgebraSignature<K>,
// >: FullLatticeSignature<Ring, K, KIsFOF, A> + RingSignature where
//     A::Basis: FiniteSetSignature,
// {
//     fn is_maximal(&self) -> bool;
// }
