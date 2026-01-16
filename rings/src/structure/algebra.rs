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
