use crate::structure::*;
use algebraeon_structures::*;

/// Algebras over `Ring`
pub trait AlgebraSignature<Ring: RingSignature>: ModuleSignature<Ring> + RingSignature {}

pub trait FiniteDimensionalAlgebraSignature<Basis: FiniteSetSignature, Field: FieldSignature>:
    AlgebraSignature<Field> + FinitelyFreeModuleSignature<Basis, Field>
{
}
