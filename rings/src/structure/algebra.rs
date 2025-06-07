use crate::structure::*;

/// Algebras over `Ring`
pub trait AlgebraSignature<Ring: RingSignature>: ModuleSignature<Ring> + RingSignature {
}
