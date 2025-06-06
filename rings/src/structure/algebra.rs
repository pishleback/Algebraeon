use crate::structure::*;

pub trait AlgebraSignature<Ring: RingSignature>: ModuleSignature<Ring> + RingSignature {}
