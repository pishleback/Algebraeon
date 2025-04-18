use super::{
    finitely_free_coset::FinitelyFreeSubmoduleCoset,
    finitely_free_modules::FinitelyFreeModuleStructure,
    matrix::{ReducedHermiteAlgorithmSignature, UniqueReducedHermiteAlgorithmSignature},
};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub enum FinitelyFreeSubmoduleAffineSubset<Ring: ReducedHermiteAlgorithmSignature> {
    Empty,
    NonEmpty(FinitelyFreeSubmoduleCoset<Ring>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeSubmoduleAffineSubsetsStructure<
    Ring: ReducedHermiteAlgorithmSignature,
    const UNIQUE: bool,
> {
    module: FinitelyFreeModuleStructure<Ring>,
}

impl<Ring: ReducedHermiteAlgorithmSignature>
    FinitelyFreeSubmoduleAffineSubsetsStructure<Ring, false>
{
    pub fn new_nonunique_reduction(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature>
    FinitelyFreeSubmoduleAffineSubsetsStructure<Ring, true>
{
    pub fn new(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> Signature
    for FinitelyFreeSubmoduleAffineSubsetsStructure<Ring, UNIQUE>
{
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> SetSignature
    for FinitelyFreeSubmoduleAffineSubsetsStructure<Ring, UNIQUE>
{
    type Set = FinitelyFreeSubmoduleAffineSubset<Ring>;

    fn is_element(&self, sm: &Self::Set) -> bool {
        todo!()
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> EqSignature
    for FinitelyFreeSubmoduleAffineSubsetsStructure<Ring, UNIQUE>
{
    fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
        todo!()
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool>
    FinitelyFreeSubmoduleAffineSubsetsStructure<Ring, UNIQUE>
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;
}
