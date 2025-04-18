use algebraeon_sets::structure::*;
use super::{
    finitely_free_coset::FinitelyFreeSubmoduleCoset,
    finitely_free_modules::FinitelyFreeModuleStructure,
    matrix::{ReducedHermiteAlgorithmSignature, UniqueReducedHermiteAlgorithmSignature},
};

#[derive(Debug, Clone)]
pub enum FinitelyFreeSubmoduleAffineSubset<Ring: ReducedHermiteAlgorithmSignature> {
    Empty,
    NonEmpty(FinitelyFreeSubmoduleCoset<Ring>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeSubmoduleAffineSubsetStructure<
    Ring: ReducedHermiteAlgorithmSignature,
    const UNIQUE: bool,
> {
    module: FinitelyFreeModuleStructure<Ring>,
}

impl<Ring: ReducedHermiteAlgorithmSignature>
    FinitelyFreeSubmoduleAffineSubsetStructure<Ring, false>
{
    pub fn new_nonunique_reduction(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature>
    FinitelyFreeSubmoduleAffineSubsetStructure<Ring, true>
{
    pub fn new(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> Signature
    for FinitelyFreeSubmoduleAffineSubsetStructure<Ring, UNIQUE>
{
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> SetSignature
    for FinitelyFreeSubmoduleAffineSubsetStructure<Ring, UNIQUE>
{
    type Set = FinitelyFreeSubmoduleAffineSubset<Ring>;

    fn is_element(&self, sm: &Self::Set) -> bool {
        todo!()
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> EqSignature
    for FinitelyFreeSubmoduleAffineSubsetStructure<Ring, UNIQUE>
{
    fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
        todo!()
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool>
    FinitelyFreeSubmoduleAffineSubsetStructure<Ring, UNIQUE>
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;
}
