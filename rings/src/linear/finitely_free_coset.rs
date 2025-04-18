use super::{
    finitely_free_modules::FinitelyFreeModuleStructure,
    finitely_free_submodule::FinitelyFreeSubmodule,
    matrix::{ReducedHermiteAlgorithmSignature, UniqueReducedHermiteAlgorithmSignature},
};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmoduleCoset<Ring: ReducedHermiteAlgorithmSignature> {
    offset: Vec<Ring::Set>,
    submodule: FinitelyFreeSubmodule<Ring>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeSubmoduleCosetsStructure<
    Ring: ReducedHermiteAlgorithmSignature,
    const UNIQUE: bool,
> {
    module: FinitelyFreeModuleStructure<Ring>,
}

impl<Ring: ReducedHermiteAlgorithmSignature> FinitelyFreeSubmoduleCosetsStructure<Ring, false> {
    pub fn new_nonunique_reduction(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: UniqueReducedHermiteAlgorithmSignature>
    FinitelyFreeSubmoduleCosetsStructure<Ring, true>
{
    pub fn new(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> Signature
    for FinitelyFreeSubmoduleCosetsStructure<Ring, UNIQUE>
{
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> SetSignature
    for FinitelyFreeSubmoduleCosetsStructure<Ring, UNIQUE>
{
    type Set = FinitelyFreeSubmoduleCoset<Ring>;

    fn is_element(&self, sm: &Self::Set) -> bool {
        todo!()
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool> EqSignature
    for FinitelyFreeSubmoduleCosetsStructure<Ring, UNIQUE>
{
    fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
        todo!()
    }
}

impl<Ring: ReducedHermiteAlgorithmSignature, const UNIQUE: bool>
    FinitelyFreeSubmoduleCosetsStructure<Ring, UNIQUE>
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;
}
