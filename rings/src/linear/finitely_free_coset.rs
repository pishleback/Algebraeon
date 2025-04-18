use super::{
    finitely_free_modules::FinitelyFreeModuleStructure,
    finitely_free_submodule::FinitelyFreeSubmodule,
    matrix::{ReducedHermiteAlgorithmSignature, UniqueReducedHermiteAlgorithmSignature},
};
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmoduleCoset<Ring: ReducedHermiteAlgorithmSignature> {
    // reduced coset offset
    offset: Vec<Ring::Set>,
    submodule: FinitelyFreeSubmodule<Ring>,
}

impl<Ring: ReducedHermiteAlgorithmSignature> FinitelyFreeSubmoduleCoset<Ring> {
    pub fn ring(&self) -> &Ring {
        self.submodule.ring()
    }

    pub fn module(&self) -> FinitelyFreeModuleStructure<Ring> {
        self.submodule.module()
    }

    pub fn module_rank(&self) -> usize {
        self.submodule.module_rank()
    }

    pub fn coset_rank(&self) -> usize {
        self.submodule.submodule_rank()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Integer;
}
