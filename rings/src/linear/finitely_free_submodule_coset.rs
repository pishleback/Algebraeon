use super::{finitely_free_modules::FinitelyFreeModuleStructure, matrix::Matrix};
use crate::{linear::matrix::MatrixStructure, structure::*};
use algebraeon_sets::structure::*;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmoduleCoset<Ring: BezoutDomainSignature> {
    ring: PhantomData<Ring>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeSubmoduleCosetStructure<Ring: BezoutDomainSignature> {
    module: FinitelyFreeModuleStructure<Ring>,
}

impl<Ring: BezoutDomainSignature> FinitelyFreeSubmoduleCosetStructure<Ring>
where
    Ring: ToStringSignature,
{
    pub fn new(module: FinitelyFreeModuleStructure<Ring>) -> Self {
        Self { module }
    }
}

impl<Ring: BezoutDomainSignature> Signature for FinitelyFreeSubmoduleCosetStructure<Ring> where
    Ring: ToStringSignature
{
}

impl<Ring: BezoutDomainSignature> SetSignature for FinitelyFreeSubmoduleCosetStructure<Ring>
where
    Ring: ToStringSignature,
{
    type Set = FinitelyFreeSubmoduleCoset<Ring>;

    fn is_element(&self, sm: &Self::Set) -> bool {
        todo!()
    }
}

impl<Ring: BezoutDomainSignature> EqSignature for FinitelyFreeSubmoduleCosetStructure<Ring>
where
    Ring: ToStringSignature,
{
    fn equal(&self, x: &Self::Set, y: &Self::Set) -> bool {
        todo!()
    }
}

impl<Ring: BezoutDomainSignature> SubModuleCosetSignature<Ring, FinitelyFreeModuleStructure<Ring>>
    for FinitelyFreeSubmoduleCosetStructure<Ring>
where
    Ring: ToStringSignature,
{
    fn ring(&self) -> &Ring {
        self.module.ring()
    }

    fn module(&self) -> &FinitelyFreeModuleStructure<Ring> {
        &self.module
    }

    fn add(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&x));
        debug_assert!(self.is_element(&y));
        todo!()
    }

    fn intersect(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&x));
        debug_assert!(self.is_element(&y));
        todo!()
    }
}
