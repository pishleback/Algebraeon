use crate::structure::*;
use algebraeon_sets::structure::*;
use std::marker::PhantomData;

use super::matrix::Matrix;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmoduleStructure<
    Ring: EqSignature + BezoutDomainSignature,
    Module: FinitelyFreeModuleSignature<Ring>,
> {
    ring: PhantomData<Ring>,
    module: Module,
    // linearly independent rows
    basis_matrix: Matrix<Ring::Set>,
}

impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
    FinitelyFreeSubmoduleStructure<Ring, Module>
{
    pub fn ring(&self) -> &Ring {
        self.module.ring()
    }

    pub fn module(&self) -> &Module {
        &self.module
    }

    // pub fn from_span(module: Module, span: Vec<Module::Set>) -> Self {
    //     for v in span {
    //         debug_assert!(module.is_element(&v));
    //     }
    //     module
    //     todo!()
    // }

    // pub fn from_basis(module: Module, basis: Vec<Module::Set>) -> Self {
    //     for v in basis {
    //         debug_assert!(module.is_element(&v));
    //     }
    //     todo!()
    // }

    // pub fn submodule(
    //     &self,
    // ) -> impl InjectiveFunction<FinitelyFreeModuleStructure<Ring>, Module>
    // + LinearTransformation<Ring, FinitelyFreeModuleStructure<Ring>, Module> {
    //     todo!()
    // }
}

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>> Signature
//     for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     SetSignature for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
//     type Set = Vec<Ring::Set>;

//     fn is_element(&self, x: &Self::Set) -> bool {
//         self.submodule().is
//     }
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     EqSignature for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
//     fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
//         todo!()
//     }
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     ModuleSignature<Ring> for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
//     fn ring(&self) -> &Ring {
//         todo!()
//     }

//     fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
//         todo!()
//     }

//     fn neg(&self, a: &Self::Set) -> Self::Set {
//         todo!()
//     }

//     fn scalar_mul(&self, x: &<Ring>::Set, a: &Self::Set) -> Self::Set {
//         todo!()
//     }
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     FreeModuleSignature<Ring> for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     FinitelyFreeModuleSignature<Ring> for FinitelyFreeSubmoduleStructure<Ring, Module>
// {
//     fn rank(&self) -> usize {
//         self.basis_matrix.rows()
//     }

//     fn basis(&self) -> Vec<Self::Set> {
//         self.submodule().basis()
//     }
// }
