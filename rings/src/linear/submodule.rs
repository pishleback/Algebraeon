use super::{finitely_free_modules::FinitelyFreeModuleStructure, matrix::Matrix};
use crate::structure::*;
use algebraeon_nzq::IntegerCanonicalStructure;
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct SubmoduleOfFinitelyFreeModule<Ring: RingSignature> {
    // rows are a basis for the submodule
    basis: Matrix<Ring::Set>,
}

pub trait SubmoduleOfFinitelyFreeModuleSignature<Ring: RingSignature>:
    SetSignature<Set = SubmoduleOfFinitelyFreeModule<Ring>>
{
    fn module(&self) -> &FinitelyFreeModuleStructure<Ring>;
    fn ring(&self) -> &Ring {
        self.module().ring()
    }
}

pub trait SubmoduleOfFinitelyFreeModuleUniqueBasisSignature<Ring: RingSignature>:
    SubmoduleOfFinitelyFreeModuleSignature<Ring>
{
    fn unique_reduce(&self, basis: Self::Set) -> Self::Set;
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SubmoduleOfFinitelyFreeModuleStructure<Ring: BezoutDomainSignature> {
    module: FinitelyFreeModuleStructure<Ring>,
}

impl<Ring: BezoutDomainSignature> Signature for SubmoduleOfFinitelyFreeModuleStructure<Ring> {}

impl<Ring: BezoutDomainSignature> SetSignature for SubmoduleOfFinitelyFreeModuleStructure<Ring> {
    type Set = SubmoduleOfFinitelyFreeModule<Ring>;

    fn is_element(&self, x: &Self::Set) -> bool {
        self.module.rank() == x.basis.cols()
    }
}

impl<Ring: BezoutDomainSignature> SubmoduleOfFinitelyFreeModuleSignature<Ring>
    for SubmoduleOfFinitelyFreeModuleStructure<Ring>
{
    fn module(&self) -> &FinitelyFreeModuleStructure<Ring> {
        &self.module
    }
}

impl<Field: FieldSignature> SubmoduleOfFinitelyFreeModuleUniqueBasisSignature<Field>
    for SubmoduleOfFinitelyFreeModuleStructure<Field>
{
    fn unique_reduce(&self, basis: Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&basis));
        todo!()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IntegralSubmoduleOfFinitelyFreeModuleStructure {
    module: FinitelyFreeModuleStructure<IntegerCanonicalStructure>,
}

impl Signature for IntegralSubmoduleOfFinitelyFreeModuleStructure {}

impl SetSignature for IntegralSubmoduleOfFinitelyFreeModuleStructure {
    type Set = SubmoduleOfFinitelyFreeModule<IntegerCanonicalStructure>;

    fn is_element(&self, x: &Self::Set) -> bool {
        self.module.rank() == x.basis.cols()
    }
}

impl SubmoduleOfFinitelyFreeModuleSignature<IntegerCanonicalStructure>
    for IntegralSubmoduleOfFinitelyFreeModuleStructure
{
    fn module(&self) -> &FinitelyFreeModuleStructure<IntegerCanonicalStructure> {
        &self.module
    }
}

impl SubmoduleOfFinitelyFreeModuleUniqueBasisSignature<IntegerCanonicalStructure>
    for IntegralSubmoduleOfFinitelyFreeModuleStructure
{
    fn unique_reduce(&self, basis: Self::Set) -> Self::Set {
        debug_assert!(self.is_element(&basis));
        todo!()
    }
}

// #[derive(Debug, Clone)]
// pub struct FinitelyFreeSubmoduleStructure<
//     Ring: EqSignature + BezoutDomainSignature,
//     Module: FinitelyFreeModuleSignature<Ring>,
// > {
//     ring: PhantomData<Ring>,
//     module: Module,
//     // linearly independent rows
//     basis_matrix: Matrix<Ring::Set>,
// }

// impl<Ring: EqSignature + BezoutDomainSignature, Module: FinitelyFreeModuleSignature<Ring>>
//     FinitelyFreeSubmoduleStructure<Ring, Module>
// {
//     pub fn ring(&self) -> &Ring {
//         self.module.ring()
//     }

//     pub fn module(&self) -> &Module {
//         &self.module
//     }

//     pub fn from_span(module: Module, span: Vec<Module::Set>) -> Self {
//         for v in span {
//             debug_assert!(module.is_element(&v));
//         }
//         module
//         todo!()
//     }

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
// }

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
