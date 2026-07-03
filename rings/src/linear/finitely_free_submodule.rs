use super::finitely_free_submodules::FinitelyFreeSubmodule;
use crate::{
    linear::finitely_free_module::FinitelyFreeModuleStructure,
    matrix::UniqueReducedHermiteAlgorithmSignature, structure::*,
};
use algebraeon_structures::*;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmoduleStructure<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> {
    _phantom: PhantomData<(Set, SetB, Ring, RingB)>,
    module: ModuleB,
    submodule: FinitelyFreeSubmodule<Ring::Elem>,
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> From<FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>>
    for FinitelyFreeSubmodule<Ring::Elem>
{
    fn from(val: FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>) -> Self {
        val.submodule
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> PartialEq for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn eq(&self, other: &Self) -> bool {
        let module = self.module.borrow();
        self.module == other.module && module.submodules().equal(&self.submodule, &other.submodule)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> Eq for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: RingSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    pub fn new(module: ModuleB, submodule: FinitelyFreeSubmodule<Ring::Elem>) -> Self {
        Self {
            _phantom: PhantomData,
            module,
            submodule,
        }
    }

    pub fn module(&self) -> &FinitelyFreeModuleStructure<Set, SetB, Ring, RingB> {
        self.module.borrow()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> Signature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> SetSignature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    type Elem = Vec<Ring::Elem>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        self.module().validate_element(x)?;
        if !self
            .module()
            .submodules()
            .contains_element(&self.submodule, x)
        {
            return Err("submodule does not contain element".to_string());
        }
        Ok(())
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> RinglikeSpecializationSignature
    for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> ZeroSignature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn zero(&self) -> Self::Elem {
        self.module().zero()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> AdditionSignature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.module().add(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> CancellativeAdditionSignature
    for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.module().try_sub(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> TryNegateSignature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.module().try_neg(a)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> AdditiveMonoidSignature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> AdditiveGroupSignature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        self.module().neg(a)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> SemiModuleSignature<Ring> for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn ring(&self) -> &Ring {
        self.module().ring()
    }

    fn scalar_mul(&self, a: &Self::Elem, x: &Ring::Elem) -> Self::Elem {
        self.module().scalar_mul(a, x)
    }
}
