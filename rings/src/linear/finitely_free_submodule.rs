use super::finitely_free_submodules::FinitelyFreeSubmodule;
use crate::{
    linear::finitely_free_module::{
        FinitelyFreeModuleStructure, RingToFinitelyFreeModuleSignature,
    },
    matrix::{RingMatricesSignature, UniqueReducedHermiteAlgorithmSignature},
    structure::*,
};
use algebraeon_sets::sets::EnumeratedFiniteSetStructure;
use algebraeon_structures::*;
use std::{cmp::Ordering, marker::PhantomData};

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

    pub fn set(&self) -> &Set {
        self.module().set()
    }

    pub fn rank(&self) -> usize {
        self.submodule.rank()
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
> EqSignature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.module().equal(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature + OrdSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> PartialOrdSignature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        self.module().partial_cmp(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature + OrdSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> OrdSignature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        self.module().cmp(a, b)
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature + EnumeratedOrdFiniteSetSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> CountableSetSignature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.list_all_elements().into_iter()
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.list_all_elements().into_iter()
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature + EnumeratedOrdFiniteSetSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> FiniteSetSignature for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn list_all_elements(&self) -> Vec<Self::Elem> {
        let row_basis = self.submodule.row_basis_matrix();
        self.ring()
            .free_module(EnumeratedFiniteSetStructure::new(row_basis.rows()))
            .generate_all_elements()
            .map(|coeffs| self.ring().matrix_structure().apply_row(row_basis, &coeffs))
            .collect()
    }

    fn size(&self) -> Natural {
        self.ring().size().pow(&self.rank().into())
    }
}

impl<
    Set: EnumeratedOrdFiniteSetSignature,
    SetB: BorrowedStructure<Set>,
    Ring: UniqueReducedHermiteAlgorithmSignature + EnumeratedOrdFiniteSetSignature,
    RingB: BorrowedStructure<Ring>,
    ModuleB: BorrowedStructure<FinitelyFreeModuleStructure<Set, SetB, Ring, RingB>>,
> EnumeratedOrdFiniteSetSignature
    for FinitelyFreeSubmoduleStructure<Set, SetB, Ring, RingB, ModuleB>
{
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.sort(self.list_all_elements())
    }

    fn element_to_enumeration(&self, elem: &Self::Elem) -> Natural {
        debug_assert!(self.is_element(elem));
        self.binary_search_index(&self.list_all_elements_ordered(), elem)
            .unwrap()
            .into()
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        if let Ok(num) = TryInto::<usize>::try_into(num) {
            self.list_all_elements_ordered().into_iter().nth(num)
        } else {
            None
        }
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
