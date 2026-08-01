use super::finitely_free_submodules::FinitelyFreeSubmodule;
use crate::{
    linear::finitely_free_module::RingToFinitelyFreeModuleSignature,
    matrix::{ReducedHermiteAlgorithmSignature, RingMatricesSignature},
    structure::*,
};
use algebraeon_sets::sets::EnumeratedFiniteSetStructure;
use algebraeon_structures::*;
use std::{cmp::Ordering, marker::PhantomData};

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmoduleStructure<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> {
    _phantom: PhantomData<(Set, Ring, Module)>,
    module: ModuleB,
    submodule: FinitelyFreeSubmodule<Ring::Elem>,
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> From<FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>>
    for FinitelyFreeSubmodule<Ring::Elem>
{
    fn from(val: FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>) -> Self {
        val.submodule
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> PartialEq for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn eq(&self, other: &Self) -> bool {
        let module = self.module.borrow();
        self.module == other.module && module.submodules().equal(&self.submodule, &other.submodule)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> Eq for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    pub fn new(module: ModuleB, submodule: FinitelyFreeSubmodule<Ring::Elem>) -> Self {
        Self {
            _phantom: PhantomData,
            module,
            submodule,
        }
    }

    pub fn module(&self) -> &Module {
        self.module.borrow()
    }

    pub fn submodule(&self) -> &FinitelyFreeSubmodule<Ring::Elem> {
        &self.submodule
    }

    pub fn rank(&self) -> usize {
        self.submodule.rank()
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> Signature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> SetSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    type Elem = Module::Elem;

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
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring> + EqSignature,
    ModuleB: BorrowedStructure<Module>,
> EqSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.module().equal(a, b)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring> + OrdSignature,
    ModuleB: BorrowedStructure<Module>,
> PartialOrdSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        self.module().partial_cmp(a, b)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring> + OrdSignature,
    ModuleB: BorrowedStructure<Module>,
> OrdSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        self.module().cmp(a, b)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature + OrderedFiniteSetSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> CountableSetSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn into_generate_all_elements(self) -> impl Iterator<Item = Self::Elem> {
        self.list_all_elements().into_iter()
    }

    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        self.list_all_elements().into_iter()
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature + OrderedFiniteSetSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> FiniteSetSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn list_all_elements(&self) -> Vec<Self::Elem> {
        let row_basis = self.submodule.row_basis_matrix();
        self.ring()
            .free_module(EnumeratedFiniteSetStructure::new(row_basis.rows()))
            .generate_all_elements()
            .map(|coeffs| {
                self.module()
                    .from_vec(self.ring().matrix_structure().apply_row(row_basis, &coeffs))
            })
            .collect()
    }

    fn size(&self) -> Natural {
        self.ring().size().pow(&self.rank().into())
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature + OrderedFiniteSetSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring> + OrdSignature,
    ModuleB: BorrowedStructure<Module>,
> OrderedFiniteSetSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
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
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> RinglikeSpecializationSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> ZeroSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn zero(&self) -> Self::Elem {
        self.module().zero()
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> AdditionSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.module().add(a, b)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> CancellativeAdditionSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.module().try_sub(a, b)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> TryNegateSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.module().try_neg(a)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> AdditiveMonoidSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> AdditiveGroupSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        self.module().neg(a)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
    ModuleB: BorrowedStructure<Module>,
> SemiModuleSignature<Ring> for FinitelyFreeSubmoduleStructure<Set, Ring, Module, ModuleB>
{
    fn ring(&self) -> &Ring {
        self.module().ring()
    }

    fn scalar_mul(&self, a: &Self::Elem, x: &Ring::Elem) -> Self::Elem {
        self.module().scalar_mul(a, x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        finite_fields::quaternary_field::QuaternaryField,
        linear::finitely_free_module::FinitelyFreeModuleStructure,
    };

    #[test]
    fn enumeration() {
        let module = FinitelyFreeModuleStructure::new(
            EnumeratedFiniteSetStructure::new(4),
            QuaternaryField::structure(),
        );

        algebraeon_structures::assert_enumerated_ord_finite_set!(
            module.submodule_structure(module.generated_submodule(vec![
                &vec![
                    QuaternaryField::Zero,
                    QuaternaryField::Alpha,
                    QuaternaryField::Beta,
                    QuaternaryField::One
                ],
                &vec![
                    QuaternaryField::Alpha,
                    QuaternaryField::Zero,
                    QuaternaryField::Zero,
                    QuaternaryField::Beta
                ]
            ])),
            16
        );
    }
}
