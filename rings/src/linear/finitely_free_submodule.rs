use super::finitely_free_submodules::FinitelyFreeSubmodule;
use crate::{
    linear::finitely_free_module::RingToFinitelyFreeModuleSignature,
    matrix::{ReducedHermiteAlgorithmSignature, RingMatricesSignature},
    structure::*,
};
use algebraeon_sets::sets::EnumeratedFiniteSetStructure;
use algebraeon_structures::*;
use std::{cmp::Ordering, marker::PhantomData, sync::Arc};

#[derive(Debug, Clone)]
pub struct FinitelyFreeSubmoduleStructure<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> {
    _phantom: PhantomData<(Set, Ring, Module)>,
    module: Arc<Module>,
    submodule: Arc<FinitelyFreeSubmodule<Ring::Elem>>,
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> From<FinitelyFreeSubmoduleStructure<Set, Ring, Module>>
    for Arc<FinitelyFreeSubmodule<Ring::Elem>>
{
    fn from(val: FinitelyFreeSubmoduleStructure<Set, Ring, Module>) -> Self {
        val.submodule
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> PartialEq for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
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
> Eq for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    pub fn new(
        module: Arc<Module>,
        submodule: Arc<FinitelyFreeSubmodule<Ring::Elem>>,
    ) -> Arc<Self> {
        Self {
            _phantom: PhantomData,
            module,
            submodule,
        }
        .into()
    }

    pub fn module(&self) -> &Arc<Module> {
        &self.module
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
> Signature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> SetSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
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
> EqSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn equal(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.module().equal(a, b)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring> + OrdSignature,
> PartialOrdSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn partial_cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        self.module().partial_cmp(a, b)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring> + OrdSignature,
> OrdSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn cmp(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        self.module().cmp(a, b)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature + OrderedFiniteSetSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> CountableSetSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn generate_all_elements(self: Arc<Self>) -> impl Iterator<Item = Self::Elem> {
        self.list_all_elements().into_iter()
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature + OrderedFiniteSetSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> FiniteSetSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn list_all_elements(self: &Arc<Self>) -> Vec<Self::Elem> {
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

    fn size(self: &Arc<Self>) -> Natural {
        self.ring().size().pow(&self.rank().into())
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature + OrderedFiniteSetSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring> + OrdSignature,
> OrderedFiniteSetSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn list_all_elements_ordered(self: &Arc<Self>) -> Vec<Self::Elem> {
        self.sort(self.list_all_elements())
    }

    fn element_to_enumeration(self: &Arc<Self>, elem: &Self::Elem) -> Natural {
        debug_assert!(self.is_element(elem));
        self.binary_search_index(&self.list_all_elements_ordered(), elem)
            .unwrap()
            .into()
    }

    fn enumeration_to_element(self: &Arc<Self>, num: &Natural) -> Option<Self::Elem> {
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
> RinglikeSpecializationSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> ZeroSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn zero(self: &Arc<Self>) -> Self::Elem {
        self.module().zero()
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> AdditionSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn add(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.module().add(a, b)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> CancellativeAdditionSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn try_sub(self: &Arc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.module().try_sub(a, b)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> TryNegateSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn try_neg(self: &Arc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        self.module().try_neg(a)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> AdditiveMonoidSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> AdditiveGroupSignature for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn neg(self: &Arc<Self>, a: &Self::Elem) -> Self::Elem {
        self.module().neg(a)
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> SemiModuleSignature<Ring> for FinitelyFreeSubmoduleStructure<Set, Ring, Module>
{
    fn ring(self: &Arc<Self>) -> Arc<Ring> {
        self.module().ring()
    }

    fn scalar_mul(self: &Arc<Self>, a: &Self::Elem, x: &Ring::Elem) -> Self::Elem {
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
