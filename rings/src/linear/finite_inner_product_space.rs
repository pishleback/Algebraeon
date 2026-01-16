use crate::{
    linear::finitely_free_module::{
        FinitelyFreeModuleStructure, RingToFinitelyFreeModuleSignature,
    },
    matrix::{SymmetricMatrix, ToSymmetrixMatricesSignature},
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, FreeModuleSignature, RealSubsetSignature, RingSignature,
        RinglikeSpecializationSignature, SemiModuleSignature, TryNegateSignature, ZeroSignature,
    },
};
use algebraeon_sets::structure::{
    BorrowedStructure, EnumeratedFiniteSetStructure, EqSignature, SetSignature, Signature,
};
use std::{
    borrow::{Borrow, Cow},
    marker::PhantomData,
};

#[derive(Debug, Clone)]
pub struct FiniteInnerProductSpaceStructure<
    Set: EqSignature + RingSignature + RealSubsetSignature,
    SetB: BorrowedStructure<Set>,
> {
    _set: PhantomData<Set>,
    set: SetB,
    basis_set: EnumeratedFiniteSetStructure,
    inner_product_matrix: SymmetricMatrix<Set::Set>,
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    FiniteInnerProductSpaceStructure<Set, SetB>
{
    pub fn set(&self) -> &Set {
        self.set.borrow()
    }

    pub fn n(&self) -> usize {
        self.inner_product_matrix.n()
    }

    pub fn coord_set_free_module_restructure(&self) -> FinitelyFreeModuleStructure<Set, &Set> {
        self.set().free_module(self.n())
    }

    pub fn into_coord_set_free_module_restructure(self) -> FinitelyFreeModuleStructure<Set, SetB> {
        let n = self.n();
        FinitelyFreeModuleStructure::new(self.set, n)
    }
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>> PartialEq
    for FiniteInnerProductSpaceStructure<Set, SetB>
{
    fn eq(&self, other: &Self) -> bool {
        let scalar_matrix_structure = self.set().symmetric_matrix_structure();
        debug_assert_eq!(
            scalar_matrix_structure,
            other.set().symmetric_matrix_structure()
        );

        self.set == other.set
            && self.basis_set == other.basis_set
            && scalar_matrix_structure
                .equal(&self.inner_product_matrix, &other.inner_product_matrix)
    }
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>> Eq
    for FiniteInnerProductSpaceStructure<Set, SetB>
{
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>> Signature
    for FiniteInnerProductSpaceStructure<Set, SetB>
{
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    SetSignature for FiniteInnerProductSpaceStructure<Set, SetB>
{
    type Set = Vec<Set::Set>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        self.coord_set_free_module_restructure().is_element(x)
    }
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    EqSignature for FiniteInnerProductSpaceStructure<Set, SetB>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.coord_set_free_module_restructure().equal(a, b)
    }
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    RinglikeSpecializationSignature for FiniteInnerProductSpaceStructure<Set, SetB>
{
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    ZeroSignature for FiniteInnerProductSpaceStructure<Set, SetB>
{
    fn zero(&self) -> Self::Set {
        self.coord_set_free_module_restructure().zero()
    }
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    AdditionSignature for FiniteInnerProductSpaceStructure<Set, SetB>
{
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.coord_set_free_module_restructure().add(a, b)
    }
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    CancellativeAdditionSignature for FiniteInnerProductSpaceStructure<Set, SetB>
{
    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        self.coord_set_free_module_restructure().try_sub(a, b)
    }
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    TryNegateSignature for FiniteInnerProductSpaceStructure<Set, SetB>
{
    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set> {
        self.coord_set_free_module_restructure().try_neg(a)
    }
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    AdditiveMonoidSignature for FiniteInnerProductSpaceStructure<Set, SetB>
{
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    AdditiveGroupSignature for FiniteInnerProductSpaceStructure<Set, SetB>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.coord_set_free_module_restructure().neg(a)
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.coord_set_free_module_restructure().sub(a, b)
    }
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    SemiModuleSignature<Set> for FiniteInnerProductSpaceStructure<Set, SetB>
{
    fn ring(&self) -> &Set {
        self.set()
    }

    fn scalar_mul(&self, v: &Self::Set, r: &Set::Set) -> Self::Set {
        self.coord_set_free_module_restructure().scalar_mul(v, r)
    }
}

impl<Set: EqSignature + RingSignature + RealSubsetSignature, SetB: BorrowedStructure<Set>>
    FreeModuleSignature<Set> for FiniteInnerProductSpaceStructure<Set, SetB>
{
    type Basis = EnumeratedFiniteSetStructure;

    fn basis_set(&self) -> impl Borrow<Self::Basis> {
        EnumeratedFiniteSetStructure::new(self.n())
    }

    fn to_component<'a>(&self, b: &usize, v: &'a Self::Set) -> Cow<'a, Set::Set> {
        self.coord_set_free_module_restructure().to_component(b, v)
    }

    fn from_component(&self, b: &usize, r: &<Set>::Set) -> Self::Set {
        self.coord_set_free_module_restructure()
            .from_component(b, r)
    }
}
