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
use std::{borrow::Cow, marker::PhantomData};

#[derive(Debug, Clone)]
pub struct FiniteInnerProductSpaceStructure<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> {
    _coords: PhantomData<Coords>,
    coords: CoordsB,
    _scalars: PhantomData<Scalars>,
    scalars: ScalarsB,
    basis_set: EnumeratedFiniteSetStructure,
    inner_product_matrix: SymmetricMatrix<Scalars::Set>,
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
    pub fn coord_set(&self) -> &Coords {
        self.coords.borrow()
    }

    pub fn scalar_set(&self) -> &Scalars {
        self.scalars.borrow()
    }

    pub fn n(&self) -> usize {
        self.inner_product_matrix.n()
    }

    pub fn coord_set_free_module_restructure(
        &self,
    ) -> FinitelyFreeModuleStructure<Coords, &Coords> {
        self.coord_set().free_module(self.n())
    }

    pub fn into_coord_set_free_module_restructure(
        self,
    ) -> FinitelyFreeModuleStructure<Coords, CoordsB> {
        let n = self.n();
        FinitelyFreeModuleStructure::new(self.coords, n)
    }
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> PartialEq for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
    fn eq(&self, other: &Self) -> bool {
        let scalar_matrix_structure = self.scalar_set().symmetric_matrix_structure();
        debug_assert_eq!(
            scalar_matrix_structure,
            other.scalar_set().symmetric_matrix_structure()
        );

        self.coords == other.coords
            && self.scalars == other.scalars
            && self.basis_set == other.basis_set
            && scalar_matrix_structure
                .equal(&self.inner_product_matrix, &other.inner_product_matrix)
    }
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> Eq for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> Signature for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> SetSignature for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
    type Set = Vec<Coords::Set>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        self.coord_set_free_module_restructure().is_element(x)
    }
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> EqSignature for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.coord_set_free_module_restructure().equal(a, b)
    }
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> RinglikeSpecializationSignature
    for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> ZeroSignature for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
    fn zero(&self) -> Self::Set {
        self.coord_set_free_module_restructure().zero()
    }
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> AdditionSignature for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.coord_set_free_module_restructure().add(a, b)
    }
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> CancellativeAdditionSignature
    for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        self.coord_set_free_module_restructure().try_sub(a, b)
    }
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> TryNegateSignature for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set> {
        self.coord_set_free_module_restructure().try_neg(a)
    }
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> AdditiveMonoidSignature for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> AdditiveGroupSignature for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.coord_set_free_module_restructure().neg(a)
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.coord_set_free_module_restructure().sub(a, b)
    }
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> SemiModuleSignature<Coords>
    for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
    fn ring(&self) -> &Coords {
        self.coord_set()
    }

    fn scalar_mul(&self, v: &Self::Set, r: &Coords::Set) -> Self::Set {
        self.coord_set_free_module_restructure().scalar_mul(v, r)
    }
}

impl<
    Coords: EqSignature + RingSignature + RealSubsetSignature,
    CoordsB: BorrowedStructure<Coords>,
    Scalars: EqSignature + RingSignature + RealSubsetSignature,
    ScalarsB: BorrowedStructure<Scalars>,
> FreeModuleSignature<Coords>
    for FiniteInnerProductSpaceStructure<Coords, CoordsB, Scalars, ScalarsB>
{
    type Basis = EnumeratedFiniteSetStructure;

    fn basis_set<'a>(&'a self) -> impl std::borrow::Borrow<Self::Basis> + 'a {
        let z = self.coord_set_free_module_restructure();
        z.basis_set()
    }

    fn to_component<'a>(&self, b: &usize, v: &'a Self::Set) -> Cow<'a, Coords::Set> {
        self.coord_set_free_module_restructure().to_component(b, v)
    }

    fn from_component(&self, b: &usize, r: &<Coords>::Set) -> Self::Set {
        self.coord_set_free_module_restructure()
            .from_component(b, r)
    }
}
