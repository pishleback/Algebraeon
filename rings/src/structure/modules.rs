use crate::{
    linear::{
        finitely_free_affine_subsets::FinitelyFreeSubmoduleAffineSubsetsStructure,
        finitely_free_cosets::FinitelyFreeSubmoduleCosetsStructure,
        finitely_free_submodule::FinitelyFreeSubmoduleStructure,
        finitely_free_submodules::{FinitelyFreeSubmodule, FinitelyFreeSubmodulesStructure},
    },
    matrix::{Matrix, MatrixStructure, ReducedHermiteAlgorithmSignature},
    structure::*,
};
use algebraeon_structures::*;
use std::borrow::{Borrow, Cow};

pub trait SemiModuleSignature<Ring: SemiRingSignature>: AdditiveMonoidSignature {
    fn ring(&self) -> &Ring;
    fn scalar_mul(&self, a: &Self::Elem, x: &Ring::Elem) -> Self::Elem;
}

pub trait MetaSemiModule<Ring: SemiRingSignature>: MetaType
where
    Self::Signature: SemiModuleSignature<Ring>,
{
    fn scalar_mul(&self, x: &Ring::Elem) -> Self {
        Self::structure().scalar_mul(self, x)
    }
}
impl<Ring: SemiRingSignature, T: MetaType> MetaSemiModule<Ring> for T where
    Self::Signature: SemiModuleSignature<Ring>
{
}

pub trait ModuleSignature<Ring: RingSignature>:
    SemiModuleSignature<Ring> + AdditiveGroupSignature
{
}
impl<Ring: RingSignature, Module: SemiModuleSignature<Ring> + AdditiveGroupSignature>
    ModuleSignature<Ring> for Module
{
}

pub trait FinitelyGeneratedModuleSignature<GeneratingSet: FiniteSetSignature, Ring: RingSignature>:
    ModuleSignature<Ring>
{
}

pub trait FreeModuleSignature<Basis: SetSignature, Ring: RingSignature>:
    ModuleSignature<Ring>
{
    fn basis_set(&self) -> impl Borrow<Basis>;

    fn to_component<'a>(&self, b: &Basis::Elem, v: &'a Self::Elem) -> Cow<'a, Ring::Elem>;

    fn from_component(&self, b: &Basis::Elem, r: &Ring::Elem) -> Self::Elem;
}

pub trait FinitelyFreeModuleSignature<Basis: FiniteSetSignature, Ring: RingSignature>:
    FreeModuleSignature<Basis, Ring> + FinitelyGeneratedModuleSignature<Basis, Ring>
{
    fn basis(&self) -> Vec<Basis::Elem> {
        self.basis_set().borrow().list_all_elements()
    }

    fn rank(&self) -> usize {
        self.basis_set()
            .borrow()
            .size()
            .try_into()
            .expect("too large")
    }

    /// The elementary basis vectors
    fn basis_vecs(&self) -> Vec<Self::Elem> {
        let zero = self.ring().zero();
        let one = self.ring().one();
        (0..self.rank())
            .map(|j| {
                self.from_vec(
                    (0..self.rank())
                        .map(|i| if i == j { &one } else { &zero })
                        .collect(),
                )
            })
            .collect()
    }

    fn to_vec(&self, a: &Self::Elem) -> Vec<Ring::Elem> {
        self.basis()
            .iter()
            .map(|b| self.to_component(b, a).as_ref().clone())
            .collect()
    }

    fn from_vec(&self, v: Vec<impl Borrow<Ring::Elem>>) -> Self::Elem {
        let n = self.rank();
        debug_assert_eq!(v.len(), n);
        let basis = self.basis();
        debug_assert_eq!(basis.len(), n);
        let mut t = self.zero();
        for i in 0..n {
            self.add_mut(
                &mut t,
                &self.scalar_mul(
                    &self.from_component(&basis[i], &self.ring().one()),
                    v[i].borrow(),
                ),
            );
        }
        t
    }

    fn to_col(&self, a: &Self::Elem) -> Matrix<Ring::Elem> {
        let basis = self.basis();
        Matrix::construct(self.rank(), 1, |r, _c| {
            self.to_component(&basis[r], a).into_owned()
        })
    }

    fn from_col(&self, v: Matrix<Ring::Elem>) -> Self::Elem {
        assert_eq!(v.cols(), 1);
        assert_eq!(v.rows(), self.rank());
        self.from_vec((0..self.rank()).map(|r| v.at(r, 0).unwrap()).collect())
    }

    fn to_row(&self, a: &Self::Elem) -> Matrix<Ring::Elem> {
        self.to_col(a).transpose()
    }

    fn from_row(&self, v: Matrix<Ring::Elem>) -> Self::Elem {
        self.from_col(v.transpose())
    }

    fn submodule_structure(
        &self,
        submodule: FinitelyFreeSubmodule<Ring::Elem>,
    ) -> FinitelyFreeSubmoduleStructure<Basis, Ring, Self, &Self>
    where
        Basis: EnumeratedOrdFiniteSetSignature,
        Ring: ReducedHermiteAlgorithmSignature,
    {
        FinitelyFreeSubmoduleStructure::new(self, submodule)
    }

    fn submodules(&self) -> FinitelyFreeSubmodulesStructure<Basis, Ring, Self, &Self>
    where
        Basis: EnumeratedOrdFiniteSetSignature,
        Ring: ReducedHermiteAlgorithmSignature,
    {
        FinitelyFreeSubmodulesStructure::new(self)
    }

    fn into_submodules(self) -> FinitelyFreeSubmodulesStructure<Basis, Ring, Self, Self>
    where
        Basis: EnumeratedOrdFiniteSetSignature,
        Ring: ReducedHermiteAlgorithmSignature,
    {
        FinitelyFreeSubmodulesStructure::new(self)
    }

    fn cosets(&self) -> FinitelyFreeSubmoduleCosetsStructure<Basis, Ring, Self, &Self>
    where
        Basis: EnumeratedOrdFiniteSetSignature,
        Ring: ReducedHermiteAlgorithmSignature,
    {
        FinitelyFreeSubmoduleCosetsStructure::new(self)
    }

    fn into_cosets(self) -> FinitelyFreeSubmoduleCosetsStructure<Basis, Ring, Self, Self>
    where
        Basis: EnumeratedOrdFiniteSetSignature,
        Ring: ReducedHermiteAlgorithmSignature,
    {
        FinitelyFreeSubmoduleCosetsStructure::new(self)
    }

    fn affine_subsets(
        &self,
    ) -> FinitelyFreeSubmoduleAffineSubsetsStructure<Basis, Ring, Self, &Self>
    where
        Basis: EnumeratedOrdFiniteSetSignature,
        Ring: ReducedHermiteAlgorithmSignature,
    {
        FinitelyFreeSubmoduleAffineSubsetsStructure::new(self)
    }

    fn into_affine_subsets(
        self,
    ) -> FinitelyFreeSubmoduleAffineSubsetsStructure<Basis, Ring, Self, Self>
    where
        Basis: EnumeratedOrdFiniteSetSignature,
        Ring: ReducedHermiteAlgorithmSignature,
    {
        FinitelyFreeSubmoduleAffineSubsetsStructure::new(self)
    }

    fn improper_submodule(&self) -> FinitelyFreeSubmodule<Ring::Elem>
    where
        Basis: EnumeratedOrdFiniteSetSignature,
        Ring: ReducedHermiteAlgorithmSignature,
    {
        self.submodules()
            .matrix_row_span(MatrixStructure::new(self.ring().clone()).ident(self.rank()))
    }

    fn generated_submodule(&self, generators: Vec<&Self::Elem>) -> FinitelyFreeSubmodule<Ring::Elem>
    where
        Basis: EnumeratedOrdFiniteSetSignature,
        Ring: ReducedHermiteAlgorithmSignature,
    {
        for generator in &generators {
            debug_assert!(self.validate_element(generator).is_ok());
        }
        let generators = generators
            .into_iter()
            .map(|g| self.to_vec(g))
            .collect::<Vec<_>>();
        let row_span = Matrix::construct(generators.len(), self.rank(), |r, c| {
            generators[r][c].clone()
        });
        self.submodules().matrix_row_span(row_span)
    }
}

impl<Basis: FiniteSetSignature, Ring: RingSignature, Module: FreeModuleSignature<Basis, Ring>>
    FinitelyGeneratedModuleSignature<Basis, Ring> for Module
{
}

impl<Basis: FiniteSetSignature, Ring: RingSignature, Module: FreeModuleSignature<Basis, Ring>>
    FinitelyFreeModuleSignature<Basis, Ring> for Module
{
}

pub trait LinearTransformation<
    Ring: RingSignature,
    Domain: ModuleSignature<Ring>,
    Range: ModuleSignature<Ring>,
>: FunctionMorphism<Domain, Range>
{
}
