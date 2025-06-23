use std::borrow::{Borrow, Cow};

use crate::{matrix::Matrix, structure::*};
use algebraeon_sets::structure::*;

pub trait SemiModuleSignature<Ring: SemiRingSignature>: AdditiveMonoidSignature {
    fn ring(&self) -> &Ring;
    fn scalar_mul(&self, a: &Self::Set, x: &Ring::Set) -> Self::Set;
}

pub trait MetaSemiModule<Ring: SemiRingSignature>: MetaType
where
    Self::Signature: SemiModuleSignature<Ring>,
{
    fn scalar_mul(&self, x: &Ring::Set) -> Self {
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

pub trait FinitelyGeneratedModuleSignature<Ring: RingSignature>: ModuleSignature<Ring> {}

pub trait FreeModuleSignature<Basis: SetSignature, Ring: RingSignature>:
    ModuleSignature<Ring>
{
    fn basis_set(&self) -> impl Borrow<Basis>;

    fn to_component<'a>(&self, b: &Basis::Set, v: &'a Self::Set) -> Cow<'a, Ring::Set>;

    fn from_component(&self, b: &Basis::Set, r: &Ring::Set) -> Self::Set;
}

pub trait FinitelyFreeModuleSignature<Basis: FiniteSetSignature, Ring: RingSignature>:
    FreeModuleSignature<Basis, Ring> + FinitelyGeneratedModuleSignature<Ring>
{
    fn basis(&self) -> Vec<Basis::Set> {
        self.basis_set().borrow().list_all_elements()
    }

    fn rank(&self) -> usize {
        self.basis_set().borrow().size()
    }

    fn basis_vecs(&self) -> Vec<Self::Set> {
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

    fn to_vec(&self, a: &Self::Set) -> Vec<Ring::Set> {
        self.basis()
            .iter()
            .map(|b| self.to_component(b, a).as_ref().clone())
            .collect()
    }

    fn from_vec(&self, v: Vec<impl Borrow<Ring::Set>>) -> Self::Set {
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

    fn to_col(&self, a: &Self::Set) -> Matrix<Ring::Set> {
        let basis = self.basis();
        Matrix::construct(self.rank(), 1, |r, _c| {
            self.to_component(&basis[r], a).into_owned()
        })
    }

    fn from_col(&self, v: Matrix<Ring::Set>) -> Self::Set {
        assert_eq!(v.cols(), 1);
        assert_eq!(v.rows(), self.rank());
        self.from_vec((0..self.rank()).map(|r| v.at(r, 0).unwrap()).collect())
    }

    fn to_row(&self, a: &Self::Set) -> Matrix<Ring::Set> {
        self.to_col(a).transpose()
    }

    fn from_row(&self, v: Matrix<Ring::Set>) -> Self::Set {
        self.from_col(v.transpose())
    }
}

pub trait LinearTransformation<
    Ring: RingSignature,
    Domain: ModuleSignature<Ring>,
    Range: ModuleSignature<Ring>,
>: Function<Domain, Range>
{
}
