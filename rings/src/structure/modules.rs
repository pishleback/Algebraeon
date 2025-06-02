use std::borrow::Borrow;

use crate::structure::*;
use algebraeon_sets::structure::*;

pub trait SemiModuleSignature<Ring: SemiRingSignature>: AdditiveMonoidSignature {
    fn ring(&self) -> &Ring;
    fn scalar_mul(&self, x: &Ring::Set, a: &Self::Set) -> Self::Set;
}

pub trait ModuleSignature<Ring: RingSignature>:
    SemiModuleSignature<Ring> + AdditiveGroupSignature
{
}
impl<Ring: RingSignature, Module: SemiModuleSignature<Ring> + AdditiveGroupSignature>
    ModuleSignature<Ring> for Module
{
}

pub trait FreeModuleSignature<Ring: RingSignature>: ModuleSignature<Ring> {
    type Basis: SetSignature;

    fn basis_set(&self) -> impl Borrow<Self::Basis>;

    fn to_component<'a>(
        &'a self,
        b: &<Self::Basis as SetSignature>::Set,
        v: &'a Self::Set,
    ) -> &'a Ring::Set;

    fn from_component(&self, b: &<Self::Basis as SetSignature>::Set, r: &Ring::Set) -> Self::Set;
}

pub trait FinitelyFreeModuleSignature<Ring: RingSignature>: FreeModuleSignature<Ring>
where
    Self::Basis: FiniteSetSignature,
{
    fn basis(&self) -> Vec<<Self::Basis as SetSignature>::Set> {
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

    fn to_vec(&self, v: &Self::Set) -> Vec<Ring::Set> {
        self.basis()
            .iter()
            .map(|b| self.to_component(b, v).clone())
            .collect()
    }

    fn from_vec(&self, v: Vec<&Ring::Set>) -> Self::Set {
        let n = self.rank();
        debug_assert_eq!(v.len(), n);
        let basis = self.basis();
        debug_assert_eq!(basis.len(), n);
        let mut t = self.zero();
        for i in 0..n {
            self.add_mut(
                &mut t,
                &self.scalar_mul(v[i], &self.from_component(&basis[i], &self.ring().one())),
            );
        }
        t
    }
}

pub trait LinearTransformation<
    Ring: RingSignature,
    Domain: ModuleSignature<Ring>,
    Range: ModuleSignature<Ring>,
>: Function<Domain, Range>
{
}

// pub trait SubModuleSignature<Ring: RingSignature, Module: ModuleSignature<Ring>>:
//     EqSignature
// {
//     fn ring(&self) -> &Ring;
//     fn module(&self) -> &Module;
//     fn zero_submodule(&self) -> Self::Set {
//         self.generated(vec![])
//     }
//     fn improper_submodule(&self) -> Self::Set;
//     fn add(&self, x: &Self::Set, y: &Self::Set) -> Self::Set;
//     fn intersect(&self, x: &Self::Set, y: &Self::Set) -> Self::Set;
//     fn generated(&self, generators: Vec<&Module::Set>) -> Self::Set;
//     /// Does x contain p
//     fn contains_element(&self, x: &Self::Set, p: &Module::Set) -> bool;
//     /// Does x contain y
//     fn contains(&self, x: &Self::Set, y: &Self::Set) -> bool;
// }

// pub trait SubModuleCosetSignature<Ring: RingSignature, Module: ModuleSignature<Ring>>:
//     EqSignature
// {
//     fn ring(&self) -> &Ring;
//     fn module(&self) -> &Module;
//     fn add(&self, x: &Self::Set, y: &Self::Set) -> Self::Set;
//     fn intersect(&self, x: &Self::Set, y: &Self::Set) -> Self::Set;
// }
