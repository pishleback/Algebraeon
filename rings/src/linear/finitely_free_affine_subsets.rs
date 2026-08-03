use super::finitely_free_cosets::FinitelyFreeSubmoduleCoset;
use crate::matrix::ReducedHermiteAlgorithmSignature;
use crate::structure::*;
use algebraeon_structures::*;
use std::borrow::Borrow;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::rc::Rc;

#[derive(Debug, Clone)]
pub enum FinitelyFreeSubmoduleAffineSubset<Set: Clone + Debug> {
    Empty,
    NonEmpty(FinitelyFreeSubmoduleCoset<Set>),
}

impl<Set: Clone + Debug> From<FinitelyFreeSubmoduleCoset<Set>>
    for FinitelyFreeSubmoduleAffineSubset<Set>
{
    fn from(submodule: FinitelyFreeSubmoduleCoset<Set>) -> Self {
        Self::NonEmpty(submodule)
    }
}

impl<Set: Clone + Debug> FinitelyFreeSubmoduleAffineSubset<Set> {
    pub fn rank(&self) -> Option<usize> {
        match self {
            FinitelyFreeSubmoduleAffineSubset::Empty => None,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => Some(coset.rank()),
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty => true,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(..) => false,
        }
    }

    pub fn unwrap_to_coset(self) -> FinitelyFreeSubmoduleCoset<Set> {
        match self {
            FinitelyFreeSubmoduleAffineSubset::Empty => {
                panic!("unwrap called on empty affine subset")
            }
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => coset,
        }
    }

    pub fn affine_rank(&self) -> usize {
        match &self {
            FinitelyFreeSubmoduleAffineSubset::Empty => 0,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => coset.submodule().rank() + 1,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinitelyFreeSubmoduleAffineSubsetsStructure<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> {
    _set: PhantomData<Set>,
    _ring: PhantomData<Ring>,
    module: Rc<Module>,
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> FinitelyFreeSubmoduleAffineSubsetsStructure<Set, Ring, Module>
{
    pub fn new(module: Rc<Module>) -> Rc<Self> {
        Self {
            _set: PhantomData,
            _ring: PhantomData,
            module,
        }
        .into()
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> Signature for FinitelyFreeSubmoduleAffineSubsetsStructure<Set, Ring, Module>
{
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> SetSignature for FinitelyFreeSubmoduleAffineSubsetsStructure<Set, Ring, Module>
{
    type Elem = FinitelyFreeSubmoduleAffineSubset<Ring::Elem>;

    fn validate_element(self: &Rc<Self>, _x: &Self::Elem) -> Result<(), String> {
        //TODO: better checks
        Ok(())
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring>,
> FinitelyFreeSubmoduleAffineSubsetsStructure<Set, Ring, Module>
{
    pub fn ring(&self) -> Rc<Ring> {
        self.module().ring()
    }

    pub fn module(&self) -> &Rc<Module> {
        &self.module
    }

    pub fn from_affine_span(
        &self,
        mut span: Vec<&Module::Elem>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Elem> {
        for v in &span {
            debug_assert!(self.module().is_element(v));
        }
        if let Some(offset) = span.pop() {
            let linear_span = span
                .into_iter()
                .map(|v| self.module().add(v, &self.module().neg(offset)))
                .collect::<Vec<_>>();
            FinitelyFreeSubmoduleAffineSubset::from(
                self.module().cosets().from_offset_and_submodule(
                    offset,
                    self.module()
                        .submodules()
                        .span(linear_span.iter().collect()),
                ),
            )
        } else {
            FinitelyFreeSubmoduleAffineSubset::Empty
        }
    }

    pub fn affine_basis(
        &self,
        subset: &FinitelyFreeSubmoduleAffineSubset<Ring::Elem>,
    ) -> Vec<Module::Elem> {
        match &subset {
            FinitelyFreeSubmoduleAffineSubset::Empty => vec![],
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => {
                let offset = self.module().from_vec(coset.offset().clone());
                let mut affine_basis = vec![offset.clone()];
                for v in coset.submodule().basis() {
                    affine_basis.push(self.module().add(&self.module().from_vec(v), &offset));
                }
                affine_basis
            }
        }
    }

    pub fn add(
        self: &Rc<Self>,
        x: FinitelyFreeSubmoduleAffineSubset<Ring::Elem>,
        y: FinitelyFreeSubmoduleAffineSubset<Ring::Elem>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Elem> {
        debug_assert!(self.validate_element(&x).is_ok());
        debug_assert!(self.validate_element(&y).is_ok());
        match (x, y) {
            (FinitelyFreeSubmoduleAffineSubset::Empty, _)
            | (_, FinitelyFreeSubmoduleAffineSubset::Empty) => {
                FinitelyFreeSubmoduleAffineSubset::Empty
            }
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(x_coset),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(y_coset),
            ) => FinitelyFreeSubmoduleAffineSubset::NonEmpty(
                self.module().cosets().add(x_coset, y_coset),
            ),
        }
    }

    pub fn intersect(
        self: &Rc<Self>,
        x: &FinitelyFreeSubmoduleAffineSubset<Ring::Elem>,
        y: &FinitelyFreeSubmoduleAffineSubset<Ring::Elem>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Elem> {
        debug_assert!(self.validate_element(x).is_ok());
        debug_assert!(self.validate_element(y).is_ok());
        match (x, y) {
            (FinitelyFreeSubmoduleAffineSubset::Empty, _)
            | (_, FinitelyFreeSubmoduleAffineSubset::Empty) => {
                FinitelyFreeSubmoduleAffineSubset::Empty
            }
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(x_coset),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(y_coset),
            ) => self.module().cosets().intersect(x_coset, y_coset),
        }
    }

    pub fn intersect_list(
        self: &Rc<Self>,
        xs: Vec<impl Borrow<FinitelyFreeSubmoduleAffineSubset<Ring::Elem>>>,
    ) -> FinitelyFreeSubmoduleAffineSubset<Ring::Elem> {
        for x in &xs {
            debug_assert!(self.validate_element(x.borrow()).is_ok());
        }

        let mut i = FinitelyFreeSubmoduleAffineSubset::NonEmpty(
            self.module()
                .cosets()
                .from_submodule(self.module().submodules().full_submodule()),
        );
        for x in xs {
            i = self.intersect(&i, x.borrow());
        }
        i
    }

    pub fn contains_element(
        self: &Rc<Self>,
        x: &FinitelyFreeSubmoduleAffineSubset<Ring::Elem>,
        p: &Module::Elem,
    ) -> bool {
        debug_assert!(self.validate_element(x).is_ok());
        debug_assert!(self.module().validate_element(p).is_ok());
        match x {
            FinitelyFreeSubmoduleAffineSubset::Empty => false,
            FinitelyFreeSubmoduleAffineSubset::NonEmpty(coset) => {
                self.module().cosets().contains_element(coset, p)
            }
        }
    }

    pub fn equal_slow(
        self: &Rc<Self>,
        x: &FinitelyFreeSubmoduleAffineSubset<Ring::Elem>,
        y: &FinitelyFreeSubmoduleAffineSubset<Ring::Elem>,
    ) -> bool {
        debug_assert!(self.validate_element(x).is_ok());
        debug_assert!(self.validate_element(y).is_ok());
        #[allow(clippy::match_same_arms)]
        match (x, y) {
            (
                FinitelyFreeSubmoduleAffineSubset::Empty,
                FinitelyFreeSubmoduleAffineSubset::Empty,
            ) => true,
            (
                FinitelyFreeSubmoduleAffineSubset::Empty,
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
                FinitelyFreeSubmoduleAffineSubset::Empty,
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(x_coset),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(y_coset),
            ) => self.module().cosets().equal_slow(x_coset, y_coset),
        }
    }
}

impl<
    Set: OrderedFiniteSetSignature,
    Ring: ReducedHermiteAlgorithmSignature,
    Module: FinitelyFreeModuleSignature<Set, Ring> + EqSignature,
> EqSignature for FinitelyFreeSubmoduleAffineSubsetsStructure<Set, Ring, Module>
{
    fn equal(
        self: &Rc<Self>,
        x: &FinitelyFreeSubmoduleAffineSubset<Ring::Elem>,
        y: &FinitelyFreeSubmoduleAffineSubset<Ring::Elem>,
    ) -> bool {
        debug_assert!(self.validate_element(x).is_ok());
        debug_assert!(self.validate_element(y).is_ok());
        #[allow(clippy::match_same_arms)]
        match (x, y) {
            (
                FinitelyFreeSubmoduleAffineSubset::Empty,
                FinitelyFreeSubmoduleAffineSubset::Empty,
            ) => true,
            (
                FinitelyFreeSubmoduleAffineSubset::Empty,
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(..),
                FinitelyFreeSubmoduleAffineSubset::Empty,
            ) => false,
            (
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(x_coset),
                FinitelyFreeSubmoduleAffineSubset::NonEmpty(y_coset),
            ) => self.module().cosets().equal(x_coset, y_coset),
        }
    }
}
