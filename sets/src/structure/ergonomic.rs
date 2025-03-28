use std::rc::Rc;

use super::structure::*;

pub trait IntoErgonomic: MetaType {
    fn into_ergonomic(self) -> StructuredElement<Self::Structure> {
        StructuredElement::new(Self::structure(), self)
    }
}
impl<T: MetaType> IntoErgonomic for T {}

#[derive(Debug, Clone)]
pub struct StructuredElement<S: Structure> {
    structure: Rc<S>,
    elem: S::Set,
}

impl<S: Structure> StructuredElement<S> {
    pub fn new(structure: Rc<S>, elem: S::Set) -> Self {
        Self { structure, elem }
    }

    pub fn structure(&self) -> Rc<S> {
        self.structure.clone()
    }

    pub fn ref_set(&self) -> &S::Set {
        &self.elem
    }

    pub fn into_verbose(self) -> S::Set {
        self.elem
    }
}

impl<S: Structure> std::fmt::Display for StructuredElement<S>
where
    S::Set: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.elem, f)
    }
}

impl<S: PartialEqStructure> PartialEq for StructuredElement<S> {
    fn eq(&self, other: &Self) -> bool {
        let structure = common_structure(self.structure(), other.structure());
        structure.equal(&self.ref_set(), &other.ref_set())
    }
}

impl<S: EqStructure> Eq for StructuredElement<S> {}
