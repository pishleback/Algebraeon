use algebraeon_structure::*;
use core::fmt::Debug;
use std::rc::Rc;

pub trait IntoRingElem: MetaType {
    fn into_ring(self) -> StructuredElement<Self::Structure> {
        StructuredElement::new(Self::structure(), self)
    }
}
impl<T: MetaType> IntoRingElem for T {}

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

    pub fn into_set(self) -> S::Set {
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
