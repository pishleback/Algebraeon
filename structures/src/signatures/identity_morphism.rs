use crate::*;
use std::rc::Rc;

/// The identity morphism X -> X
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdentityMorphism<X: Signature> {
    x: Rc<X>,
}

impl<X: Signature> IdentityMorphism<X> {
    pub fn new(x: Rc<X>) -> Self {
        Self { x }
    }
}

impl<X: Signature> Morphism<X, X> for IdentityMorphism<X> {
    fn domain(self: &Rc<Self>) -> Rc<X> {
        self.x.clone()
    }

    fn range(self: &Rc<Self>) -> Rc<X> {
        self.x.clone()
    }
}

impl<X: SetSignature> FunctionMorphism<X, X> for IdentityMorphism<X> {
    fn image(self: &Rc<Self>, x: &X::Elem) -> X::Elem {
        x.clone()
    }
}

impl<X: SetSignature> InjectiveFunctionMorphism<X, X> for IdentityMorphism<X> {
    fn try_preimage(self: &Rc<Self>, x: &X::Elem) -> Option<X::Elem> {
        Some(x.clone())
    }
}

impl<X: SetSignature> BijectiveFunctionMorphism<X, X> for IdentityMorphism<X> {}
