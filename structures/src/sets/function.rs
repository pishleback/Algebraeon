use std::rc::Rc;

use crate::*;

pub trait FunctionMorphism<Domain: SetSignature, Range: SetSignature>:
    Morphism<Domain, Range>
{
    fn image(self: &Rc<Self>, x: &Domain::Elem) -> Range::Elem;
}

/// A function from a set into itself
pub trait EndofunctionMorphism<X: SetSignature + EqSignature>: FunctionMorphism<X, X> {
    // TODO: remove EqSignature requirement and use specialization once it is stable.
    /// check if an element is fixed
    fn is_fixed_point(self: &Rc<Self>, x: X::Elem) -> bool {
        self.domain().equal(&self.image(&x), &x)
    }
}

pub trait InjectiveFunctionMorphism<Domain: SetSignature, Range: SetSignature>:
    FunctionMorphism<Domain, Range>
{
    fn try_preimage(self: &Rc<Self>, y: &Range::Elem) -> Option<Domain::Elem>;
}

pub trait BijectiveFunctionMorphism<Domain: SetSignature, Range: SetSignature>:
    InjectiveFunctionMorphism<Domain, Range>
{
    fn preimage(self: &Rc<Self>, y: &Range::Elem) -> Domain::Elem {
        self.try_preimage(y).unwrap()
    }
}
