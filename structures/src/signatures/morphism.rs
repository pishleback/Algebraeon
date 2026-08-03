use crate::*;
use std::rc::Rc;

pub trait Morphism<Domain: Signature, Range: Signature>: Debug + Clone {
    fn domain(self: &Rc<Self>) -> Rc<Domain>;
    fn range(self: &Rc<Self>) -> Rc<Range>;
}

/// A morphism from an object to itself
pub trait Endomorphism<X: Signature>: Morphism<X, X> {}
impl<X: Signature, T: Morphism<X, X>> Endomorphism<X> for T {}
