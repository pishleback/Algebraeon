use crate::*;
use algebraeon_macros::signature_meta_trait;
use std::rc::Rc;

#[signature_meta_trait]
pub trait FunctionsSignature<Domain: SetSignature, Range: SetSignature>: SetSignature {
    fn function(self: &Rc<Self>, f: impl Fn(&Domain::Elem) -> Range::Elem) -> Option<Self::Elem>;
    fn image<'a>(self: &Rc<Self>, f: &'a Self::Elem, x: &Domain::Elem) -> &'a Range::Elem;
}
