use crate::*;
use algebraeon_macros::signature_meta_trait;
use std::sync::Arc;

#[signature_meta_trait]
pub trait FunctionsSignature<Domain: SetSignature, Range: SetSignature>: SetSignature {
    fn function(self: &Arc<Self>, f: impl Fn(&Domain::Elem) -> Range::Elem) -> Option<Self::Elem>;
    fn image<'a>(self: &Arc<Self>, f: &'a Self::Elem, x: &Domain::Elem) -> &'a Range::Elem;
}
