use crate::{SetSignature, Signature};
use std::sync::Arc;

impl<A: Signature, B: Signature> Signature for (Arc<A>, Arc<B>) {}

impl<A: SetSignature, B: SetSignature> SetSignature for (Arc<A>, Arc<B>) {
    type Elem = (A::Elem, B::Elem);

    fn validate_element(self: &Arc<Self>, x: &Self::Elem) -> Result<(), String> {
        self.0.validate_element(&x.0)?;
        self.1.validate_element(&x.1)?;
        Ok(())
    }
}
