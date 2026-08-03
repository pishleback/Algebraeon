use crate::{SetSignature, Signature};
use std::rc::Rc;

impl<A: Signature, B: Signature> Signature for (Rc<A>, Rc<B>) {}

impl<A: SetSignature, B: SetSignature> SetSignature for (Rc<A>, Rc<B>) {
    type Elem = (A::Elem, B::Elem);

    fn validate_element(self: &Rc<Self>, x: &Self::Elem) -> Result<(), String> {
        self.0.validate_element(&x.0)?;
        self.1.validate_element(&x.1)?;
        Ok(())
    }
}
