use crate::{SetSignature, Signature};

impl<A: Signature, B: Signature> Signature for (A, B) {}

impl<A: SetSignature, B: SetSignature> SetSignature for (A, B) {
    type Elem = (A::Elem, B::Elem);

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        self.0.validate_element(&x.0)?;
        self.1.validate_element(&x.1)?;
        Ok(())
    }
}
