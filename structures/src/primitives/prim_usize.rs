use crate::*;
use std::cmp::Ordering;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimitiveUsizeCanonicalStructure {}

impl MetaType for usize {
    type Signature = PrimitiveUsizeCanonicalStructure;

    fn structure() -> Self::Signature {
        PrimitiveUsizeCanonicalStructure {}
    }
}

impl Signature for PrimitiveUsizeCanonicalStructure {}

impl SetSignature for PrimitiveUsizeCanonicalStructure {
    type Elem = usize;

    fn validate_element(&self, _x: &usize) -> Result<(), String> {
        Ok(())
    }
}

impl EqSignature for PrimitiveUsizeCanonicalStructure {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        a == b
    }
}

impl PartialOrdSignature for PrimitiveUsizeCanonicalStructure {
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        Some(self.cmp(a, b))
    }
}

impl OrdSignature for PrimitiveUsizeCanonicalStructure {
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        a.cmp(b)
    }
}
