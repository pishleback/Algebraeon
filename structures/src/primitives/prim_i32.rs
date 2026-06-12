use crate::*;
use std::cmp::Ordering;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimitiveI32CanonicalStructure {}

impl MetaType for i32 {
    type Signature = PrimitiveI32CanonicalStructure;

    fn structure() -> Self::Signature {
        PrimitiveI32CanonicalStructure {}
    }
}

impl Signature for PrimitiveI32CanonicalStructure {}

impl SetSignature for PrimitiveI32CanonicalStructure {
    type Elem = i32;

    fn validate_element(&self, _x: &i32) -> Result<(), String> {
        Ok(())
    }
}

impl EqSignature for PrimitiveI32CanonicalStructure {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        a == b
    }
}

impl PartialOrdSignature for PrimitiveI32CanonicalStructure {
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        Some(self.cmp(a, b))
    }
}

impl OrdSignature for PrimitiveI32CanonicalStructure {
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        a.cmp(b)
    }
}
