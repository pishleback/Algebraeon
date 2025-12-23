use crate::open_box::OpenRationalBox;

use super::*;
use std::fmt::Debug;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ComplexPlaneIsolatingBoxesStructure {}

impl Signature for ComplexPlaneIsolatingBoxesStructure {}

impl SetSignature for ComplexPlaneIsolatingBoxesStructure {
    type Set = OpenRationalBox;

    fn is_element(&self, _: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl TopologicalSpaceOpenSubsetsSignature for ComplexPlaneIsolatingBoxesStructure {}
