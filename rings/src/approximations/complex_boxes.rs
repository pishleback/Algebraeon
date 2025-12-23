use algebraeon_nzq::Rational;
use algebraeon_sets::{
    approximations::SubsetsSignature,
    structure::{SetSignature, Signature},
};

use crate::approximations::open_box::OpenRationalBox;

use super::*;
use std::fmt::Debug;

#[derive(Debug, Clone)]
pub enum Subset {
    Singleton { real: Rational, imag: Rational },
    Box(OpenRationalBox),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BoxesStructure {}

impl Signature for BoxesStructure {}

impl SetSignature for BoxesStructure {
    type Set = Subset;

    fn is_element(&self, _: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl SubsetsSignature for BoxesStructure {}
