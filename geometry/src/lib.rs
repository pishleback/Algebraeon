#![allow(dead_code)]
#![allow(
    clippy::uninlined_format_args,
    clippy::to_string_in_format_args,
    clippy::must_use_candidate,
    clippy::return_self_not_must_use,
    clippy::missing_panics_doc,
    clippy::missing_errors_doc,
    clippy::doc_markdown,
    clippy::too_many_lines,
    clippy::similar_names,
    clippy::many_single_char_names,
    clippy::wrong_self_convention,
    clippy::from_over_into,
    clippy::wildcard_imports
)]

use algebraeon_rings::structure::{FieldSignature, OrderedRingSignature};
use std::borrow::Borrow;
use std::hash::Hash;

mod coordinates;
pub use coordinates::*;

mod ambient_space;
pub use ambient_space::*;

mod affine_subspace;
pub use affine_subspace::*;

pub mod simplexes;
