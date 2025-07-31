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
    clippy::wildcard_imports,
    clippy::type_complexity
)]

use algebraeon_rings::structure::{FieldSignature, OrderedRingSignature};
use std::borrow::Borrow;
use std::hash::Hash;

pub mod affine_subspace;
pub mod ambient_space;
pub mod boolean_operations;
pub mod convex_hull;
pub mod coordinates;
pub mod partial_simplicial_complex;
pub mod simplex;
pub mod simplex_collection;
pub mod simplicial_complex;
pub mod simplicial_disjoint_union;
