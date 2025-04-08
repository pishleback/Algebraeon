#![allow(dead_code)]

use algebraeon_rings::structure::{FieldStructure, OrderedRingStructure};
use std::borrow::Borrow;
use std::hash::Hash;

mod coordinates;
pub use coordinates::*;

mod ambient_space;
pub use ambient_space::*;

mod affine_subspace;
pub use affine_subspace::*;

pub mod simplexes;
