//! Abtractions over sets with certain structure.
//!
//! The structure framework used by algebraeon_rings is established here.

mod structure;

pub use structure::{
    CannonicalStructure, EqStructure, MetaType, PartialEqStructure, Structure, ToStringStructure,
    common_structure,
};
