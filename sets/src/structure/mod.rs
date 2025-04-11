//! Abtractions over sets with certain structure.
//!
//! The structure framework used by algebraeon_rings is established here.

mod morphism;
mod structure;

pub use algebraeon_canonical_structure_derive::CanonicalStructure;
pub use morphism::{
    BijectiveFunction, CompositionMorphism, Function, Functions, IdentityMorphism,
    InjectiveFunction, Morphism, MorphismsStructure,
};
pub use structure::{
    CountableSetStructure, EqStructure, FiniteSetStructure, MetaType, SetStructure, Structure,
    ToStringStructure, common_structure,
};
