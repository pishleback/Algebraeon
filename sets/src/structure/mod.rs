//! Abtractions over sets with certain structure.
//!
//! The structure framework used by algebraeon_rings is established here.

mod morphism;
mod structure;

pub use algebraeon_cannonical_structure_derive::CannonicalStructure;
pub use morphism::{
    BijectiveFunctionStructure, CompositionMorphism, FunctionStructure, IdentityMorphism,
    InjectiveFunctionStructure, Morphism, MorphismStructure,
};
pub use structure::{
    CannonicalStructure, CountableSetStructure, EqStructure, FiniteSetStructure, MetaType,
    PartialEqStructure, SetStructure, Structure, ToStringStructure, common_structure,
};
