//! Abtractions over sets with certain structure.
//!
//! The structure framework used by algebraeon_rings is established here.

mod morphism;
mod structure;

pub use morphism::{
    BijectiveFunctionStructure, CompositionMorphism, FunctionStructure, IdentityMorphism,
    InjectiveFunctionStructure, Morphism, MorphismStructure,
};
pub use structure::{
    CannonicalStructure, EqStructure, MetaType, PartialEqStructure, SetStructure, Structure,
    ToStringStructure, common_structure,
};
