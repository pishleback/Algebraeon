//! Abtractions over sets with certain structure.
//!
//! The structure framework used by `algebraeon_rings` is established here.

mod empty_set;
mod morphism;
mod orderings;
mod structure;

pub use algebraeon_macros::CanonicalStructure;
pub use empty_set::EmptySetStructure;
pub use morphism::{
    BijectiveFunction, CompositionMorphism, Function, Functions, IdentityMorphism,
    InjectiveFunction, Morphism,
};
pub use orderings::OrdSignature;
pub use structure::{
    BorrowedStructure, CountableSetSignature, EqSignature, FiniteSetSignature, MetaType,
    SetSignature, Signature, ToStringSignature, common_structure,
};
