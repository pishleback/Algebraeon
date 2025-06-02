//! Abstractions over sets with certain structure.
//!
//! The structure framework used by algebraeon_rings is established here.

mod morphism;
mod structure;

pub use algebraeon_macros::CanonicalStructure;
pub use morphism::{
    BijectiveFunction, CompositionMorphism, Function, Functions, IdentityMorphism,
    InjectiveFunction, Morphism,
};
pub use structure::{
    BorrowedStructure, CountableSetSignature, EqSignature, FiniteSetSignature, MetaType,
    SetSignature, Signature, ToStringSignature, common_structure,
};
