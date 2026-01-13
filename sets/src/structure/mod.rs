//! Abstractions over sets with certain structure.
//!
//! The structure framework used by `algebraeon_rings` is established here.

mod empty_set;
mod finite_set;
mod morphism;
mod orderings;
mod pairs;
mod singleton_set;
#[allow(clippy::module_inception)]
mod structure;

pub use algebraeon_macros::CanonicalStructure;
pub use empty_set::EmptySetStructure;
pub use finite_set::EnumeratedFiniteSetStructure;
pub use morphism::{
    BijectiveFunction, BorrowedMorphism, CompositionMorphism, Endofunction, Endomorphism,
    FiniteSetEndofunctions, Function, Functions, IdentityMorphism, InjectiveFunction, Morphism,
    Permutation,
};
pub use orderings::{OrdSignature, PartialOrdSignature};
pub use pairs::{PairsStructure, UnorderedPair, UnorderedPairs};
pub use singleton_set::SingletonSetStructure;
pub use structure::{
    BorrowedSet, BorrowedStructure, CountableSetSignature, EqSignature, FiniteSetSignature,
    MaybeFiniteSetSignature, MetaType, SetSignature, Signature, ToStringSignature,
};
