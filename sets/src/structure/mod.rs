//! Abstractions over sets with certain structure.
//!
//! The structure framework used by `algebraeon_rings` is established here.

mod empty_set;
mod finite_set;
mod morphism;
mod pairs;
mod singleton_set;

pub use empty_set::EmptySetStructure;
pub use finite_set::EnumeratedFiniteSetStructure;
pub use morphism::{FiniteSetEndofunctions, Functions};
pub use pairs::{PairsStructure, UnorderedPair, UnorderedPairs};
pub use singleton_set::SingletonSetStructure;
