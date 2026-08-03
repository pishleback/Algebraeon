//! Given a 6-element set make the following definitions
//!  - A duad is a 2-element subset
//!  - A syntheme is an unordered set of 3 disjoint duads
//!  - A pentad is an unordered set of 5 disjoint synthemes
//!
//! Then there are 6 pentads, and any permutation of the 6 points induces a permutation of the 6 pentads via an outer automorphism of S6

mod duads;
mod pentads;
mod synthemes;

pub use duads::*;
pub use pentads::*;
pub use synthemes::*;
