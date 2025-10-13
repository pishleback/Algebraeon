#![allow(
    clippy::must_use_candidate,
    clippy::uninlined_format_args,
    clippy::cast_lossless,
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::implied_bounds_in_impls,
    clippy::needless_lifetimes,
    clippy::unreadable_literal
)]
#![allow(dead_code)]

pub mod combinatorics;
mod integer;
mod natural;
mod random;
mod rational;
pub mod traits;

pub use integer::{Integer, IntegerCanonicalStructure};
pub use natural::choose;
pub use natural::gcd;
pub use natural::lcm;
pub use natural::primes;
pub use natural::{Natural, NaturalCanonicalStructure};
pub use random::Rng;
pub use rational::{Rational, RationalCanonicalStructure};
