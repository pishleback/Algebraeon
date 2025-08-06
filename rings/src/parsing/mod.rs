#![allow(clippy::match_same_arms, clippy::unnested_or_patterns)]
mod ast;
mod polynomial;

pub use polynomial::parse_integer_polynomial;
pub use polynomial::parse_multivariate_integer_polynomial;
pub use polynomial::parse_multivariate_rational_polynomial;
pub use polynomial::parse_rational_polynomial;

#[cfg(test)]
mod polynomial_tests;
