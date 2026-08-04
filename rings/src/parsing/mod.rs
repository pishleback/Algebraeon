#![allow(clippy::match_same_arms, clippy::unnested_or_patterns)]
mod algebra;
mod ast;
mod polynomial;
mod quaternion;

pub use algebra::parse_algebra_expression;
pub use polynomial::parse_integer_polynomial;
pub use polynomial::parse_multivariate_integer_polynomial;
pub use polynomial::parse_multivariate_rational_polynomial;
pub use polynomial::parse_rational_polynomial;
pub use quaternion::parse_quaternion;

#[cfg(test)]
mod algebra_tests;
#[cfg(test)]
mod polynomial_tests;
#[cfg(test)]
mod quaternion_tests;
