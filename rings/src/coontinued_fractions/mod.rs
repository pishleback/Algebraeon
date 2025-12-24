use algebraeon_nzq::Integer;
use std::fmt::Debug;

pub trait SimpleContinuedFraction: Debug + Clone + Send + Sync {
    /// First value returned can be any integer
    /// All subsequent values must be >= 1
    fn next(&mut self) -> Integer;
}

mod eulers_constant;

pub fn eulers_constant() -> impl SimpleContinuedFraction {
    eulers_constant::EulersConstantSimpleContinuedFraction::new()
}
