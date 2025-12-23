use algebraeon_nzq::Rational;

/// The open interval (a, b) in the real line.
/// Must have a < b
#[derive(Debug, Clone)]
pub struct OpenRationalInterval {
    a: Rational,
    b: Rational,
}
