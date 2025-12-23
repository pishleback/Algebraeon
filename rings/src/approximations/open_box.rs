use algebraeon_nzq::Rational;

/// The open box in the complex plane of complex numbers z such that a < Re(z) < b and c < Im(z) < d
/// Must have a < b and c < d
#[derive(Debug, Clone)]
pub struct OpenRationalBox {
    a: Rational,
    b: Rational,
    c: Rational,
    d: Rational,
}
