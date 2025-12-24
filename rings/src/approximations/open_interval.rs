use algebraeon_nzq::Rational;

/// The open interval (a, b) in the real line.
/// Must have a < b
#[derive(Debug, Clone)]
pub struct OpenRationalInterval {
    a: Rational,
    b: Rational,
}

impl OpenRationalInterval {
    pub fn a(&self) -> &Rational {
        &self.a
    }

    pub fn b(&self) -> &Rational {
        &self.b
    }

    pub fn new_unchecked(a: Rational, b: Rational) -> Self {
        debug_assert!(a < b);
        Self { a, b }
    }
}
