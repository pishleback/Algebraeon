use algebraeon_nzq::Rational;

/// The open interval (a, b) in the real line.
/// Must have a < b
#[derive(Debug, Clone)]
pub struct RationalInterval {
    a: Rational,
    b: Rational,
}

impl RationalInterval {
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

    #[allow(clippy::should_implement_trait)]
    pub fn neg(self) -> RationalInterval {
        Self {
            a: -self.b,
            b: -self.a,
        }
    }
}
