use crate::number::rational::*;

#[derive(Debug, Clone)]
pub enum LowerBound {
    Inf,
    Finite(Rational),
}

#[derive(Debug, Clone)]
pub enum UpperBound {
    Inf,
    Finite(Rational),
}

impl LowerBound {
    pub fn neg(self) -> UpperBound {
        match self {
            LowerBound::Inf => UpperBound::Inf,
            LowerBound::Finite(a) => UpperBound::Finite(-a),
        }
    }
}

impl UpperBound {
    pub fn neg(self) -> LowerBound {
        match self {
            UpperBound::Inf => LowerBound::Inf,
            UpperBound::Finite(a) => LowerBound::Finite(-a),
        }
    }
}

impl std::hash::Hash for LowerBound {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            LowerBound::Inf => {}
            LowerBound::Finite(x) => {
                x.hash(state);
            }
        }
    }
}

impl std::hash::Hash for UpperBound {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        match self {
            UpperBound::Inf => {}
            UpperBound::Finite(x) => {
                x.hash(state);
            }
        }
    }
}

impl PartialEq<UpperBound> for LowerBound {
    fn eq(&self, other: &UpperBound) -> bool {
        match (self, other) {
            (LowerBound::Finite(a), UpperBound::Finite(b)) => a == b,
            _ => false,
        }
    }
}

impl PartialOrd<UpperBound> for LowerBound {
    fn partial_cmp(&self, other: &UpperBound) -> Option<std::cmp::Ordering> {
        match (self, other) {
            (LowerBound::Finite(a), UpperBound::Finite(b)) => a.partial_cmp(b),
            _ => Some(std::cmp::Ordering::Less),
        }
    }
}

impl PartialEq<Rational> for LowerBound {
    fn eq(&self, b: &Rational) -> bool {
        match self {
            LowerBound::Finite(a) => a == b,
            _ => false,
        }
    }
}

impl PartialOrd<Rational> for LowerBound {
    fn partial_cmp(&self, b: &Rational) -> Option<std::cmp::Ordering> {
        match self {
            LowerBound::Finite(a) => a.partial_cmp(b),
            _ => Some(std::cmp::Ordering::Less),
        }
    }
}

impl PartialEq<UpperBound> for Rational {
    fn eq(&self, other: &UpperBound) -> bool {
        match other {
            UpperBound::Finite(b) => self == b,
            _ => false,
        }
    }
}

impl PartialOrd<UpperBound> for Rational {
    fn partial_cmp(&self, other: &UpperBound) -> Option<std::cmp::Ordering> {
        match other {
            UpperBound::Finite(b) => self.partial_cmp(b),
            _ => Some(std::cmp::Ordering::Less),
        }
    }
}
