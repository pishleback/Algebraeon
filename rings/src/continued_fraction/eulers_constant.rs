use crate::continued_fraction::IrrationalSimpleContinuedFractionGenerator;
use algebraeon_nzq::Integer;

#[derive(Debug, Clone)]
pub struct EulersConstantSimpleContinuedFractionGenerator {
    n: usize,
}
#[allow(clippy::new_without_default)]
impl EulersConstantSimpleContinuedFractionGenerator {
    pub fn new() -> Self {
        Self { n: 0 }
    }
}
impl IrrationalSimpleContinuedFractionGenerator for EulersConstantSimpleContinuedFractionGenerator {
    fn next(&mut self) -> Integer {
        let c = if self.n == 0 {
            Integer::TWO
        } else if self.n % 3 == 2 {
            Integer::from(2 * (self.n / 3 + 1))
        } else {
            Integer::ONE
        };
        self.n += 1;
        c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eulers_constant_continued_fraction() {
        let mut e = EulersConstantSimpleContinuedFractionGenerator::new();
        assert_eq!(e.next(), Integer::from(2));
        assert_eq!(e.next(), Integer::from(1));
        assert_eq!(e.next(), Integer::from(2));
        assert_eq!(e.next(), Integer::from(1));
        assert_eq!(e.next(), Integer::from(1));
        assert_eq!(e.next(), Integer::from(4));
        assert_eq!(e.next(), Integer::from(1));
        assert_eq!(e.next(), Integer::from(1));
        assert_eq!(e.next(), Integer::from(6));
    }
}
