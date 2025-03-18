//Given a pair of rational numbers a and b, how to iterator over possible simple middle rational numbers m?

use std::collections::BTreeMap;

use algebraeon_nzq::rational::*;

#[derive(Debug, Clone)]
pub struct RationalSimpleBetweenGenerator {
    //store intervals to check, stored by length
    regions: BTreeMap<Rational, Vec<(Rational, Rational)>>,
}

impl RationalSimpleBetweenGenerator {
    pub fn new_between(a: Rational, b: Rational) -> Self {
        assert!(a < b);
        Self {
            regions: BTreeMap::from([(&b - &a, vec![(a, b)])]),
        }
    }

    //between the middle f fraction of a and b
    pub fn new_between_middle(a: &Rational, b: &Rational, f: &Rational) -> Self {
        assert!(a < b);
        assert!(&Rational::ZERO < f);
        assert!(f <= &Rational::ONE);
        let d = &(b - a);
        let mid_a = a + (d - d * f) / Rational::TWO;
        let mid_b = b - (d - d * f) / Rational::TWO;
        Self::new_between(mid_a, mid_b)
    }

    fn add_region(&mut self, a: Rational, b: Rational) {
        let key = &b - &a;
        match self.regions.get_mut(&key) {
            Some(r) => r.push((a, b)),
            None => {
                self.regions.insert(key, vec![(a, b)]);
            }
        }
    }
}

impl Iterator for RationalSimpleBetweenGenerator {
    type Item = Rational;

    fn next(&mut self) -> Option<Self::Item> {
        Some({
            let key = self.regions.last_entry().unwrap().key().clone();
            let (a, b) = self.regions.get_mut(&key).unwrap().pop().unwrap();
            if self.regions.get(&key).unwrap().is_empty() {
                self.regions.remove(&key);
            }
            let m = Rational::simplest_rational_in_open_interval(&a, &b);
            self.add_region(a, m.clone());
            self.add_region(m.clone(), b);
            m
        })
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_bisection_gen() {
        let a = Rational::from_str("1/4").unwrap();
        let b = Rational::from_str("1").unwrap();
        let mut g = RationalSimpleBetweenGenerator::new_between(a, b);

        assert_eq!(g.next(), Some(Rational::from_str("1/2").unwrap()));
        assert_eq!(g.next(), Some(Rational::from_str("2/3").unwrap()));
        assert_eq!(g.next(), Some(Rational::from_str("3/4").unwrap()));
        assert_eq!(g.next(), Some(Rational::from_str("4/5").unwrap()));
        assert_eq!(g.next(), Some(Rational::from_str("1/3").unwrap()));
        assert_eq!(g.next(), Some(Rational::from_str("5/6").unwrap()));
        assert_eq!(g.next(), Some(Rational::from_str("6/7").unwrap()));
        assert_eq!(g.next(), Some(Rational::from_str("2/5").unwrap()));
    }
}
