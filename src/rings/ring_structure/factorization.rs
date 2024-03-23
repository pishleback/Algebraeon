use std::{fmt::Display, rc::Rc};

use malachite_base::num::basic::traits::Zero;
use malachite_nz::natural::Natural;

use crate::sets::set::Set;

use super::super::structure::*;
use super::structure::*;

#[derive(Debug)]
pub struct Factored<RS: UniqueFactorizationStructure> {
    ring: Rc<RS>,
    unit: RS::Set,
    factors: Vec<(RS::Set, Natural)>,
}

impl<RS: UniqueFactorizationStructure> Display for Factored<RS>
where
    RS::Set: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.unit);
        for (factor, k) in &self.factors {
            write!(f, " * (");
            write!(f, "{}", factor);
            write!(f, ")");
            if k != &Natural::from(1u8) {
                write!(f, "^");
                write!(f, "{}", k);
            }
        }
        Ok(())
    }
}

impl<RS: UniqueFactorizationStructure> Factored<RS> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        if !self.ring.is_unit(&self.unit) {
            return Err("unit must be a unit");
        }
        for (p, k) in &self.factors {
            if k == &0 {
                return Err("prime powers must not be zero");
            }
            if !self.ring.is_fav_assoc(p) {
                return Err("prime factor must be their favoriate associate");
            }
            if !self.ring.is_irreducible(p) {
                return Err("prime factor must be irreducible");
            }

            let mut i = Natural::from(0u8);
            while &i < k {
                i += Natural::from(1u8);
            }
        }
        Ok(())
    }

    //change to a new ring structure type
    //the new ring structure type must represent the same ring structure
    pub fn change_ring_unsafe<NewRS: UniqueFactorizationStructure<Set = RS::Set>>(
        self,
        ring: Rc<NewRS>,
    ) -> Factored<NewRS> {
        Factored {
            ring,
            unit: self.unit,
            factors: self.factors,
        }
    }

    pub fn expand(&self) -> RS::Set {
        let mut ans = self.unit.clone();
        for (p, k) in &self.factors {
            self.ring.mul_mut(&mut ans, &self.ring.nat_pow(p, k));
        }
        ans
    }

    pub fn new_unchecked(ring: Rc<RS>, unit: RS::Set, factors: Vec<(RS::Set, Natural)>) -> Self {
        Self {
            ring,
            unit,
            factors,
        }
    }

    pub fn equal(a: &Self, b: &Self) -> bool {
        let ring = common_structure(&a.ring, &b.ring);
        if !ring.equal(a.unit(), b.unit()) {
            false
        } else {
            //every a_factor is a b_factor
            for (a_factor, a_power) in a.factors() {
                debug_assert!(ring.is_fav_assoc(a_factor));
                if !b
                    .factors()
                    .into_iter()
                    .any(|(b_factor, b_power)| ring.equal(a_factor, b_factor))
                {
                    return false;
                }
            }
            //every b_factor is an a_factor
            for (b_factor, b_power) in b.factors() {
                debug_assert!(ring.is_fav_assoc(b_factor));
                if !a
                    .factors()
                    .into_iter()
                    .any(|(a_factor, a_power)| ring.equal(a_factor, b_factor))
                {
                    return false;
                }
            }
            //the powers of the factors are equal
            for (a_factor, a_power) in a.factors() {
                for (b_factor, b_power) in b.factors() {
                    if ring.equal(a_factor, b_factor) {
                        if a_power != b_power {
                            return false;
                        }
                    }
                }
            }
            true
        }
    }

    pub fn unit(&self) -> &RS::Set {
        &self.unit
    }

    pub fn factors(&self) -> &Vec<(RS::Set, Natural)> {
        &self.factors
    }

    pub fn unit_and_factors(self) -> (RS::Set, Vec<(RS::Set, Natural)>) {
        (self.unit, self.factors)
    }

    pub fn is_irreducible(&self) -> bool {
        if self.factors.len() == 1 {
            let (_p, k) = self.factors.iter().next().unwrap();
            if k == &1 {
                return true;
            }
        }
        false
    }

    fn mul_by_unchecked(&mut self, p: RS::Set, k: Natural) {
        for (q, t) in &mut self.factors {
            match (self.ring.div(&p, q), self.ring.div(q, &p)) {
                (Ok(u), Ok(v)) => {
                    if self.ring.is_unit(&u) && self.ring.is_unit(&v) {
                        //q = v*p so q^k = v^kp^k and this is what we are multiplying by
                        self.ring
                            .mul_mut(&mut self.unit, &self.ring.nat_pow(&v, &k));
                        *t += k;
                        return;
                    }
                }
                _ => {}
            }
        }
        self.factors.push((p, k));
    }

    pub fn mul_mut(&mut self, other: Self) {
        self.ring.mul_mut(&mut self.unit, &other.unit);
        for (p, k) in other.factors {
            self.mul_by_unchecked(p, k);
        }
    }

    pub fn mul(mut a: Self, b: Self) -> Self {
        let ring = common_structure(&a.ring, &b.ring);
        ring.mul_mut(&mut a.unit, &b.unit);
        for (p, k) in b.factors {
            a.mul_by_unchecked(p, k)
        }
        a
    }

    pub fn pow(mut self, k: &Natural) -> Self {
        for (factor, power) in self.factors.iter_mut() {
            *power *= k;
        }
        self
    }

    pub fn factored_one(ring: Rc<RS>) -> Self {
        Factored {
            unit: ring.one(),
            ring,
            factors: vec![],
        }
    }

    pub fn factored_irreducible_unchecked(ring: Rc<RS>, elem: RS::Set) -> Self {
        let (unit, assoc) = ring.factor_fav_assoc(&elem);
        Factored {
            ring,
            unit,
            factors: vec![(assoc, Natural::from(1u8))],
        }
    }

    pub fn factored_irreducible_power_unchecked(ring: Rc<RS>, elem: RS::Set, k: Natural) -> Self {
        debug_assert!(k >= 1);
        let (unit, assoc) = ring.factor_fav_assoc(&elem);
        Factored {
            ring,
            unit,
            factors: vec![(assoc, k)],
        }
    }

    pub fn factored_unit_unchecked(ring: Rc<RS>, unit: RS::Set) -> Self {
        Factored {
            ring,
            unit,
            factors: vec![],
        }
    }

    pub fn divisors<'a>(&'a self) -> Box<dyn Iterator<Item = RS::Set> + 'a> {
        if self.factors.len() == 0 {
            Box::new(vec![self.ring.one()].into_iter())
        } else {
            let mut factor_powers = vec![];
            for (p, k) in &self.factors {
                let j = factor_powers.len();
                factor_powers.push(vec![]);
                let mut p_pow = self.ring.one();
                let mut i = Natural::from(0u8);
                while &i <= k {
                    factor_powers[j].push(p_pow.clone());
                    p_pow = self.ring.mul(&p_pow, &p);
                    i += Natural::from(1u8);
                }
            }

            Box::new(
                itertools::Itertools::multi_cartesian_product(
                    factor_powers.into_iter().map(|p_pows| p_pows.into_iter()),
                )
                .map(|prime_power_factors| {
                    self.ring
                        .product(prime_power_factors.iter().collect())
                        .clone()
                }),
            )
        }
    }

    pub fn count_divisors(&self) -> Option<Natural> {
        let mut count = Natural::from(1u8);
        for (_p, k) in &self.factors {
            count *= k + Natural::from(1u8);
        }
        Some(count)
    }
}

pub fn factorize_by_find_factor<RS: UniqueFactorizationStructure>(
    ring: &RS,
    elem: RS::Set,
    partial_factor: &impl Fn(RS::Set) -> Option<(RS::Set, RS::Set)>,
) -> Factored<RS> {
    debug_assert!(!ring.is_zero(&elem));
    if ring.is_unit(&elem) {
        Factored::factored_unit_unchecked(ring.clone().into(), elem)
    } else {
        debug_assert!(!ring.is_unit(&elem));
        match partial_factor(elem.clone()) {
            Some((g, h)) => Factored::mul(
                factorize_by_find_factor(ring, g, partial_factor),
                factorize_by_find_factor(ring, h, partial_factor),
            ),
            None => {
                //f is irreducible
                Factored::factored_irreducible_unchecked(ring.clone().into(), elem)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use super::super::cannonical::*;
    use super::*;

    #[test]
    fn factorization_invariants() {
        let f = Factored::new_unchecked(
            Integer::structure(),
            Integer::from(-1),
            vec![
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        f.check_invariants().unwrap();

        let f = Factored::new_unchecked(Integer::structure(), Integer::from(1), vec![]);
        f.check_invariants().unwrap();

        let f = Factored::new_unchecked(
            Integer::structure(),
            Integer::from(-1),
            vec![
                (Integer::from(2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
                (Integer::from(5), Natural::from(0u8)),
            ],
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "can't have a power of zero"
        );

        let f = Factored::new_unchecked(
            Integer::structure(),
            Integer::from(3),
            vec![(Integer::from(2), Natural::from(2u8))],
        );
        assert_eq!(f.check_invariants().is_ok(), false, "unit should be a unit");

        let f = Factored::new_unchecked(
            Integer::structure(),
            Integer::from(1),
            vec![
                (Integer::from(0), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "prime factors must not be zero"
        );

        let f = Factored::new_unchecked(
            Integer::structure(),
            Integer::from(-1),
            vec![
                (Integer::from(4), Natural::from(1u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "prime factors must be prime"
        );

        let f = Factored::new_unchecked(
            Integer::structure(),
            Integer::from(-1),
            vec![
                (Integer::from(-2), Natural::from(2u8)),
                (Integer::from(3), Natural::from(1u8)),
            ],
        );
        assert_eq!(
            f.check_invariants().is_ok(),
            false,
            "prime factors must be favoriate associate"
        );
    }

    #[test]
    fn test_count_divisors() {
        for a in 1..25 {
            let b = Integer::from(a);
            let fs = Integer::structure().factor(&b).unwrap();
            println!("{}", fs);
            assert_eq!(
                fs.count_divisors().unwrap(),
                Natural::from(fs.divisors().collect::<Vec<Integer>>().len())
            );
        }
    }
}
