use crate::{
    number::algebraic::number_field::new_anf, polynomial::polynomial::*,
    ring_structure::structure::*,
};
use algebraeon_sets::structure::*;
use bounds::*;
use interval::*;
use malachite_base::num::basic::traits::{One, Two, Zero};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;
use polynomial::*;
use std::{fmt::Display, rc::Rc};

use super::{bisection_gen::RationalSimpleBetweenGenerator, poly_tools::*, rat_to_string};

mod bounds;
mod interval;
pub mod polynomial;

#[derive(Debug, Clone)]
pub struct RealAlgebraicRoot {
    poly: Polynomial<Integer>, //a primitive irreducible polynomial of degree >= 2 with a unique real root between a and b
    //an arbitrarily small interval containing the root. May be mutated
    tight_a: Rational, //tight lower bound
    tight_b: Rational, //tight upper bound
    //a heuristically large interval containing the root. Should not shrink
    wide_a: LowerBound, //wide lower bound. None means -inf
    wide_b: UpperBound, //wide upper bound. None means +inf
    //false : decreasing i.e. poly(a) > poly(b), true : increasing i.e. poly(a) < poly(b)
    dir: bool,
}

impl std::hash::Hash for RealAlgebraicRoot {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        debug_assert!(self.poly.leading_coeff().unwrap() > 0);
        self.poly.hash(state);
    }
}

impl RealAlgebraicRoot {
    pub fn poly(&self) -> &Polynomial<Integer> {
        &self.poly
    }
    pub fn tight_a(&self) -> &Rational {
        &self.tight_a
    }
    pub fn tight_b(&self) -> &Rational {
        &self.tight_b
    }
}

impl Display for RealAlgebraicRoot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.poly.num_coeffs() == 3 {
            //quadratic
            let a = self.poly.coeff(2);
            let b = self.poly.coeff(1);
            let c = self.poly.coeff(0);
            debug_assert!(a > Integer::ZERO);

            let d = &b * &b - Integer::from(4) * &a * &c;
            let mut d_sq = Integer::ONE;
            let mut d_sqfreee = Integer::ONE;
            let (d_sign, d_factors) = d.factor().unwrap().unit_and_factors();
            for (d_factor, k) in d_factors {
                d_sq *= d_factor.nat_pow(&(&k / Natural::TWO));
                if k % Natural::TWO == Natural::ONE {
                    d_sqfreee *= d_factor;
                }
            }
            debug_assert_eq!(d_sign, Integer::ONE); //because we are a real number
            debug_assert_eq!(d, &d_sqfreee * &d_sq * &d_sq);

            let two_a = Integer::TWO * a;

            let x = Rational::from_integers(-b, two_a.clone());
            let y = Rational::from_integers(d_sq, two_a);
            debug_assert!(y > Rational::ZERO);
            let r = d_sqfreee;

            let mut tight_a_abs = self.tight_a.clone();
            if tight_a_abs < Rational::ZERO {
                tight_a_abs = -tight_a_abs;
            }

            let mut tight_b_abs = self.tight_b.clone();
            if tight_b_abs < Rational::ZERO {
                tight_b_abs = -tight_b_abs;
            }

            let sign = tight_a_abs < tight_b_abs;
            if x == Rational::ZERO {
                if y == Rational::ONE {
                    write!(
                        f,
                        "{}√{}",
                        match sign {
                            true => "",
                            false => "-",
                        },
                        r
                    )?;
                } else {
                    write!(
                        f,
                        "{}{}√{}",
                        match sign {
                            true => "",
                            false => "-",
                        },
                        y,
                        r
                    )?;
                }
            } else {
                if y == Rational::ONE {
                    write!(
                        f,
                        "{}{}√{}",
                        x,
                        match sign {
                            true => "+",
                            false => "-",
                        },
                        r
                    )?;
                } else {
                    write!(
                        f,
                        "{}{}{}√{}",
                        x,
                        match sign {
                            true => "+",
                            false => "-",
                        },
                        y,
                        r
                    )?;
                }
            }
        } else {
            let mut root = self.clone();
            root.refine_to_accuracy(&Rational::from_integers(
                Integer::from(1),
                Integer::from(100),
            ));
            let m = (&root.tight_a + &root.tight_b) / Rational::TWO;

            write!(f, "≈")?;
            write!(f, "{}", rat_to_string(m))?;
            // write!(f, "±");
            // write!(f, "{}", rat_to_string(self.accuracy() / Rational::TWO));
        }
        Ok(())
    }
}

impl RealAlgebraicRoot {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        if !(self.tight_a < self.tight_b) {
            return Err("tight a should be strictly less than b");
        }
        if !(self.wide_a.clone() < self.wide_b.clone()) {
            return Err("wide a should be strictly less than b");
        }
        if self.poly
            != self
                .poly
                .clone()
                .primitive_squarefree_part()
                .factor_fav_assoc()
                .1
        {
            return Err("poly should be primitive and favoriate associate");
        }
        if !self.poly.is_irreducible() {
            return Err("poly should be irreducible");
        }
        if self.poly.degree().unwrap() < 2 {
            return Err("poly should have degree at least 2");
        }
        let at_a = self.evaluate(&self.tight_a);
        let at_b = self.evaluate(&self.tight_b);
        assert_ne!(at_a, Rational::from(0));
        assert_ne!(at_b, Rational::from(0));
        let sign_a = &at_a > &Rational::from(0);
        let sign_b = &at_b > &Rational::from(0);
        if sign_a == sign_b {
            return Err("sign at a and b should be different");
        }
        if self.dir != (sign_a == false) {
            return Err("dir is incorrect");
        }
        Ok(())
    }

    fn new_wide_bounds(poly: Polynomial<Integer>, wide_a: Rational, wide_b: Rational) -> Self {
        let dir = poly.apply_map(|x| Rational::from(x)).evaluate(&wide_a) < Rational::from(0);
        let x = Self {
            poly,
            tight_a: wide_a.clone(),
            tight_b: wide_b.clone(),
            wide_a: LowerBound::Finite(wide_a),
            wide_b: UpperBound::Finite(wide_b),
            dir,
        };
        debug_assert!(x.check_invariants().is_ok());
        x
    }

    fn evaluate(&self, val: &Rational) -> Rational {
        evaluate_at_rational(&self.poly, val)
    }

    pub fn accuracy(&self) -> Rational {
        &self.tight_b - &self.tight_a
    }

    pub fn refine(&mut self) {
        let m = RationalSimpleBetweenGenerator::new_between_middle(
            &self.tight_a,
            &self.tight_b,
            &Rational::from_integers(Integer::from(1), Integer::from(3)),
        )
        .next()
        .unwrap();
        let m_sign = self.evaluate(&m) > Rational::from(0);
        match self.dir == m_sign {
            true => {
                self.tight_b = m;
            }
            false => {
                self.tight_a = m;
            }
        }
    }

    pub fn refine_to_accuracy(&mut self, accuracy: &Rational) {
        while &self.accuracy() > accuracy {
            self.refine();
        }
    }

    pub fn cmp_mut(&mut self, other: &mut Self) -> std::cmp::Ordering {
        let polys_are_eq = self.poly == other.poly; //polys should be irreducible primitive fav-assoc so this is valid
        loop {
            //test for equality: If the tight bounds on one are within the wide bounds of the other
            if polys_are_eq {
                if other.wide_a <= self.tight_a && self.tight_b <= other.wide_b {
                    return std::cmp::Ordering::Equal;
                }
                if self.wide_a <= other.tight_a && other.tight_b <= self.wide_b {
                    return std::cmp::Ordering::Equal;
                }
            }

            //test for inequality: If the tight bounds are disjoint
            if self.tight_b <= other.tight_a {
                return std::cmp::Ordering::Less;
            }
            if other.tight_b <= self.tight_a {
                return std::cmp::Ordering::Greater;
            }

            //refine
            self.refine();
            other.refine();
        }
    }

    pub fn cmp_rat_mut(&mut self, other: &Rational) -> std::cmp::Ordering {
        loop {
            // println!("cmp_rat_mut {:?}", self);

            //test for inequality: other is outside the tight bounds
            if &self.tight_b <= other {
                return std::cmp::Ordering::Less;
            }
            if other <= &self.tight_a {
                return std::cmp::Ordering::Greater;
            }

            // println!("refine");

            //refine
            self.refine();
        }
    }

    fn neg_mut(&mut self) {
        let (unit, fav_assoc) = Polynomial::compose(
            &self.poly,
            &Polynomial::from_coeffs(vec![Integer::from(0), Integer::from(-1)]),
        )
        .factor_fav_assoc();
        if unit == Polynomial::one() {
            self.poly = fav_assoc;
            self.dir = !self.dir;
        } else if unit == Polynomial::neg(&Polynomial::one()) {
            self.poly = fav_assoc;
        } else {
            panic!();
        }
        (self.tight_a, self.tight_b) = (-self.tight_b.clone(), -self.tight_a.clone());
        (self.wide_a, self.wide_b) = (self.wide_b.clone().neg(), self.wide_a.clone().neg());
    }

    fn neg(mut self) -> Self {
        self.neg_mut();
        self
    }

    pub fn min_poly(&self) -> Polynomial<Rational> {
        self.poly.apply_map(|c| Rational::from(c)).fav_assoc()
    }

    pub fn apply_poly(&mut self, poly: &Polynomial<Rational>) -> RealAlgebraic {
        let poly = Polynomial::rem(poly, &self.min_poly());
        match poly.as_constant() {
            Some(rat) => RealAlgebraic::Rational(rat),
            None => {
                let ans_poly = new_anf(self.min_poly())
                    .min_poly(&poly)
                    .primitive_part_fof();

                identify_real_root(
                    ans_poly,
                    (0..).map(|i| {
                        if i != 0 {
                            self.refine();
                        }

                        // eg: c + bx + ax^2 = c + x(b + x(a))
                        let mut coeffs = poly.coeffs().into_iter().rev();
                        let lc = coeffs.next().unwrap();
                        let mut ans = mul_interval_rat((&self.tight_a, &self.tight_b), lc);
                        for (i, c) in coeffs.enumerate() {
                            if i != 0 {
                                ans =
                                    mul_intervals((&ans.0, &ans.1), (&self.tight_a, &self.tight_b));
                            }
                            ans = add_interval_rat((&ans.0, &ans.1), c);
                        }

                        ans
                    }),
                )
            }
        }
    }
}

#[derive(Debug, Clone, Hash)]
pub enum RealAlgebraic {
    Rational(Rational),
    Real(RealAlgebraicRoot),
}

impl RealAlgebraic {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        match self {
            RealAlgebraic::Rational(_x) => {}
            RealAlgebraic::Real(x) => match x.check_invariants() {
                Ok(()) => {}
                Err(e) => {
                    return Err(e);
                }
            },
        }
        Ok(())
    }

    pub fn cmp_mut(&mut self, other: &mut Self) -> std::cmp::Ordering {
        {
            match self {
                RealAlgebraic::Rational(self_rep) => match other {
                    RealAlgebraic::Rational(other_rep) => self_rep.cmp(&other_rep),
                    RealAlgebraic::Real(other_rep) => other_rep.cmp_rat_mut(self_rep).reverse(),
                },
                RealAlgebraic::Real(self_rep) => match other {
                    RealAlgebraic::Rational(other_rep) => self_rep.cmp_rat_mut(other_rep),
                    RealAlgebraic::Real(other_rep) => self_rep.cmp_mut(other_rep),
                },
            }
        }
    }

    pub fn min_poly(&self) -> Polynomial<Rational> {
        match self {
            RealAlgebraic::Rational(rat) => Polynomial::from_coeffs(vec![-rat, Rational::ONE]),
            RealAlgebraic::Real(real_root) => real_root.min_poly(),
        }
    }

    pub fn apply_poly(&mut self, poly: &Polynomial<Rational>) -> Self {
        match self {
            RealAlgebraic::Rational(rat) => RealAlgebraic::Rational(poly.evaluate(rat)),
            RealAlgebraic::Real(real_root) => real_root.apply_poly(poly),
        }
    }

    pub fn degree(&self) -> usize {
        self.min_poly().degree().unwrap()
    }
}

impl PositiveRealNthRootStructure for CannonicalStructure<RealAlgebraic> {
    fn nth_root(&self, x: &Self::Set, n: usize) -> Result<Self::Set, ()> {
        nth_root(x, n)
    }
}

impl PartialEq for RealAlgebraic {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for RealAlgebraic {}

impl PartialOrd for RealAlgebraic {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.clone().cmp_mut(&mut other.clone()))
    }
}

impl Ord for RealAlgebraic {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl MetaType for RealAlgebraic {
    type Structure = CannonicalStructure<RealAlgebraic>;

    fn structure() -> Rc<Self::Structure> {
        CannonicalStructure::new().into()
    }
}

impl Display for RealAlgebraic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RealAlgebraic::Rational(a) => write!(f, "{}", a),
            RealAlgebraic::Real(a) => write!(f, "{}", a),
        }
    }
}

// impl PartialEqStructure for <RealAlgebraic as MetaType>::Structure {
//     fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
//         a == b
//     }
// }

// impl EqStructure for CannonicalStructure<RealAlgebraic> {}

impl RingStructure for CannonicalStructure<RealAlgebraic> {
    fn zero(&self) -> Self::Set {
        RealAlgebraic::Rational(Rational::from(0))
    }

    fn one(&self) -> Self::Set {
        RealAlgebraic::Rational(Rational::from(1))
    }

    fn neg(&self, a: &Self::Set) -> Self::Set {
        match a {
            RealAlgebraic::Rational(a) => RealAlgebraic::Rational(-a),
            RealAlgebraic::Real(root) => RealAlgebraic::Real(root.clone().neg()),
        }
    }

    fn add(&self, alg1: &Self::Set, alg2: &Self::Set) -> Self::Set {
        // println!("add {:?} {:?}", alg1, alg2);

        fn add_rat(mut elem: RealAlgebraicRoot, rat: &Rational) -> RealAlgebraicRoot {
            elem.tight_a += rat;
            elem.tight_b += rat;
            match &elem.wide_a {
                LowerBound::Inf => {}
                LowerBound::Finite(a) => {
                    elem.wide_a = LowerBound::Finite(a + rat);
                }
            }
            match &elem.wide_b {
                UpperBound::Inf => {}
                UpperBound::Finite(b) => {
                    elem.wide_b = UpperBound::Finite(b + rat);
                }
            }

            //compose with x - n/d = dx - n
            elem.poly = Polynomial::compose(
                &elem.poly.apply_map(|c| Rational::from(c)),
                &Polynomial::from_coeffs(vec![-rat, Rational::ONE]),
            )
            .primitive_part_fof();

            #[cfg(debug_assertions)]
            elem.check_invariants().unwrap();
            elem
        }

        match (alg1, alg2) {
            (RealAlgebraic::Rational(rat1), RealAlgebraic::Rational(rat2)) => {
                RealAlgebraic::Rational(rat1 + rat2)
            }
            (RealAlgebraic::Rational(rat1), RealAlgebraic::Real(alg2)) => {
                RealAlgebraic::Real(add_rat(alg2.clone(), rat1))
            }
            (RealAlgebraic::Real(alg1), RealAlgebraic::Rational(rat2)) => {
                RealAlgebraic::Real(add_rat(alg1.clone(), rat2))
            }
            (RealAlgebraic::Real(alg1), RealAlgebraic::Real(alg2)) => {
                let mut alg1 = alg1.clone();
                let mut alg2 = alg2.clone();

                identify_real_root(
                    root_sum_poly(&alg1.poly, &alg2.poly),
                    (0..).map(|i| {
                        if i != 0 {
                            alg1.refine();
                            alg2.refine();
                        }

                        add_intervals(
                            (&alg1.tight_a, &alg1.tight_b),
                            (&alg2.tight_a, &alg2.tight_b),
                        )
                    }),
                )
            }
        }
    }

    fn mul(&self, elem1: &Self::Set, elem2: &Self::Set) -> Self::Set {
        match elem1.cmp(&self.zero()) {
            std::cmp::Ordering::Less => {
                return self.neg(&self.mul(&self.neg(elem1), elem2));
            }
            std::cmp::Ordering::Equal => {
                return self.zero();
            }
            std::cmp::Ordering::Greater => {}
        }

        match elem2.cmp(&self.zero()) {
            std::cmp::Ordering::Less => {
                return self.neg(&self.mul(elem1, &self.neg(elem2)));
            }
            std::cmp::Ordering::Equal => {
                return self.zero();
            }
            std::cmp::Ordering::Greater => {}
        }

        debug_assert!(elem1 > &self.zero());
        debug_assert!(elem2 > &self.zero());

        fn mul_pos_rat(mut elem: RealAlgebraicRoot, rat: &Rational) -> RealAlgebraicRoot {
            debug_assert!(rat > &Rational::from(0));
            elem.tight_a *= rat;
            elem.tight_b *= rat;
            match &elem.wide_a {
                LowerBound::Inf => {}
                LowerBound::Finite(a) => {
                    elem.wide_a = LowerBound::Finite(a * rat);
                }
            }
            match &elem.wide_b {
                UpperBound::Inf => {}
                UpperBound::Finite(b) => {
                    elem.wide_b = UpperBound::Finite(b * rat);
                }
            }
            elem.poly = root_pos_rat_mul_poly(elem.poly, rat);
            #[cfg(debug_assertions)]
            elem.check_invariants().unwrap();
            elem
        }

        match (elem1, elem2) {
            (RealAlgebraic::Rational(rat1), RealAlgebraic::Rational(rat2)) => {
                RealAlgebraic::Rational(rat1 * rat2)
            }
            (RealAlgebraic::Rational(rat1), RealAlgebraic::Real(alg2)) => {
                RealAlgebraic::Real(mul_pos_rat(alg2.clone(), rat1))
            }
            (RealAlgebraic::Real(alg1), RealAlgebraic::Rational(rat2)) => {
                RealAlgebraic::Real(mul_pos_rat(alg1.clone(), rat2))
            }
            (RealAlgebraic::Real(alg1), RealAlgebraic::Real(alg2)) => {
                let mut alg1 = alg1.clone();
                let mut alg2 = alg2.clone();

                identify_real_root(
                    root_prod_poly(&alg1.poly, &alg2.poly),
                    (0..).map(|i| {
                        if i != 0 {
                            alg1.refine();
                            alg2.refine();
                        }
                        let ans_tight_a = &alg1.tight_a * &alg2.tight_a;
                        let ans_tight_b = &alg1.tight_b * &alg2.tight_b;
                        (ans_tight_a, ans_tight_b)
                    }),
                )
            }
        }
    }
}

impl IntegralDomainStructure for CannonicalStructure<RealAlgebraic> {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        match self.inv(b) {
            Ok(b_inv) => Ok(self.mul(a, &b_inv)),
            Err(err) => Err(err),
        }
    }

    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        let mut a = a.clone();
        match RealAlgebraic::cmp_mut(&mut a, &mut self.zero()) {
            std::cmp::Ordering::Less => match self.inv(&self.neg(&a)) {
                Ok(neg_elem_inv) => Ok(self.neg(&neg_elem_inv)),
                Err(err) => Err(err),
            },
            std::cmp::Ordering::Equal => Err(RingDivisionError::DivideByZero),
            std::cmp::Ordering::Greater => match a {
                RealAlgebraic::Rational(x) => {
                    Ok(RealAlgebraic::Rational(Rational::inv(&x).unwrap()))
                }
                RealAlgebraic::Real(mut root) => {
                    debug_assert!(root.tight_a >= Rational::from(0));
                    while root.tight_a == Rational::from(0) {
                        root.refine();
                    }
                    debug_assert!(Rational::from(0) < root.tight_a);
                    (root.tight_a, root.tight_b) = (
                        Rational::inv(&root.tight_b).unwrap(),
                        Rational::inv(&root.tight_a).unwrap(),
                    );
                    (root.wide_a, root.wide_b) = (
                        {
                            match root.wide_b {
                                UpperBound::Inf => LowerBound::Finite(Rational::from(0)),
                                UpperBound::Finite(x) => match Rational::inv(&x) {
                                    Ok(x_inv) => LowerBound::Finite(x_inv),
                                    Err(RingDivisionError::DivideByZero) => panic!("wide upper bound of strictly positive root should be strictly positive i.e. non-zero"),
                                    Err(_) => panic!()
                                },
                            }
                        },
                        {
                            match root.wide_a {
                                LowerBound::Inf => UpperBound::Inf,
                                LowerBound::Finite(x) => match x.cmp(&Rational::from(0)) {
                                    std::cmp::Ordering::Less => UpperBound::Inf,
                                    std::cmp::Ordering::Equal => UpperBound::Inf,
                                    std::cmp::Ordering::Greater => {
                                        UpperBound::Finite(Rational::inv(&x).unwrap())
                                    }
                                },
                            }
                        },
                    );
                    let (unit, fav_assoc) = Polynomial::from_coeffs(
                        root.poly.into_coeffs().into_iter().rev().collect(),
                    )
                    .factor_fav_assoc();
                    if unit == Polynomial::one() {
                        root.poly = fav_assoc;
                        root.dir = !root.dir;
                    } else if unit == Polynomial::neg(&Polynomial::one()) {
                        root.poly = fav_assoc;
                    } else {
                        panic!();
                    }
                    let ans = RealAlgebraic::Real(root);
                    #[cfg(debug_assertions)]
                    ans.check_invariants().unwrap();
                    Ok(ans)
                }
            },
        }
    }
}

impl FieldStructure for CannonicalStructure<RealAlgebraic> {}

impl CharZeroStructure for CannonicalStructure<RealAlgebraic> {}

impl ComplexSubsetStructure for CannonicalStructure<RealAlgebraic> {}

impl RealSubsetStructure for CannonicalStructure<RealAlgebraic> {}

impl OrderedRingStructure for CannonicalStructure<RealAlgebraic> {
    fn ring_cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
        a.cmp(b)
    }
}

impl RealToFloatStructure for CannonicalStructure<RealAlgebraic> {
    fn as_f64(&self, x: &Self::Set) -> f64 {
        match x {
            RealAlgebraic::Rational(x) => x.as_f64(),
            RealAlgebraic::Real(x) => {
                let mut x = x.clone();
                x.refine_to_accuracy(&Rational::from_integers(
                    Integer::from(1),
                    Integer::from(1000000000000000i64),
                ));
                ((x.tight_a + x.tight_b) / Rational::from(2)).as_f64()
            }
        }
    }
}

impl RealRoundingStructure for CannonicalStructure<RealAlgebraic> {
    fn floor(&self, _x: &Self::Set) -> Integer {
        todo!()
    }
    fn ceil(&self, _x: &Self::Set) -> Integer {
        todo!()
    }
    fn round(&self, _x: &Self::Set) -> Integer {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    

    use crate::ring_structure::elements::IntoErgonomic;

    use super::*;

    #[test]
    fn test_real_neg() {
        {
            let f = Polynomial::from_coeffs(vec![
                Integer::from(-2),
                Integer::from(0),
                Integer::from(1),
            ]);
            let roots = Polynomial::all_real_roots(&f);

            assert_eq!(roots.len(), 2);
            let a = &roots[0];
            let b = &roots[1];

            let a_neg = RealAlgebraic::neg(a);
            let b_neg = RealAlgebraic::neg(b);

            a_neg.check_invariants().unwrap();
            b_neg.check_invariants().unwrap();

            println!("a = {}", a.to_string());
            println!("b = {}", b.to_string());
            println!("a_neg = {}", a_neg.to_string());
            println!("b_neg = {}", b_neg.to_string());

            assert_ne!(a, b);
            assert_eq!(a, &b_neg);
            assert_eq!(b, &a_neg);
        }
        {
            let f = Polynomial::from_coeffs(vec![
                Integer::from(-1),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(3),
                Integer::from(1),
            ]);
            let roots = Polynomial::all_real_roots(&f);

            assert_eq!(roots.len(), 3);
            for root in roots {
                RealAlgebraic::neg(&root).check_invariants().unwrap();
            }
        }
        {
            //example where f(g(x)) is not primitive even though f and g are
            let f = Polynomial::from_coeffs(vec![
                Integer::from(-4),
                Integer::from(-1),
                Integer::from(1),
            ]);
            let roots = Polynomial::all_real_roots(&f);
            for root in roots {
                let root2 = RealAlgebraic::add(
                    &root,
                    &RealAlgebraic::from_rat(&Rational::from_signeds(1, 2)).unwrap(),
                );
                root2.check_invariants().unwrap();
            }
        }
    }

    #[test]
    fn test_real_add() {
        let f =
            Polynomial::from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(3)]);
        let roots = Polynomial::all_real_roots(&f);
        let a = RealAlgebraic::sum(roots.iter().collect());
        assert_eq!(a, RealAlgebraic::zero());

        let f = Polynomial::from_coeffs(vec![
            Integer::from(-7),
            Integer::from(0),
            Integer::from(100),
        ]);
        let roots = Polynomial::all_real_roots(&f);
        let a = RealAlgebraic::sum(roots.iter().collect());
        assert_eq!(a, RealAlgebraic::zero());

        let f = Polynomial::from_coeffs(vec![
            Integer::from(-100),
            Integer::from(0),
            Integer::from(7),
        ]);
        let roots = Polynomial::all_real_roots(&f);
        let a = RealAlgebraic::sum(roots.iter().collect());
        assert_eq!(a, RealAlgebraic::zero());
    }

    #[test]
    fn test_real_mul() {
        let f = Polynomial::from_coeffs(vec![
            Integer::from(-100),
            Integer::from(0),
            Integer::from(7),
        ]);
        // (x-a)(x-b) = x^2 - 100/7
        // so ab=-100/7
        let roots = Polynomial::all_real_roots(&f);
        let a = RealAlgebraic::product(roots.iter().collect());
        assert_eq!(
            a,
            RealAlgebraic::from_rat(&Rational::from_signeds(-100, 7)).unwrap()
        );
    }

    #[test]
    fn test_real_nth_root() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();
        let f = ((4 * x.pow(5) - 12 * x.pow(3) + 8 * x + 1)
            * (x + 1)
            * (x)
            * (x - 1)
            * (x - 2)
            * (x - 3)
            * (x - 4)
            * (x - 5)
            * (x - 144)
            * (x.pow(2) - 3))
            .into_verbose();
        let n = 2;

        for root in f.all_real_roots() {
            println!();
            println!("root = {}", root);
            match root.nth_root(n) {
                Ok(nth_root) => {
                    println!("YES {}-root = {}", n, nth_root);
                    debug_assert!(RealAlgebraic::zero() <= root);
                    debug_assert!(RealAlgebraic::zero() <= nth_root);
                    debug_assert_eq!(nth_root.nat_pow(&Natural::from(n)), root);
                }
                Err(()) => {
                    println!("NO {}-root", n);
                    debug_assert!(RealAlgebraic::zero() > root);
                }
            }
        }
    }

    #[test]
    fn test_real_algebraic_ordering() {
        let mut all_roots = vec![];
        for f in vec![
            Polynomial::from_coeffs(vec![
                Integer::from(-2),
                Integer::from(-4),
                Integer::from(-2),
            ]),
            Polynomial::from_coeffs(vec![Integer::from(6), Integer::from(0), Integer::from(-3)]),
            Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-3), Integer::from(1)]),
            Polynomial::from_coeffs(vec![
                Integer::from(2),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ]),
            Polynomial::from_coeffs(vec![
                Integer::from(1),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ]),
            Polynomial::from_coeffs(vec![
                Integer::from(-1),
                Integer::from(12),
                Integer::from(-4),
                Integer::from(-15),
                Integer::from(5),
                Integer::from(3),
                Integer::from(-1),
            ]),
        ] {
            for root in Polynomial::real_roots(&f, None, None, false, false) {
                all_roots.push(root.clone());
            }
        }

        all_roots.sort();

        for mut root in &mut all_roots {
            root.check_invariants().unwrap();
            match &mut root {
                RealAlgebraic::Rational(_a) => {}
                RealAlgebraic::Real(a) => {
                    a.refine_to_accuracy(&Rational::from_signeds(1, i64::MAX))
                }
            }
            println!("    {} {:?}", root.to_string(), root);
        }

        let mut all_roots_sorted_by_lower_tight_bound = all_roots.clone();
        all_roots_sorted_by_lower_tight_bound.sort_by_key(|root| match root {
            RealAlgebraic::Rational(a) => a.clone(),
            RealAlgebraic::Real(r) => r.tight_a.clone(),
        });
        assert_eq!(all_roots, all_roots_sorted_by_lower_tight_bound);
    }
}
