#![allow(dead_code)]

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;
use malachite_q::Rational;

use super::multipoly::*;
use super::nzq::*;
use super::poly::*;
use super::ring::ComRing;
use super::ring::*;

pub const QQ_BAR_REAL: RealAlgebraicField = RealAlgebraicField {};
pub const QQ_BAR: ComplexAlgebraicField = ComplexAlgebraicField {};

fn root_sum_poly(p: &Polynomial<Integer>, q: &Polynomial<Integer>) -> Polynomial<Integer> {
    let zz_poly_of_multi_poly = PolynomialRing::new(&ZZ_MULTIPOLY);

    let x = Variable::new(String::from("x"));
    let z = Variable::new(String::from("z"));

    let p = ZZ_POLY.apply_map(&ZZ_MULTIPOLY, p, |c| ZZ_MULTIPOLY.constant(c.clone()));
    let q = ZZ_POLY.apply_map(&ZZ_MULTIPOLY, q, |c| ZZ_MULTIPOLY.constant(c.clone()));
    let r = ZZ_MULTIPOLY.expand(
        &zz_poly_of_multi_poly.evaluate(
            &q,
            &ZZ_MULTIPOLY.add(
                ZZ_MULTIPOLY.var(z.clone()),
                ZZ_MULTIPOLY.neg(ZZ_MULTIPOLY.var(x.clone())),
            ),
        ),
        &x,
    );

    let root_sum_poly = zz_poly_of_multi_poly.apply_map(
        &ZZ,
        &ZZ_MULTIPOLY.expand(&zz_poly_of_multi_poly.resultant(p, r), &z),
        |c| ZZ_MULTIPOLY.as_constant(c).unwrap(),
    );
    if root_sum_poly == ZZ_POLY.zero() {
        root_sum_poly
    } else {
        ZZ_POLY.primitive_squarefree_part(root_sum_poly).unwrap()
    }
}

fn root_prod_poly(p: &Polynomial<Integer>, q: &Polynomial<Integer>) -> Polynomial<Integer> {
    let zz_poly_of_multi_poly = PolynomialRing::new(&ZZ_MULTIPOLY);

    let x = Variable::new(String::from("x"));
    let t = Variable::new(String::from("t"));

    let p = ZZ_POLY.apply_map(&ZZ_MULTIPOLY, p, |c| ZZ_MULTIPOLY.constant(c.clone()));
    let q = zz_poly_of_multi_poly.evaluate(
        &ZZ_POLY.apply_map(&ZZ_MULTIPOLY, q, |c| ZZ_MULTIPOLY.constant(c.clone())),
        &ZZ_MULTIPOLY.var(x.clone()),
    );
    let r = ZZ_MULTIPOLY.expand(&ZZ_MULTIPOLY.homogenize(&q, &t), &t);
    //x ** q.degree() * q(t * x ** -1)

    let root_prod_poly = zz_poly_of_multi_poly.apply_map(
        &ZZ,
        &ZZ_MULTIPOLY.expand(&zz_poly_of_multi_poly.resultant(p, r), &x),
        |c| ZZ_MULTIPOLY.as_constant(c).unwrap(),
    );
    if root_prod_poly == ZZ_POLY.zero() {
        root_prod_poly
    } else {
        ZZ_POLY.primitive_squarefree_part(root_prod_poly).unwrap()
    }
}

impl<'a> PolynomialRing<'a, IntegerRing> {
    fn sign_variations(&self, poly: &Polynomial<Integer>) -> usize {
        //https://en.wikipedia.org/wiki/Descartes'_rule_of_signs
        //equals the number of strictly positive real roots modulo 2
        //and number of positive real roots is less than this number
        let nonzero_coeffs: Vec<_> = poly
            .coeffs()
            .into_iter()
            .filter(|c| c != &self.ring().zero())
            .collect();
        let mut v = 0;
        for i in 0..nonzero_coeffs.len() - 1 {
            if (nonzero_coeffs[i] < 0) != (nonzero_coeffs[i + 1] < 0) {
                v += 1;
            }
        }
        v
    }

    //Collins and Akritas algorithm %https://en.wikipedia.org/wiki/Real-root_isolation
    fn isolate_real_roots_by_collin_akritas(
        &self,
        poly: &Polynomial<Integer>,
    ) -> Vec<(Natural, usize, bool)> {
        //input: p(x), a square-free polynomial, such that p(0) p(1) ≠ 0, for which the roots in the interval [0, 1] are searched
        //output: a list of triples (c, k, h) representing isolating intervals of the form [c/2^k, (c+h)/2^k]
        debug_assert_ne!(self.evaluate(poly, &self.ring().zero()), self.ring().zero());
        debug_assert_ne!(self.evaluate(poly, &self.ring().one()), self.ring().zero());
        debug_assert_eq!(
            self.degree(&self.primitive_squarefree_part(poly.clone()).unwrap())
                .unwrap(),
            self.degree(poly).unwrap()
        );

        let mut n = self.degree(&poly).unwrap();
        let mut l = vec![(Natural::from(0u8), 0, poly.clone())];
        let mut isol = vec![];
        while l.len() != 0 {
            let (c, k, mut q) = l.pop().unwrap();
            if ZZ_POLY.evaluate(&q, &Integer::from(0)) == Integer::from(0) {
                //q = q/x
                q = self.div(q, self.var()).unwrap();
                n -= 1;
                isol.push((c.clone(), k.clone(), false)); //rational root
            }
            let v = self.sign_variations(&self.compose(
                &self.reversed(&q),
                &self.from_coeffs(vec![Integer::from(1), Integer::from(1)]),
            ));
            if v == 1 {
                isol.push((c, k, true)); //root
            } else if v > 1 {
                //bisect
                //q_small(x) = 2^n q(x/2)
                let q_small = self.apply_map_with_powers(&ZZ, &q, |(i, coeff)| {
                    coeff * Integer::from(2) << (n - i)
                });
                l.push((c.clone() << 1, k + 1, q_small.clone()));
                l.push((
                    (c << 1) + Natural::from(1u8),
                    k + 1,
                    self.compose(
                        &q_small,
                        &self.from_coeffs(vec![Integer::from(1), Integer::from(1)]),
                    ),
                ));
            }
        }
        isol
    }

    //isolate all real roots of the irreducible poly in the open interval (a, b)
    fn real_roots_irreducible(
        &self,
        poly: &Polynomial<Integer>,
        a: &Option<Rational>,
        b: &Option<Rational>,
    ) -> Vec<RealAlgebraicNumber> {
        assert_ne!(poly, &self.zero());
        debug_assert!(self.is_irreducible(&poly).unwrap());

        let (mut a, mut b) = (a.clone(), b.clone());

        let d = ZZ_POLY.degree(&poly).unwrap();
        if d == 0 {
            //constant polynomial has no roots
            vec![]
        } else if d == 1 {
            //poly = a+bx
            //root = -a/b
            let root = -Rational::from(self.coeff(&poly, 0)) / Rational::from(self.coeff(&poly, 1));
            if a.is_some() {
                if !(a.unwrap() < root) {
                    return vec![];
                }
            }
            if b.is_some() {
                if !(root < b.unwrap()) {
                    return vec![];
                }
            }
            vec![RealAlgebraicNumber::Rational(root)]
        } else {
            //compute a bound M on the absolute value of any root
            if a.is_none() || b.is_none() {
                //m = (Cauchy's bound + 1) https://captainblack.wordpress.com/2009/03/08/cauchys-upper-bound-for-the-roots-of-a-polynomial/
                let m = Rational::from(2)
                    + Rational::from_integers(
                        Integer::from(
                            itertools::max(
                                (0..d).map(|i| ZZ_POLY.coeff(&poly, i).unsigned_abs_ref().clone()),
                            )
                            .unwrap(),
                        ),
                        ZZ_POLY.coeff(&poly, d),
                    );
                if a.is_none() {
                    a = Some(-m.clone());
                }
                if b.is_none() {
                    b = Some(m.clone());
                }
            }
            let (a, b) = (a.unwrap(), b.unwrap());
            assert!(a < b);

            //there are not roots A < r < B if A >= B
            if a >= b {
                return vec![];
            }

            //apply a transformation to p so that its roots in (a, b) are moved to roots in (0, 1)
            let (_, trans_poly) = QQ_POLY.factor_primitive_fof(&QQ_POLY.compose(
                &ZZ_POLY.apply_map(&QQ, &poly, |c| Rational::from(c)),
                &QQ_POLY.from_coeffs(vec![a.clone(), b.clone() - a.clone()]),
            ));

            let mut roots = vec![];
            for (c, k, h) in self.isolate_real_roots_by_collin_akritas(&trans_poly) {
                assert!(h); //should not isolate any rational roots since poly is irreducible with degree >= 2
                let d = Natural::from(1u8) << k;
                roots.push(RealAlgebraicNumber::Real(
                    RealAlgebraicRoot::new_wide_bounds(
                        poly.clone(),
                        (&b - &a) * Rational::from_naturals(c.clone(), d.clone()) + &a,
                        (&b - &a) * Rational::from_naturals(&c + Natural::from(1u8), d.clone())
                            + &a,
                    ),
                ));
            }
            roots
        }
    }

    //get the real roots with multiplicity of poly
    pub fn real_roots(
        &self,
        poly: &Polynomial<Integer>,
        a: &Option<Rational>,
        b: &Option<Rational>,
    ) -> Vec<RealAlgebraicNumber> {
        assert_ne!(poly, &self.zero());
        let factors = self.factor(&poly).unwrap();
        let mut roots = vec![];
        for (factor, k) in factors.factors() {
            for root in self.real_roots_irreducible(factor, a, b) {
                let mut i = Natural::from(0u8);
                while &i < k {
                    roots.push(root.clone());
                    i += Natural::from(1u8);
                }
            }
        }
        roots
    }
}

#[derive(Debug, Clone, Hash)]
pub struct RealAlgebraicRoot {
    poly: Polynomial<Integer>, //a primitive irreducible polynomial of degree >= 2 with a unique real root between a and b
    //an arbitrarily small interval containing the root. May be mutated
    tight_a: Rational, //tight lower bound
    tight_b: Rational, //tight upper bound
    //a heuristically large interval containing the root. Should not shrink
    wide_a: Rational, //wide lower bound
    wide_b: Rational, //wide upper bound
    //false : decreasing i.e. poly(a) > poly(b), true : increasing i.e. poly(a) < poly(b)
    dir: bool,
}

impl RealAlgebraicRoot {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        if !(self.tight_a < self.tight_b) {
            return Err("tight a should be strictly less than b");
        }
        if !(self.wide_a < self.wide_b) {
            return Err("wide a should be strictly less than b");
        }
        if self.poly
            != ZZ_POLY
                .factor_fav_assoc(
                    ZZ_POLY
                        .primitive_squarefree_part(self.poly.clone())
                        .unwrap(),
                )
                .1
        {
            return Err("poly should be primitive and favoriate associate");
        }
        match ZZ_POLY.is_irreducible(&self.poly) {
            Some(is_irr) => {
                if !is_irr {
                    return Err("poly should be irreducible");
                }
            }
            None => {
                return Err("poly should be non-zero");
            }
        }
        if ZZ_POLY.degree(&self.poly).unwrap() < 2 {
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
        let dir = QQ_POLY.evaluate(
            &ZZ_POLY.apply_map(&QQ, &poly, |x| Rational::from(x)),
            &wide_a,
        ) < Rational::from(0);
        let x = Self {
            poly,
            tight_a: wide_a.clone(),
            tight_b: wide_b.clone(),
            wide_a,
            wide_b,
            dir,
        };
        debug_assert!(x.check_invariants().is_ok());
        x
    }

    fn evaluate(&self, val: &Rational) -> Rational {
        QQ_POLY.evaluate(
            &ZZ_POLY.apply_map(&QQ, &self.poly, |x| Rational::from(x)),
            &val,
        )
    }

    pub fn accuracy(&self) -> Rational {
        &self.tight_b - &self.tight_a
    }

    pub fn refine(&mut self) {
        let m = (&self.tight_a + &self.tight_b) / Rational::from(2);
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
        let eq_poly = self.poly == other.poly; //polys should be irreducible primitive fav-assoc so this is valid
        loop {
            //test for equality: if the tight bounds on one are within the wide bounds of the other
            if eq_poly {
                if other.wide_a <= self.tight_a && self.tight_b <= other.wide_b {
                    return std::cmp::Ordering::Equal;
                }
                if self.wide_a <= other.tight_a && other.tight_b <= self.wide_b {
                    return std::cmp::Ordering::Equal;
                }
            }

            //test for inequality: if the tight bounds are disjoint
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
            //test for inequality: other is outside the tight bounds
            if &self.tight_b <= other {
                return std::cmp::Ordering::Less;
            }
            if other <= &self.tight_a {
                return std::cmp::Ordering::Greater;
            }

            //refine
            self.refine();
        }
    }
}

impl ToString for RealAlgebraicRoot {
    fn to_string(&self) -> String {
        let m = (&self.tight_a + &self.tight_b) / Rational::from(2);

        fn rat_to_string(a: Rational) -> String {
            let neg = a < Rational::from(0);
            let (mant, exp): (f64, _) = a
                .sci_mantissa_and_exponent_with_rounding(
                    malachite_base::rounding_modes::RoundingMode::Nearest,
                )
                .unwrap();
            let mut b = (2.0 as f64).powf(exp as f64) * mant;
            if neg {
                b = -b;
            }
            b.to_string()
        }

        "≈".to_owned()
            + rat_to_string(m).as_str()
            + "±"
            + rat_to_string(self.accuracy() / Rational::from(2)).as_str()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ComplexAlgebraicRoot {
    tight_a: Rational, //tight lower bound for the real part
    tight_b: Rational, //tight upper bound for the real part
    tight_c: Rational, //tight lower bound for the imaginary part
    tight_d: Rational, //tight upper bound for the imaginary part

    wide_a: Rational, //wide lower bound for the real part
    wide_b: Rational, //wide upper bound for the real part
    wide_c: Rational, //wide lower bound for the imaginary part
    wide_d: Rational, //wide upper bound for the imaginary part

    poly: Polynomial<Integer>, //a primitive irreducible polynomial of degree >= 2 with a unique non-real complex root in the box defined by (a, b, c, d)
}

impl ComplexAlgebraicRoot {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        if !(self.tight_a < self.tight_b) {
            return Err("tight a should be strictly less than b");
        }
        if !(self.tight_c < self.tight_d) {
            return Err("tight c should be strictly less than d");
        }
        if !(self.wide_a < self.wide_b) {
            return Err("wide a should be strictly less than b");
        }
        if !(self.wide_c < self.wide_d) {
            return Err("wide c should be strictly less than d");
        }
        match ZZ_POLY.is_irreducible(&self.poly) {
            Some(is_irr) => {
                if !is_irr {
                    return Err("poly should be irreducible");
                }
            }
            None => {
                return Err("poly should be non-zero");
            }
        }
        if ZZ_POLY.degree(&self.poly).unwrap() < 2 {
            return Err("poly should have degree at least 2");
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Hash)]
pub enum RealAlgebraicNumber {
    Rational(Rational),
    Real(RealAlgebraicRoot),
}

impl RealAlgebraicNumber {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        match self {
            RealAlgebraicNumber::Rational(x) => {}
            RealAlgebraicNumber::Real(x) => match x.check_invariants() {
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
                RealAlgebraicNumber::Rational(self_rep) => match other {
                    RealAlgebraicNumber::Rational(other_rep) => self_rep.cmp(&other_rep),
                    RealAlgebraicNumber::Real(other_rep) => {
                        other_rep.cmp_rat_mut(self_rep).reverse()
                    }
                },
                RealAlgebraicNumber::Real(self_rep) => match other {
                    RealAlgebraicNumber::Rational(other_rep) => {
                        self_rep.cmp_rat_mut(other_rep)
                    }
                    RealAlgebraicNumber::Real(other_rep) => self_rep.cmp_mut(other_rep),
                },
            }
        }
    }
}

impl PartialEq for RealAlgebraicNumber {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for RealAlgebraicNumber {}

impl PartialOrd for RealAlgebraicNumber {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.clone().cmp_mut(&mut other.clone()))
    }
}

impl Ord for RealAlgebraicNumber {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ComplexAlgebraicNumber {
    Real(RealAlgebraicNumber),
    Complex(ComplexAlgebraicRoot),
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RealAlgebraicField {}

impl ComRing for RealAlgebraicField {
    type ElemT = RealAlgebraicNumber;

    fn to_string(&self, elem: &Self::ElemT) -> String {
        match elem {
            RealAlgebraicNumber::Rational(a) => a.to_string(),
            RealAlgebraicNumber::Real(a) => a.to_string(),
        }
    }

    fn zero(&self) -> Self::ElemT {
        RealAlgebraicNumber::Rational(Rational::from(0))
    }

    fn one(&self) -> Self::ElemT {
        RealAlgebraicNumber::Rational(Rational::from(1))
    }

    fn neg_mut(&self, elem: &mut Self::ElemT) {
        todo!()
    }

    fn add_mut(&self, elem: &mut Self::ElemT, offset: &Self::ElemT) {
        todo!()
    }

    fn mul_mut(&self, elem: &mut Self::ElemT, mul: &Self::ElemT) {
        todo!()
    }

    fn div(&self, a: Self::ElemT, b: Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        todo!()
    }
}

impl IntegralDomain for RealAlgebraicField {}

impl Field for RealAlgebraicField {}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ComplexAlgebraicField {}

impl ComRing for ComplexAlgebraicField {
    type ElemT = ComplexAlgebraicNumber;

    fn to_string(&self, elem: &Self::ElemT) -> String {
        match elem {
            ComplexAlgebraicNumber::Real(a) => RealAlgebraicField {}.to_string(a),
            ComplexAlgebraicNumber::Complex(a) => todo!(),
        }
    }

    fn zero(&self) -> Self::ElemT {
        ComplexAlgebraicNumber::Real(RealAlgebraicNumber::Rational(Rational::from(0)))
    }

    fn one(&self) -> Self::ElemT {
        ComplexAlgebraicNumber::Real(RealAlgebraicNumber::Rational(Rational::from(1)))
    }

    fn neg_mut(&self, elem: &mut Self::ElemT) {
        todo!()
    }

    fn add_mut(&self, elem: &mut Self::ElemT, offset: &Self::ElemT) {
        todo!()
    }

    fn mul_mut(&self, elem: &mut Self::ElemT, mul: &Self::ElemT) {
        todo!()
    }

    fn div(&self, a: Self::ElemT, b: Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        todo!()
    }
}

impl IntegralDomain for ComplexAlgebraicField {}

impl Field for ComplexAlgebraicField {}

#[cfg(test)]
mod tests {
    use super::super::poly::*;
    use super::*;

    #[test]
    fn test_root_sum_poly() {
        for (f, g, exp) in vec![
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(0)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(0)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(0)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(0)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(-3), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-5), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-8), Integer::from(1)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-7), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(-1), Integer::from(2)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-1), Integer::from(3)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-5), Integer::from(6)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(-1), Integer::from(-2), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![
                    Integer::from(-7),
                    Integer::from(5),
                    Integer::from(3),
                    Integer::from(-1),
                ]),
            ),
        ] {
            println!();
            let rsp = root_sum_poly(&f, &g);
            println!("f = {}", ZZ_POLY.to_string(&f));
            println!("g = {}", ZZ_POLY.to_string(&g));
            println!(
                "exp = {}    exp_factored = {:?}",
                ZZ_POLY.to_string(&exp),
                ZZ_POLY.factorize_by_kroneckers_method(&exp)
            );
            println!(
                "rsp = {}    rsp_factored = {:?}",
                ZZ_POLY.to_string(&rsp),
                ZZ_POLY.factorize_by_kroneckers_method(&rsp)
            );
            assert!(ZZ_POLY.are_associate(&exp, &rsp));
        }
    }

    #[test]
    fn test_root_prod_poly() {
        for (f, g, exp) in vec![
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(0)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(0)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(0)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(0)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(-3), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-5), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-15), Integer::from(1)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-7), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(1)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(-1), Integer::from(2)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-1), Integer::from(3)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-1), Integer::from(6)]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(-1), Integer::from(-2), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![
                    Integer::from(4),
                    Integer::from(0),
                    Integer::from(-12),
                    Integer::from(0),
                    Integer::from(1),
                ]),
            ),
            (
                ZZ_POLY.from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-2), Integer::from(0), Integer::from(1)]),
                ZZ_POLY.from_coeffs(vec![Integer::from(-4), Integer::from(0), Integer::from(1)]),
            ),
        ] {
            println!();
            let rpp = root_prod_poly(&f, &g);
            println!("f = {}", ZZ_POLY.to_string(&f));
            println!("g = {}", ZZ_POLY.to_string(&g));
            println!(
                "exp = {}    exp_factored = {:?}",
                ZZ_POLY.to_string(&exp),
                ZZ_POLY.factorize_by_kroneckers_method(&exp)
            );
            println!(
                "rpp = {}    rpp_factored = {:?}",
                ZZ_POLY.to_string(&rpp),
                ZZ_POLY.factorize_by_kroneckers_method(&rpp)
            );
            assert!(ZZ_POLY.are_associate(&exp, &rpp));
        }
    }

    #[test]
    fn test_real_root_irreducible_count() {
        assert_eq!(
            ZZ_POLY
                .real_roots_irreducible(
                    &ZZ_POLY.from_coeffs(vec![
                        Integer::from(3),
                        Integer::from(-3),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(1)
                    ]),
                    &None,
                    &None
                )
                .len(),
            1
        );
        assert_eq!(
            ZZ_POLY
                .real_roots_irreducible(
                    &ZZ_POLY.from_coeffs(vec![
                        Integer::from(1),
                        Integer::from(-3),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(1)
                    ]),
                    &None,
                    &None
                )
                .len(),
            3
        );
    }
}
