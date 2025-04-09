use super::bisection_gen::RationalSimpleBetweenGenerator;
use super::poly_tools::*;
use super::real::RealAlgebraic;
use crate::polynomial::*;
use crate::structure::*;
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use boxes::*;
use std::{collections::HashSet, fmt::Display, str::FromStr};

mod boxes;
pub mod polynomial;

fn bisect_box(
    poly: &Polynomial<Integer>,
    n: usize,
    a: &Rational,
    b: &Rational,
    c: &Rational,
    d: &Rational,
) -> (
    (usize, Rational, Rational, Rational, Rational),
    (usize, Rational, Rational, Rational, Rational),
) {
    let ba = b - a;
    let dc = d - c;
    if ba >= dc {
        //bisect (a, b)
        for m in RationalSimpleBetweenGenerator::new_between_middle(
            a,
            b,
            &Rational::from_integers(Integer::from(1), Integer::from(3)),
        ) {
            match (
                poly.count_complex_roots(&a, &m, &c, &d),
                poly.count_complex_roots(&m, &b, &c, &d),
            ) {
                (Some(n1), Some(n2)) => {
                    debug_assert_eq!(n1 + n2, n);
                    return (
                        (n1, a.clone(), m.clone(), c.clone(), d.clone()),
                        (n2, m, b.clone(), c.clone(), d.clone()),
                    );
                }
                _ => {}
            }
        }
    } else {
        //bisect (c, d)
        for m in RationalSimpleBetweenGenerator::new_between_middle(
            c,
            d,
            &Rational::from_integers(Integer::from(1), Integer::from(3)),
        ) {
            match (
                poly.count_complex_roots(&a, &b, &c, &m),
                poly.count_complex_roots(&a, &b, &m, &d),
            ) {
                (Some(n1), Some(n2)) => {
                    debug_assert_eq!(n1 + n2, n);
                    return (
                        (n1, a.clone(), b.clone(), c.clone(), m.clone()),
                        (n2, a.clone(), b.clone(), m, d.clone()),
                    );
                }
                _ => {}
            }
        }
    }
    unreachable!()
}

//box gen must yield boxes which converge to some root of the polynomial
fn identify_complex_root(
    poly: Polynomial<Integer>,
    mut box_gen: impl Iterator<Item = (Rational, Rational, Rational, Rational)>,
) -> ComplexAlgebraic {
    let poly = poly.primitive_squarefree_part();
    //find the irreducible factor poly which contains the root being converged to
    let (mut a, mut b, mut c, mut d) = box_gen.next().unwrap();

    let irr_poly = {
        let (_unit, factors) = poly.factor().unwrap().unit_and_factors();
        let irr_polys = factors.into_iter().map(|(f, _k)| f).collect::<Vec<_>>();
        let mut possible_irr_poly_idxs: HashSet<_> = (0..irr_polys.len()).collect();
        loop {
            debug_assert!(!possible_irr_poly_idxs.is_empty());
            possible_irr_poly_idxs = possible_irr_poly_idxs
                .into_iter()
                .filter(|idx| irr_polys[*idx].count_complex_roots(&a, &b, &c, &d) != Some(0))
                .collect();
            if possible_irr_poly_idxs.len() == 1 {
                break;
            }
            (a, b, c, d) = box_gen.next().unwrap();
        }
        debug_assert_eq!(possible_irr_poly_idxs.len(), 1);
        irr_polys
            .into_iter()
            .nth(possible_irr_poly_idxs.into_iter().next().unwrap())
            .unwrap()
    };

    let mut roots = irr_poly.all_complex_roots_irreducible();
    let mut possible_roots: HashSet<_> = (0..roots.len()).collect();
    loop {
        debug_assert!(!possible_roots.is_empty());
        possible_roots = possible_roots
            .into_iter()
            .filter(|idx| match &roots[*idx] {
                ComplexAlgebraic::Real(RealAlgebraic::Rational(root)) => {
                    &a < root && root < &b && c < Rational::ZERO && Rational::ZERO < d
                }
                ComplexAlgebraic::Real(RealAlgebraic::Real(root)) => {
                    &a < root.tight_b()
                        && root.tight_a() < &b
                        && c < Rational::ZERO
                        && Rational::ZERO < d
                }
                ComplexAlgebraic::Complex(root) => {
                    a < root.tight_b && root.tight_a < b && c < root.tight_d && root.tight_c < d
                }
            })
            .collect();
        if possible_roots.len() == 1 {
            break;
        }
        (a, b, c, d) = box_gen.next().unwrap();
        for idx in &possible_roots {
            match &mut roots[*idx] {
                ComplexAlgebraic::Real(RealAlgebraic::Rational(_root)) => {}
                ComplexAlgebraic::Real(RealAlgebraic::Real(root)) => {
                    root.refine();
                }
                ComplexAlgebraic::Complex(root) => {
                    root.refine();
                }
            }
        }
    }
    debug_assert_eq!(possible_roots.len(), 1);
    let ans = roots
        .into_iter()
        .nth(possible_roots.into_iter().next().unwrap())
        .unwrap();
    #[cfg(debug_assertions)]
    ans.check_invariants().unwrap();
    ans
}

#[derive(Debug, Clone)]
pub struct ComplexAlgebraicRoot {
    tight_a: Rational, //tight lower bound for the real part
    tight_b: Rational, //tight upper bound for the real part
    tight_c: Rational, //tight lower bound for the imaginary part
    tight_d: Rational, //tight upper bound for the imaginary part

    poly: Polynomial<Integer>, //a primitive irreducible polynomial of degree >= 2 with a unique non-real complex root in the box defined by (a, b, c, d)
}

impl Display for ComplexAlgebraicRoot {
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
            debug_assert_eq!(d_sign, -Integer::ONE); //because we are a real number
            debug_assert_eq!(-d, &d_sqfreee * &d_sq * &d_sq);

            let two_a = Integer::TWO * a;

            let x = Rational::from_integers(-b, two_a.clone());
            let y = Rational::from_integers(d_sq, two_a);
            debug_assert!(y > Rational::ZERO);
            let r = d_sqfreee;
            let r_str = {
                if r == Integer::ONE {
                    String::from("i")
                } else {
                    String::from("i√") + &r.to_string()
                }
            };

            let mut tight_c_abs = self.tight_c.clone();
            if tight_c_abs < Rational::ZERO {
                tight_c_abs = -tight_c_abs;
            }

            let mut tight_d_abs = self.tight_d.clone();
            if tight_d_abs < Rational::ZERO {
                tight_d_abs = -tight_d_abs;
            }

            let sign = tight_c_abs < tight_d_abs;
            if x == Rational::ZERO {
                if y == Rational::ONE {
                    write!(
                        f,
                        "{}{}",
                        match sign {
                            true => "",
                            false => "-",
                        },
                        r_str
                    )?;
                } else {
                    write!(
                        f,
                        "{}{}{}",
                        match sign {
                            true => "",
                            false => "-",
                        },
                        y,
                        r_str
                    )?;
                }
            } else {
                if y == Rational::ONE {
                    write!(
                        f,
                        "{}{}{}",
                        x,
                        match sign {
                            true => "+",
                            false => "-",
                        },
                        r_str
                    )?;
                } else {
                    write!(
                        f,
                        "{}{}{}{}",
                        x,
                        match sign {
                            true => "+",
                            false => "-",
                        },
                        y,
                        r_str
                    )?;
                }
            }
        } else {
            let mut root = self.clone();
            root.refine_to_accuracy(&Rational::from_integers(
                Integer::from(1),
                Integer::from(100),
            ));

            let m_re = (&root.tight_a + &root.tight_b) / Rational::TWO;
            let m_im = (&root.tight_c + &root.tight_d) / Rational::TWO;

            write!(f, "≈")?;
            write!(f, "{}", m_re.decimal_string_approx())?;
            // write!(f, "±");
            // write!(f, "{}", decimal_string_approx(self.accuracy_re() / Rational::TWO));
            if m_im >= Rational::ZERO {
                write!(f, "+")?;
            }
            write!(f, "{}", m_im.decimal_string_approx())?;
            // write!(f, "±");
            // write!(f, "{}", decimal_string_approx(self.accuracy_im() / Rational::TWO));
            write!(f, "i")?;
        }
        Ok(())
    }
}

impl ComplexAlgebraicRoot {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        if !(self.tight_a < self.tight_b) {
            return Err("tight a should be strictly less than b");
        }
        if !(self.tight_c < self.tight_d) {
            return Err("tight c should be strictly less than d");
        }
        // if !(self.wide_a < self.wide_b) {
        //     return Err("wide a should be strictly less than b");
        // }
        // if !(self.wide_c < self.wide_d) {
        //     return Err("wide c should be strictly less than d");
        // }

        if !self.poly.is_irreducible() {
            return Err("Isolated complex root minimal polynomial should be irreducible");
        }

        if self.poly.degree().unwrap() < 2 {
            return Err("Isolated complex root minimal polynomial should have degree at least 2");
        }

        match self.poly.count_complex_roots(
            &self.tight_a,
            &self.tight_b,
            &self.tight_c,
            &self.tight_d,
        ) {
            Some(1) => {}
            Some(_) => {
                return Err("Isolated complex root must exactly 1 root with none on the boundary");
            }
            None => {
                return Err(
                    "Isolated complex root must contain exactly 1 root with none on the boundary",
                );
            }
        }

        let real_roots_in_box = if self.tight_c < Rational::ZERO && Rational::ZERO < self.tight_d {
            self.poly
                .all_real_roots()
                .into_iter()
                .filter(|x| {
                    &RealAlgebraic::Rational(self.tight_a.clone()) < x
                        && x < &RealAlgebraic::Rational(self.tight_b.clone())
                })
                .collect::<Vec<_>>()
        } else {
            vec![]
        };
        if !real_roots_in_box.is_empty() {
            return Err("Isolated complex root must not be a real root");
        }

        Ok(())
    }

    pub fn accuracy_re(&self) -> Rational {
        &self.tight_b - &self.tight_a
    }

    pub fn accuracy_im(&self) -> Rational {
        &self.tight_d - &self.tight_c
    }

    fn conj(mut self) -> Self {
        (self.tight_c, self.tight_d) = (-self.tight_d, -self.tight_c);
        self
    }

    pub fn neg_mut(&mut self) {
        self.poly = Polynomial::compose(
            &self.poly,
            &Polynomial::from_coeffs(vec![Integer::from(0), Integer::from(-1)]),
        )
        .fav_assoc();
        (self.tight_a, self.tight_b) = (-self.tight_b.clone(), -self.tight_a.clone());
        (self.tight_c, self.tight_d) = (-self.tight_d.clone(), -self.tight_c.clone());
    }

    pub fn neg(mut self) -> Self {
        self.neg_mut();
        self
    }

    pub fn refine(&mut self) {
        let ((n1, a1, b1, c1, d1), (n2, a2, b2, c2, d2)) = bisect_box(
            &self.poly,
            1,
            &self.tight_a,
            &self.tight_b,
            &self.tight_c,
            &self.tight_d,
        );

        match (n1, n2) {
            (1, 0) => {
                self.tight_a = a1;
                self.tight_b = b1;
                self.tight_c = c1;
                self.tight_d = d1;
            }
            (0, 1) => {
                self.tight_a = a2;
                self.tight_b = b2;
                self.tight_c = c2;
                self.tight_d = d2;
            }
            _ => {
                unreachable!();
            }
        }

        #[cfg(debug_assertions)]
        self.check_invariants().unwrap();
    }

    pub fn refine_to_accuracy(&mut self, accuracy: &Rational) {
        while &self.accuracy_re() > accuracy || &self.accuracy_im() > accuracy {
            self.refine();
        }
    }

    pub fn equal(&self, other: &Self) -> bool {
        let poly = &self.poly;
        //polys should be irreducible primitive fav-assoc so this is valid
        if poly == &other.poly {
            //find the overlap of the two isolating boxes

            /*
                There are only two possibilities:
                 - A and C contain a root and B is empty, so self != other
                 - B contains a root and both A and C are empty, so self == other

                          +----------------+
                          |                |
                          |           C    |
                +---------+--------+       |
                |         |   B    |       |
                |         +--------+-------+
                |   A              |
                |                  |
                +------------------+

            */

            //is there any overlap at all?
            if self.tight_a < other.tight_b
                && other.tight_a < self.tight_b
                && self.tight_c < other.tight_d
                && other.tight_c < self.tight_d
            {
                let overlap_a = std::cmp::max(&self.tight_a, &other.tight_a);
                let overlap_b = std::cmp::min(&self.tight_b, &other.tight_b);
                let overlap_c = std::cmp::max(&self.tight_c, &other.tight_c);
                let overlap_d = std::cmp::min(&self.tight_d, &other.tight_d);

                let n = poly
                    .count_complex_roots(overlap_a, overlap_b, overlap_c, overlap_d)
                    .unwrap();

                return match n {
                    0 => false,
                    1 => true,
                    _ => panic!(),
                };
            }
        }
        false
    }

    pub fn apply_poly(&mut self, poly: &Polynomial<Rational>) -> ComplexAlgebraic {
        let poly = Polynomial::rem(poly, &self.min_poly());
        match poly.as_constant() {
            Some(rat) => ComplexAlgebraic::Real(RealAlgebraic::Rational(rat)),
            None => {
                let ans_poly = self
                    .min_poly()
                    .algebraic_number_field()
                    .min_poly(&poly)
                    .primitive_part_fof();

                identify_complex_root(
                    ans_poly,
                    (0..).map(|i| {
                        if i != 0 {
                            self.refine();
                        }

                        // eg: c + bx + ax^2 = c + x(b + x(a))
                        let mut coeffs = poly.coeffs().into_iter().rev();
                        let lc = coeffs.next().unwrap();
                        let mut ans = mul_box_rat(
                            (&self.tight_a, &self.tight_b, &self.tight_c, &self.tight_d),
                            lc,
                        );
                        for (i, c) in coeffs.enumerate() {
                            if i != 0 {
                                ans = mul_boxes(
                                    (&ans.0, &ans.1, &ans.2, &ans.3),
                                    (&self.tight_a, &self.tight_b, &self.tight_c, &self.tight_d),
                                );
                            }
                            ans = add_box_rat((&ans.0, &ans.1, &ans.2, &ans.3), c);
                        }

                        ans
                    }),
                )
            }
        }
    }

    pub fn min_poly(&self) -> Polynomial<Rational> {
        self.poly.apply_map(|c| Rational::from(c)).fav_assoc()
    }
}

#[derive(Debug, Clone, CannonicalStructure)]
pub enum ComplexAlgebraic {
    Real(RealAlgebraic),
    Complex(ComplexAlgebraicRoot),
}

impl ComplexAlgebraic {
    pub fn i() -> Self {
        let i = ComplexAlgebraic::Complex(ComplexAlgebraicRoot {
            tight_a: Rational::from_integers(Integer::from(-1), Integer::from(1)),
            tight_b: Rational::from_integers(Integer::from(1), Integer::from(1)),
            tight_c: Rational::from_integers(Integer::from(0), Integer::from(1)),
            tight_d: Rational::from_integers(Integer::from(2), Integer::from(1)),
            poly: Polynomial::from_coeffs(vec![Integer::ONE, Integer::ZERO, Integer::ONE]),
        });
        #[cfg(debug_assertions)]
        i.check_invariants().unwrap();
        i
    }
}

impl ComplexAlgebraic {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        match self {
            ComplexAlgebraic::Real(x) => match x.check_invariants() {
                Ok(()) => {}
                Err(e) => {
                    return Err(e);
                }
            },
            ComplexAlgebraic::Complex(x) => match x.check_invariants() {
                Ok(()) => {}
                Err(e) => {
                    return Err(e);
                }
            },
        }
        Ok(())
    }

    pub fn eq_mut(&mut self, other: &mut Self) -> bool {
        match (self, other) {
            (ComplexAlgebraic::Real(a), ComplexAlgebraic::Real(b)) => {
                RealAlgebraic::cmp_mut(a, b).is_eq()
            }
            (ComplexAlgebraic::Real(_a), ComplexAlgebraic::Complex(_b)) => false,
            (ComplexAlgebraic::Complex(_a), ComplexAlgebraic::Real(_b)) => false,
            (ComplexAlgebraic::Complex(a), ComplexAlgebraic::Complex(b)) => a.equal(b),
        }
    }

    pub fn apply_poly(&mut self, poly: &Polynomial<Rational>) -> Self {
        match self {
            ComplexAlgebraic::Real(reat_alg) => ComplexAlgebraic::Real(reat_alg.apply_poly(poly)),
            ComplexAlgebraic::Complex(cpx_root) => cpx_root.apply_poly(poly),
        }
    }

    pub fn degree(&self) -> usize {
        self.min_poly().degree().unwrap()
    }
}

impl PartialEq for ComplexAlgebraic {
    fn eq(&self, other: &Self) -> bool {
        Self::eq_mut(&mut self.clone(), &mut other.clone())
    }
}

impl Eq for ComplexAlgebraic {}

impl Display for ComplexAlgebraic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ComplexAlgebraic::Real(a) => write!(f, "{}", a),
            ComplexAlgebraic::Complex(a) => write!(f, "{}", a),
        }
    }
}

impl PartialEqStructure for ComplexAlgebraicCannonicalStructure {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        a == b
    }
}

impl EqStructure for ComplexAlgebraicCannonicalStructure {}

impl ToStringStructure for ComplexAlgebraicCannonicalStructure {
    fn to_string(&self, elem: &Self::Set) -> String {
        format!("{}", elem)
    }
}

impl SemiRingStructure for ComplexAlgebraicCannonicalStructure {
    fn zero(&self) -> Self::Set {
        ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::zero()))
    }

    fn one(&self) -> Self::Set {
        ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::one()))
    }

    fn add(&self, alg1: &Self::Set, alg2: &Self::Set) -> Self::Set {
        // println!("add {:?} {:?}", alg1, alg2);
        // alg1.check_invariants().unwrap();
        // alg2.check_invariants().unwrap();

        fn add_real(mut cpx: ComplexAlgebraicRoot, real: RealAlgebraic) -> ComplexAlgebraic {
            match real {
                RealAlgebraic::Rational(rat) => {
                    cpx.tight_a += &rat;
                    cpx.tight_b += &rat;
                    //compose with x - n/d = dx - n
                    cpx.poly = Polynomial::compose(
                        &cpx.poly.apply_map(|c| Rational::from(c)),
                        &Polynomial::from_coeffs(vec![-rat, Rational::ONE]),
                    )
                    .primitive_part_fof();

                    let cpx = ComplexAlgebraic::Complex(cpx);

                    #[cfg(debug_assertions)]
                    cpx.check_invariants().unwrap();
                    cpx
                }
                RealAlgebraic::Real(mut real) => identify_complex_root(
                    root_sum_poly(&cpx.poly, real.poly()),
                    (0..).map(|i| {
                        if i != 0 {
                            cpx.refine();
                            real.refine();
                        }
                        let ans_tight_a = &cpx.tight_a + real.tight_a();
                        let ans_tight_b = &cpx.tight_b + real.tight_b();
                        let ans_tight_c = cpx.tight_c.clone();
                        let ans_tight_d = cpx.tight_d.clone();
                        (ans_tight_a, ans_tight_b, ans_tight_c, ans_tight_d)
                    }),
                ),
            }
        }

        match (alg1, alg2) {
            (ComplexAlgebraic::Real(real1), ComplexAlgebraic::Real(real2)) => {
                ComplexAlgebraic::Real(RealAlgebraic::add(real1, real2))
            }
            (ComplexAlgebraic::Real(real1), ComplexAlgebraic::Complex(cpx2)) => {
                add_real(cpx2.clone(), real1.clone())
            }
            (ComplexAlgebraic::Complex(cpx1), ComplexAlgebraic::Real(real2)) => {
                add_real(cpx1.clone(), real2.clone())
            }
            (ComplexAlgebraic::Complex(cpx1), ComplexAlgebraic::Complex(cpx2)) => {
                let mut cpx1 = cpx1.clone();
                let mut cpx2 = cpx2.clone();

                identify_complex_root(
                    root_sum_poly(&cpx1.poly, &cpx2.poly),
                    (0..).map(|i| {
                        if i != 0 {
                            cpx1.refine();
                            cpx2.refine();
                        }
                        add_boxes(
                            (&cpx1.tight_a, &cpx1.tight_b, &cpx1.tight_c, &cpx1.tight_d),
                            (&cpx2.tight_a, &cpx2.tight_b, &cpx2.tight_c, &cpx2.tight_d),
                        )
                    }),
                )
            }
        }
    }

    fn mul(&self, alg1: &Self::Set, alg2: &Self::Set) -> Self::Set {
        // println!("mul {:?} {:?}", alg1, alg2);
        // alg1.check_invariants().unwrap();
        // alg2.check_invariants().unwrap();

        fn mul_real(mut cpx: ComplexAlgebraicRoot, real: RealAlgebraic) -> ComplexAlgebraic {
            match real {
                RealAlgebraic::Rational(rat) => match rat.cmp(&Rational::ZERO) {
                    std::cmp::Ordering::Equal => ComplexAlgebraic::zero(),
                    std::cmp::Ordering::Less => mul_real(cpx.neg(), RealAlgebraic::Rational(-rat)),
                    std::cmp::Ordering::Greater => {
                        debug_assert!(rat > Rational::ZERO);
                        cpx.tight_a *= &rat;
                        cpx.tight_b *= &rat;
                        cpx.tight_c *= &rat;
                        cpx.tight_d *= &rat;
                        cpx.poly = root_rat_mul_poly(cpx.poly, &rat);
                        #[cfg(debug_assertions)]
                        assert!(cpx.check_invariants().is_ok());
                        ComplexAlgebraic::Complex(cpx)
                    }
                },
                RealAlgebraic::Real(mut real) => {
                    identify_complex_root(
                        root_product_poly(&cpx.poly, real.poly()),
                        (0..).map(|i| {
                            if i != 0 {
                                cpx.refine();
                                real.refine();
                            }

                            let mut pts_re = vec![];
                            let mut pts_im = vec![];
                            for (re, im) in [
                                (&cpx.tight_a, &cpx.tight_c),
                                (&cpx.tight_a, &cpx.tight_d),
                                (&cpx.tight_b, &cpx.tight_c),
                                (&cpx.tight_b, &cpx.tight_d),
                            ] {
                                for t in [real.tight_a(), real.tight_b()] {
                                    pts_re.push(t * re);
                                    pts_im.push(t * im);
                                }
                            }

                            let ans_tight_a = pts_re.iter().min().unwrap().clone();
                            let ans_tight_b = pts_re.into_iter().max().unwrap();
                            let ans_tight_c = pts_im.iter().min().unwrap().clone();
                            let ans_tight_d = pts_im.into_iter().max().unwrap();

                            //allow the bounds to expand a little bit so that the bounds are simpler
                            let diff_re = &ans_tight_b - &ans_tight_a;
                            let diff_im = &ans_tight_d - &ans_tight_c;
                            let ans_tight_a = Rational::simplest_rational_in_closed_interval(
                                &(&ans_tight_a - (&diff_re * Rational::from_str("1/10").unwrap())),
                                &ans_tight_a,
                            );
                            let ans_tight_b = Rational::simplest_rational_in_closed_interval(
                                &ans_tight_b,
                                &(&ans_tight_b + (&diff_re * Rational::from_str("1/10").unwrap())),
                            );
                            let ans_tight_c = Rational::simplest_rational_in_closed_interval(
                                &(&ans_tight_c - (&diff_im * Rational::from_str("1/10").unwrap())),
                                &ans_tight_c,
                            );
                            let ans_tight_d = Rational::simplest_rational_in_closed_interval(
                                &ans_tight_d,
                                &(&ans_tight_d + (&diff_im * Rational::from_str("1/10").unwrap())),
                            );
                            (ans_tight_a, ans_tight_b, ans_tight_c, ans_tight_d)
                        }),
                    )
                }
            }
        }

        match (alg1, alg2) {
            (ComplexAlgebraic::Real(real1), ComplexAlgebraic::Real(real2)) => {
                ComplexAlgebraic::Real(RealAlgebraic::mul(real1, real2))
            }
            (ComplexAlgebraic::Real(real1), ComplexAlgebraic::Complex(cpx2)) => {
                mul_real(cpx2.clone(), real1.clone())
            }
            (ComplexAlgebraic::Complex(cpx1), ComplexAlgebraic::Real(real2)) => {
                mul_real(cpx1.clone(), real2.clone())
            }
            (ComplexAlgebraic::Complex(cpx1), ComplexAlgebraic::Complex(cpx2)) => {
                let mut cpx1 = cpx1.clone();
                let mut cpx2 = cpx2.clone();

                identify_complex_root(
                    root_product_poly(&cpx1.poly, &cpx2.poly),
                    (0..).map(|i| {
                        if i != 0 {
                            cpx1.refine();
                            cpx2.refine();
                        }

                        let (ans_tight_a, ans_tight_b, ans_tight_c, ans_tight_d) = mul_boxes(
                            (&cpx1.tight_a, &cpx1.tight_b, &cpx1.tight_c, &cpx1.tight_d),
                            (&cpx2.tight_a, &cpx2.tight_b, &cpx2.tight_c, &cpx2.tight_d),
                        );

                        //allow the bounds to expand a little bit so that the bounds are simpler
                        let diff_re = &ans_tight_b - &ans_tight_a;
                        let diff_im = &ans_tight_d - &ans_tight_c;
                        let ans_tight_a = Rational::simplest_rational_in_closed_interval(
                            &(&ans_tight_a - (&diff_re * Rational::from_str("1/10").unwrap())),
                            &ans_tight_a,
                        );
                        let ans_tight_b = Rational::simplest_rational_in_closed_interval(
                            &ans_tight_b,
                            &(&ans_tight_b + (&diff_re * Rational::from_str("1/10").unwrap())),
                        );
                        let ans_tight_c = Rational::simplest_rational_in_closed_interval(
                            &(&ans_tight_c - (&diff_im * Rational::from_str("1/10").unwrap())),
                            &ans_tight_c,
                        );
                        let ans_tight_d = Rational::simplest_rational_in_closed_interval(
                            &ans_tight_d,
                            &(&ans_tight_d + (&diff_im * Rational::from_str("1/10").unwrap())),
                        );

                        (ans_tight_a, ans_tight_b, ans_tight_c, ans_tight_d)
                    }),
                )
            }
        }
    }
}

impl RingStructure for ComplexAlgebraicCannonicalStructure {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        match a {
            ComplexAlgebraic::Real(root) => ComplexAlgebraic::Real(RealAlgebraic::neg(root)),
            ComplexAlgebraic::Complex(root) => ComplexAlgebraic::Complex(root.clone().neg()),
        }
    }
}

impl UnitsStructure for ComplexAlgebraicCannonicalStructure {
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        // println!("inv {:?}", a);
        // a.check_invariants().unwrap();

        match a {
            ComplexAlgebraic::Real(a) => Ok(ComplexAlgebraic::Real(a.inv()?)),
            ComplexAlgebraic::Complex(a) => {
                let mut root = a.clone();

                let inv_poly = Polynomial::from_coeffs(
                    root.poly.clone().into_coeffs().into_iter().rev().collect(),
                )
                .fav_assoc();
                debug_assert!(inv_poly.is_irreducible());

                // println!("root = {} = {:?}", root, root);
                // println!("inv_poly = {}", inv_poly);
                // for root in inv_poly.all_complex_roots() {
                //     println!("root = {}", root);
                // }

                // z = the root represented by self
                // a = the center of the approximation
                // eps = the radius of the approximation

                //result used to get region which contains the reciprocal
                // if for some 0 < lambda < 1 we have |z - a| < eps <= (1 - lambda) |a| then |1/z - 1/a| < eps / lambda |a|^2

                loop {
                    //estimate eps such that |z - a| < eps
                    let eps = ((&root.tight_b - &root.tight_a) + (&root.tight_d - &root.tight_c))
                        / Rational::TWO;

                    let w_re = (&root.tight_a + &root.tight_b) / Rational::TWO;
                    let w_im = (&root.tight_c + &root.tight_d) / Rational::TWO;
                    let w_mag_sq = &w_re * &w_re + &w_im * &w_im;

                    // println!(
                    //     "w = {} + {} i",
                    //     decimal_string_approx(w_re.clone()),
                    //     decimal_string_approx(w_im.clone())
                    // );

                    //refine until eps < |a|
                    if &eps * &eps > w_mag_sq {
                        root.refine();
                        continue;
                    }

                    // find 0 < lambda < 1 such that eps < (1 - lambda) * |a|
                    let mut lambda = Rational::ONE_HALF;
                    while &eps * &eps >= (Rational::ONE - &lambda) * &w_mag_sq {
                        lambda = lambda * &Rational::ONE_HALF;
                    }
                    let lambda = lambda;
                    // println!("lambda = {}", lambda);

                    // |1/z - 1/a| < eps / lambda |a|^2

                    let w_recip_re = w_re / &w_mag_sq;
                    let w_recip_im = -w_im / &w_mag_sq;
                    let delta = &eps / (lambda * w_mag_sq);

                    let inv_tight_a = &w_recip_re - &delta;
                    let inv_tight_b = &w_recip_re + &delta;
                    let inv_tight_c = &w_recip_im - &delta;
                    let inv_tight_d = &w_recip_im + &delta;

                    // println!(
                    //     "w_inv = {} + {} i",
                    //     decimal_string_approx(w_recip_re),
                    //     decimal_string_approx(w_recip_im)
                    // );
                    // println!(
                    //     "eps = {}  delta = {}",
                    //     decimal_string_approx(eps),
                    //     decimal_string_approx(delta)
                    // );

                    match inv_poly.count_complex_roots(
                        &inv_tight_a,
                        &inv_tight_b,
                        &inv_tight_c,
                        &inv_tight_d,
                    ) {
                        Some(count) => {
                            if count == 0 {
                                unreachable!();
                            } else if count == 1 {
                                let ans = ComplexAlgebraic::Complex(ComplexAlgebraicRoot {
                                    tight_a: inv_tight_a,
                                    tight_b: inv_tight_b,
                                    tight_c: inv_tight_c,
                                    tight_d: inv_tight_d,
                                    poly: inv_poly,
                                });
                                #[cfg(debug_assertions)]
                                ans.check_invariants().unwrap();
                                return Ok(ans);
                            }
                        }
                        None => {}
                    }

                    root.refine();
                }
            }
        }
    }
}

impl IntegralDomainStructure for ComplexAlgebraicCannonicalStructure {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        Ok(self.mul(a, &self.inv(b)?))
    }
}

impl FieldStructure for ComplexAlgebraicCannonicalStructure {}

impl CharZeroStructure for ComplexAlgebraicCannonicalStructure {
    fn try_to_int(&self, alg: &Self::Set) -> Option<Integer> {
        match alg {
            ComplexAlgebraic::Real(real_alg) => real_alg.try_to_int(),
            ComplexAlgebraic::Complex(_) => None,
        }
    }
}

impl ComplexSubsetStructure for ComplexAlgebraicCannonicalStructure {}

impl ComplexConjugateStructure for ComplexAlgebraicCannonicalStructure {
    fn conjugate(&self, x: &Self::Set) -> Self::Set {
        match x {
            ComplexAlgebraic::Real(x) => ComplexAlgebraic::Real(x.clone()),
            ComplexAlgebraic::Complex(x) => ComplexAlgebraic::Complex(x.clone().conj()),
        }
    }
}

impl PositiveRealNthRootStructure for ComplexAlgebraicCannonicalStructure {
    fn nth_root(&self, x: &Self::Set, n: usize) -> Result<Self::Set, ()> {
        match x {
            ComplexAlgebraic::Real(x) => Ok(ComplexAlgebraic::Real(x.nth_root(n)?)),
            ComplexAlgebraic::Complex(_) => Err(()),
        }
    }
}

impl ComplexAlgebraic {
    pub fn min_poly(&self) -> Polynomial<Rational> {
        match self {
            ComplexAlgebraic::Real(real_algebraic) => real_algebraic.min_poly(),
            ComplexAlgebraic::Complex(complex_root) => complex_root.min_poly(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_apply_poly() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();

        let f = (2 * x.pow(2) - 4 * x - 3).into_verbose();
        let g = (3 * x.pow(3) + 7 * x - 1).into_verbose();
        //h(x) = f(g(x))
        let h = Polynomial::compose(&f, &g);

        println!("f = {}", f);
        println!("g = {}", g);
        println!("h = {}", h);

        for x in h.primitive_part_fof().all_complex_roots() {
            println!("");
            println!("x = {} deg = {}", x, x.min_poly().degree().unwrap());
            let gx = x.clone().apply_poly(&g);
            println!("gx = {}", gx);
            let fgx = gx.clone().apply_poly(&f);
            println!("fgx = {}", fgx);
            debug_assert_eq!(fgx, ComplexAlgebraic::zero());
        }

        let i = &ComplexAlgebraic::i().into_ergonomic();
        let a = (2 + 3 * i).into_verbose();
        let f = (x.pow(10)).into_verbose();
        assert_eq!(
            a.clone().apply_poly(&f),
            (-341525 - 145668 * i).into_verbose()
        );
    }

    #[test]
    fn test_all_complex_roots() {
        let f = Polynomial::<Integer>::from_coeffs(vec![-1, -1, 0, 0, 0, 1]);
        let roots = f.all_complex_roots();
        assert_eq!(roots.len(), Polynomial::degree(&f).unwrap());
        for root in &roots {
            root.check_invariants().unwrap();
        }
    }

    #[test]
    fn test_complex_equal() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();
        let f = (x.pow(5) - x + 1).into_verbose();
        assert!(f.is_irreducible());
        let mut count = 0;
        for a in f.all_complex_roots() {
            for b in f.all_complex_roots() {
                if a == b {
                    count += 1;
                }
            }
        }
        assert_eq!(count, f.degree().unwrap());
    }

    #[test]
    fn test_complex_root_sum() {
        let f = Polynomial::<Integer>::from_coeffs(vec![1, 0, 0, 1]);
        let roots = f.all_complex_roots();
        let s = ComplexAlgebraic::sum(roots.iter().collect());
        println!("{:?}", s);
        assert_eq!(s, ComplexAlgebraic::zero());

        let f = Polynomial::<Integer>::from_coeffs(vec![7, -3, 42, 9]);
        println!("f = {}", f);
        let roots = f.all_complex_roots();
        for root in &roots {
            println!("root = {}", root);
        }
        let s = ComplexAlgebraic::sum(roots.iter().collect());
        println!("root sum = {:?}", s);
        assert_eq!(
            s,
            ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from_integers(
                Integer::from(-14),
                Integer::from(3)
            )))
        );
    }

    #[test]
    fn test_complex_mul() {
        let f = Polynomial::<Integer>::from_coeffs(vec![1, 0, 0, 1]);
        let roots = f.all_complex_roots();
        let s = ComplexAlgebraic::product(roots.iter().collect());
        println!("{:?}", s);
        assert_eq!(s, ComplexAlgebraic::one().neg());

        let f = Polynomial::<Integer>::from_coeffs(vec![7, -3, 42, 9]);
        println!("f = {}", f);
        let roots = f.all_complex_roots();
        for root in &roots {
            println!("root = {}", root);
        }
        let s = ComplexAlgebraic::product(roots.iter().collect());
        println!("root prod = {:?}", s);
        assert_eq!(
            s,
            ComplexAlgebraic::Real(RealAlgebraic::Rational(Rational::from_integers(
                Integer::from(-7),
                Integer::from(9)
            )))
        );
    }

    #[test]
    fn test_complex_inv() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();
        let f = (x.pow(4) - x + 1).into_verbose();

        for root in f.all_complex_roots() {
            assert_eq!(
                ComplexAlgebraic::mul(&root.inv().unwrap(), &root),
                ComplexAlgebraic::one()
            );
        }
    }
}
