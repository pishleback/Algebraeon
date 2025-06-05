use super::*;
use algebraeon_nzq::traits::{Abs, Fraction};

fn unique_linear_root(poly: &Polynomial<Integer>) -> Rational {
    debug_assert_eq!(poly.degree().unwrap(), 1);
    -Rational::from_integers(poly.coeff(0), poly.coeff(1))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Interleave {
    First,
    Second,
}

#[derive(Debug, Clone)]
pub enum SquarefreePolyRealRootInterval {
    Rational(Rational),
    //lower bound, upper bound, increasing
    //increasing = false : decreasing i.e. poly(a) > poly(b), true : increasing i.e. poly(a) < poly(b)
    Real(Rational, Rational, bool),
}

#[derive(Debug, Clone)]
pub struct SquarefreePolyRealRoots {
    poly_sqfr: Polynomial<Integer>,
    //an ordered list of isolating intervals for the squarefree polynomial
    //e.g. if r represents a real algebraic number and | represents a rational root
    //        (      r    )      |  ( r     )   |   |   (        r   )
    //note: it is allowed that some r might actually be rational but not known to be
    intervals: Vec<SquarefreePolyRealRootInterval>,
}

impl SquarefreePolyRealRoots {
    pub fn interval(&self, idx: usize) -> &SquarefreePolyRealRootInterval {
        &self.intervals[idx]
    }

    pub fn intervals_len(&self) -> usize {
        self.intervals.len()
    }
}

impl SquarefreePolyRealRoots {
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        //poly should be squarefree
        if self
            .poly_sqfr
            .clone()
            .primitive_squarefree_part()
            .degree()
            .unwrap()
            != self.poly_sqfr.degree().unwrap()
        {
            return Err("poly should be squarefree");
        }

        //check the isolating intervals
        if !self.intervals.is_empty() {
            for i in 0..self.intervals.len() - 1 {
                let int1 = &self.intervals[i];
                let int2 = &self.intervals[i + 1];
                match (int1, int2) {
                    (
                        SquarefreePolyRealRootInterval::Rational(a),
                        SquarefreePolyRealRootInterval::Rational(x),
                    ) => {
                        if a >= x {
                            return Err("interval values should be strictly increasing");
                        }
                    }
                    (
                        SquarefreePolyRealRootInterval::Rational(a),
                        SquarefreePolyRealRootInterval::Real(x, y, _),
                    ) => {
                        if a >= x {
                            return Err("interval values should be strictly increasing");
                        }
                        if x >= y {
                            return Err("interval values should be strictly increasing");
                        }
                    }
                    (
                        SquarefreePolyRealRootInterval::Real(a, b, _),
                        SquarefreePolyRealRootInterval::Rational(x),
                    ) => {
                        if a >= b {
                            return Err("interval values should be strictly increasing");
                        }
                        if b >= x {
                            return Err("interval values should be strictly increasing");
                        }
                    }
                    (
                        SquarefreePolyRealRootInterval::Real(a, b, _),
                        SquarefreePolyRealRootInterval::Real(x, y, _),
                    ) => {
                        if a >= b {
                            return Err("interval values should be strictly increasing");
                        }
                        if b > x {
                            return Err("interval values should be increasing");
                        }
                        if x >= y {
                            return Err("interval values should be strictly increasing");
                        }
                    }
                }
            }
        }

        for interval in &self.intervals {
            match interval {
                SquarefreePolyRealRootInterval::Rational(a) => {
                    if evaluate_at_rational(&self.poly_sqfr, a) != Rational::from(0) {
                        return Err("poly should be zero at a rational root");
                    }
                }
                SquarefreePolyRealRootInterval::Real(a, b, incr) => {
                    let at_a = evaluate_at_rational(&self.poly_sqfr, a);
                    let at_b = evaluate_at_rational(&self.poly_sqfr, b);

                    if at_a == Rational::from(0) || at_b == Rational::from(0) {
                        return Err("poly should not be zero at boundary of isolating interval");
                    }

                    if (at_a > Rational::from(0)) == (at_b > Rational::from(0)) {
                        return Err("sign of poly should be different at a and at b");
                    }

                    #[allow(clippy::collapsible_else_if)]
                    if *incr {
                        if !((at_a < Rational::from(0)) && (at_b > Rational::from(0))) {
                            return Err("sign of poly should go from neg to pos here");
                        }
                    } else {
                        if !((at_a > Rational::from(0)) && (at_b < Rational::from(0))) {
                            return Err("sign of poly should go from pos to neg here");
                        }
                    }
                }
            }
        }

        Ok(())
    }

    fn get_wide_interval(&self, idx: usize) -> (LowerBound, UpperBound) {
        assert!(idx < self.intervals.len());

        let wide_a = {
            if idx == 0 {
                LowerBound::Inf
            } else {
                LowerBound::Finite({
                    match &self.intervals[idx - 1] {
                        SquarefreePolyRealRootInterval::Rational(a) => a.clone(),
                        SquarefreePolyRealRootInterval::Real(_, prev_b, _) => prev_b.clone(),
                    }
                })
            }
        };
        let wide_b = {
            if idx == self.intervals.len() - 1 {
                UpperBound::Inf
            } else {
                UpperBound::Finite({
                    match &self.intervals[idx + 1] {
                        SquarefreePolyRealRootInterval::Rational(a) => a.clone(),
                        SquarefreePolyRealRootInterval::Real(prev_a, _, _) => prev_a.clone(),
                    }
                })
            }
        };
        debug_assert!(wide_a.clone() < wide_b.clone());
        (wide_a, wide_b)
    }

    pub fn refine(&mut self, idx: usize) {
        assert!(idx < self.intervals.len());

        match &mut self.intervals[idx] {
            SquarefreePolyRealRootInterval::Rational(_a) => {}
            SquarefreePolyRealRootInterval::Real(a, b, dir) => {
                let m = RationalSimpleBetweenGenerator::new_between_middle(
                    a,
                    b,
                    &Rational::from_integers(Integer::from(1), Integer::from(3)),
                )
                .next()
                .unwrap();

                match evaluate_at_rational(&self.poly_sqfr, &m).cmp(&Rational::from(0)) {
                    std::cmp::Ordering::Less => match dir {
                        true => {
                            *a = m;
                        }
                        false => {
                            *b = m;
                        }
                    },
                    std::cmp::Ordering::Equal => {
                        self.intervals[idx] = SquarefreePolyRealRootInterval::Rational(m);
                    }
                    std::cmp::Ordering::Greater => match dir {
                        true => {
                            *b = m;
                        }
                        false => {
                            *a = m;
                        }
                    },
                }
            }
        }
    }

    pub fn refine_all(&mut self) {
        for idx in 0..self.intervals.len() {
            self.refine(idx);
        }
    }

    fn get_real_root(&self, idx: usize) -> RealAlgebraic {
        // println!("get_real_root of {:?} at {}", self, idx);
        // println!("{}", self.poly_sqfr);
        assert!(idx < self.intervals.len());
        match &self.intervals[idx] {
            SquarefreePolyRealRootInterval::Rational(rat) => RealAlgebraic::Rational(rat.clone()),
            SquarefreePolyRealRootInterval::Real(a, b, _dir) => {
                let (_unit, factors) = self.poly_sqfr.factor().unwrap().into_unit_and_powers();
                for (factor, k) in factors {
                    // println!("factor = {}", factor);
                    debug_assert_eq!(k, Natural::ONE); //square free
                    let deg = factor.degree().unwrap();
                    debug_assert_ne!(deg, 0);
                    let at_a = evaluate_at_rational(&factor, a);
                    let at_b = evaluate_at_rational(&factor, b);
                    debug_assert_ne!(at_a, Rational::ZERO);
                    debug_assert_ne!(at_b, Rational::ZERO);
                    let sign_a = at_a >= Rational::ZERO;
                    let sign_b = at_b >= Rational::ZERO;
                    // println!("at_a = {}", at_a);
                    // println!("at_b = {}", at_b);
                    #[allow(clippy::collapsible_else_if)]
                    if deg == 1 {
                        if sign_a != sign_b {
                            return RealAlgebraic::Rational(unique_linear_root(&factor));
                        }
                    } else {
                        if sign_a && !sign_b {
                            let (wide_a, wide_b) = self.get_wide_interval(idx);
                            return RealAlgebraic::Real(RealAlgebraicRoot {
                                poly: factor,
                                tight_a: a.clone(),
                                tight_b: b.clone(),
                                wide_a,
                                wide_b,
                                dir: false,
                            });
                        } else if !sign_a && sign_b {
                            let (wide_a, wide_b) = self.get_wide_interval(idx);
                            return RealAlgebraic::Real(RealAlgebraicRoot {
                                poly: factor,
                                tight_a: a.clone(),
                                tight_b: b.clone(),
                                wide_a,
                                wide_b,
                                dir: true,
                            });
                        }
                    }
                    debug_assert_eq!(sign_a, sign_b);
                }
                unreachable!()
            }
        }
    }

    pub fn to_real_roots(self) -> Vec<RealAlgebraic> {
        debug_assert!(self.poly_sqfr.is_irreducible());
        debug_assert!(self.poly_sqfr.is_fav_assoc());
        let deg = self.poly_sqfr.degree().unwrap();
        if deg == 0 {
            vec![]
        } else if deg == 1 {
            if self.intervals.is_empty() {
                vec![]
            } else if self.intervals.len() == 1 {
                match self.intervals.into_iter().next().unwrap() {
                    SquarefreePolyRealRootInterval::Rational(a) => {
                        vec![RealAlgebraic::Rational(a)]
                    }
                    SquarefreePolyRealRootInterval::Real(_, _, _) => {
                        // panic!("degree 1 polynomial should have had rational root found exactly");
                        vec![RealAlgebraic::Rational(unique_linear_root(&self.poly_sqfr))]
                    }
                }
            } else {
                panic!();
            }
        } else {
            let mut roots = vec![];
            for (idx, interval) in self.intervals.iter().enumerate() {
                roots.push({
                    match interval {
                        SquarefreePolyRealRootInterval::Rational(a) => {
                            RealAlgebraic::Rational(a.clone())
                        }
                        SquarefreePolyRealRootInterval::Real(tight_a, tight_b, dir) => {
                            let (wide_a, wide_b) = self.get_wide_interval(idx);
                            RealAlgebraic::Real(RealAlgebraicRoot {
                                poly: self.poly_sqfr.clone(),
                                tight_a: tight_a.clone(),
                                tight_b: tight_b.clone(),
                                wide_a,
                                wide_b,
                                dir: *dir,
                            })
                        }
                    }
                });
            }
            roots
        }
    }

    //separate the isolating intervals of the roots in roots1 and roots2
    //return Err if a root in roots1 and a root in roots2 are equal and thus cant be separated
    //ends of real root intervals should not equal rational roots in the other one
    pub fn separate(roots1: &mut Self, roots2: &mut Self) -> Result<Vec<(Interleave, usize)>, ()> {
        let poly_gcd_sqfr = Polynomial::subresultant_gcd(&roots1.poly_sqfr, &roots2.poly_sqfr);
        let (_, poly_gcd_sqfr) = poly_gcd_sqfr.factor_primitive().unwrap();

        //compute which of the roots are equals to some root from the other one
        let is_gcdroot_1: Vec<_> = roots1
            .intervals
            .iter()
            .map(|root| match root {
                SquarefreePolyRealRootInterval::Rational(x) => {
                    evaluate_at_rational(&poly_gcd_sqfr, x) == Rational::from(0)
                }
                SquarefreePolyRealRootInterval::Real(a, b, _dir) => {
                    debug_assert_ne!(evaluate_at_rational(&poly_gcd_sqfr, a), Rational::from(0));
                    debug_assert_ne!(evaluate_at_rational(&poly_gcd_sqfr, b), Rational::from(0));
                    (evaluate_at_rational(&poly_gcd_sqfr, a) > Rational::from(0))
                        != (evaluate_at_rational(&poly_gcd_sqfr, b) > Rational::from(0))
                }
            })
            .collect();
        let is_gcdroot_2: Vec<_> = roots2
            .intervals
            .iter()
            .map(|root| match root {
                SquarefreePolyRealRootInterval::Rational(x) => {
                    evaluate_at_rational(&poly_gcd_sqfr, x) == Rational::from(0)
                }
                SquarefreePolyRealRootInterval::Real(a, b, _dir) => {
                    debug_assert_ne!(evaluate_at_rational(&poly_gcd_sqfr, a), Rational::from(0));
                    debug_assert_ne!(evaluate_at_rational(&poly_gcd_sqfr, b), Rational::from(0));
                    (evaluate_at_rational(&poly_gcd_sqfr, a) > Rational::from(0))
                        != (evaluate_at_rational(&poly_gcd_sqfr, b) > Rational::from(0))
                }
            })
            .collect();

        //do the separation
        let mut all_roots = vec![];

        let mut idx1 = 0;
        let mut idx2 = 0;
        while idx1 < roots1.intervals.len() && idx2 < roots2.intervals.len() {
            let (wide1_a, wide1_b) = roots1.get_wide_interval(idx1);
            let (wide2_a, wide2_b) = roots2.get_wide_interval(idx2);

            let root1 = &roots1.intervals[idx1];
            let root2 = &roots2.intervals[idx2];

            //check if the roots are equal
            if is_gcdroot_1[idx1] && is_gcdroot_2[idx2] {
                match root1 {
                    SquarefreePolyRealRootInterval::Rational(x1) => {
                        if &wide2_a < x1 && x1 < &wide2_b {
                            return Err(());
                        }
                    }
                    SquarefreePolyRealRootInterval::Real(a1, b1, _dir1) => {
                        if &wide2_a < a1 && b1 < &wide2_b {
                            return Err(());
                        }
                    }
                }
                match root2 {
                    SquarefreePolyRealRootInterval::Rational(x2) => {
                        if &wide1_a < x2 && x2 < &wide1_b {
                            return Err(());
                        }
                    }
                    SquarefreePolyRealRootInterval::Real(a2, b2, _dir2) => {
                        if &wide1_a < a2 && b2 < &wide1_b {
                            return Err(());
                        }
                    }
                }
            }

            //check if one is bigger than the other
            match (root1, root2) {
                (
                    SquarefreePolyRealRootInterval::Rational(x1),
                    SquarefreePolyRealRootInterval::Rational(x2),
                ) => match x1.cmp(x2) {
                    std::cmp::Ordering::Less => {
                        all_roots.push((Interleave::First, idx1));
                        idx1 += 1;
                        continue;
                    }
                    std::cmp::Ordering::Equal => panic!(),
                    std::cmp::Ordering::Greater => {
                        all_roots.push((Interleave::Second, idx2));
                        idx2 += 1;
                        continue;
                    }
                },
                (
                    SquarefreePolyRealRootInterval::Rational(x1),
                    SquarefreePolyRealRootInterval::Real(a2, b2, _dir2),
                ) => {
                    if x1 < a2 {
                        all_roots.push((Interleave::First, idx1));
                        idx1 += 1;
                        continue;
                    }
                    if b2 < x1 {
                        all_roots.push((Interleave::Second, idx2));
                        idx2 += 1;
                        continue;
                    }
                }
                (
                    SquarefreePolyRealRootInterval::Real(a1, b1, _dir1),
                    SquarefreePolyRealRootInterval::Rational(x2),
                ) => {
                    if x2 < a1 {
                        all_roots.push((Interleave::Second, idx2));
                        idx2 += 1;
                        continue;
                    }
                    if b1 < x2 {
                        all_roots.push((Interleave::First, idx1));
                        idx1 += 1;
                        continue;
                    }
                }
                (
                    SquarefreePolyRealRootInterval::Real(a1, b1, _dir1),
                    SquarefreePolyRealRootInterval::Real(a2, b2, _dir2),
                ) => {
                    if b2 < a1 {
                        all_roots.push((Interleave::Second, idx2));
                        idx2 += 1;
                        continue;
                    }
                    if b1 < a2 {
                        all_roots.push((Interleave::First, idx1));
                        idx1 += 1;
                        continue;
                    }
                }
            }

            //refine and try again
            roots1.refine(idx1);
            roots2.refine(idx2);
        }

        debug_assert!(idx1 == roots1.intervals.len() || idx2 == roots2.intervals.len());

        while idx1 < roots1.intervals.len() {
            all_roots.push((Interleave::First, idx1));
            idx1 += 1;
        }

        while idx2 < roots2.intervals.len() {
            all_roots.push((Interleave::Second, idx2));
            idx2 += 1;
        }

        for r1 in &roots1.intervals {
            for r2 in &roots2.intervals {
                match (r1, r2) {
                    (
                        SquarefreePolyRealRootInterval::Rational(a),
                        SquarefreePolyRealRootInterval::Rational(x),
                    ) => {
                        debug_assert!(a != x);
                    }
                    (
                        SquarefreePolyRealRootInterval::Rational(a),
                        SquarefreePolyRealRootInterval::Real(x, y, _),
                    ) => {
                        debug_assert!(a < x || y < a);
                    }
                    (
                        SquarefreePolyRealRootInterval::Real(a, b, _),
                        SquarefreePolyRealRootInterval::Rational(x),
                    ) => {
                        debug_assert!(x < a || b < x);
                    }
                    (
                        SquarefreePolyRealRootInterval::Real(a, b, _),
                        SquarefreePolyRealRootInterval::Real(x, y, _),
                    ) => {
                        debug_assert!(b < x || y < a);
                    }
                }
            }
        }

        debug_assert_eq!(
            all_roots.len(),
            roots1.intervals.len() + roots2.intervals.len()
        );

        Ok(all_roots)
    }
}

//interval gen must yield intervals which converge to some real root of the polynomial
pub fn identify_real_root(
    poly: Polynomial<Integer>,
    mut interval_gen: impl Iterator<Item = (Rational, Rational)>,
) -> RealAlgebraic {
    let poly = poly.primitive_squarefree_part();
    let factored_poly = poly.factor().unwrap();
    let polys = Polynomial::<Integer>::structure()
        .into_factorizations()
        .to_powers(&factored_poly)
        .into_iter()
        .map(|(f, _k)| f)
        .collect::<Vec<_>>();
    //the root we are after is exactly one of the roots of the irreducible polynomials in polys
    //the task now is to refine alg1 and alg2 until the root is identified

    #[allow(clippy::redundant_closure_for_method_calls)]
    let mut root_groups: Vec<_> = polys
        .into_iter()
        .map(|p| p.all_real_roots_squarefree())
        .collect();

    //store indices of possible roots
    let mut possible = std::collections::HashSet::new();
    for i in 0..root_groups.len() {
        for j in 0..root_groups[i].intervals.len() {
            possible.insert((i, j));
        }
    }

    while possible.len() > 1 {
        let (ans_tight_a, ans_tight_b) = interval_gen.next().unwrap();
        //filter out roots which dont overlap with the known range for the sum root
        possible.retain(|(i, j)| match &root_groups[*i].intervals[*j] {
            SquarefreePolyRealRootInterval::Rational(x) => &ans_tight_a < x && x < &ans_tight_b,
            SquarefreePolyRealRootInterval::Real(ta, tb, _dir) => {
                ta < &ans_tight_b && &ans_tight_a < tb
            }
        });

        for (i, j) in &possible {
            root_groups[*i].refine(*j);
        }
    }
    assert_eq!(possible.len(), 1);
    let (i, j) = possible.into_iter().next().unwrap();
    let ans = root_groups
        .into_iter()
        .nth(i)
        .unwrap()
        .to_real_roots()
        .into_iter()
        .nth(j)
        .unwrap();
    #[cfg(debug_assertions)]
    ans.check_invariants().unwrap();
    ans
}

impl Polynomial<Integer> {
    fn sign_variations(&self) -> usize {
        //https://en.wikipedia.org/wiki/Descartes'_rule_of_signs
        //equals the number of strictly positive real roots modulo 2
        //and number of positive real roots is less than this number
        let nonzero_coeffs = self
            .coeffs()
            .into_iter()
            .filter(|c| c != &&Integer::zero())
            .collect::<Vec<_>>();
        let mut v = 0;
        for i in 0..nonzero_coeffs.len() - 1 {
            if (nonzero_coeffs[i] < &Integer::ZERO) != (nonzero_coeffs[i + 1] < &Integer::ZERO) {
                v += 1;
            }
        }
        v
    }

    //Collins and Akritas algorithm %https://en.wikipedia.org/wiki/Real-root_isolation
    fn isolate_real_roots_by_collin_akritas(&self) -> Vec<(Natural, usize, bool)> {
        //input: p(x), a square-free polynomial, such that p(0) p(1) ≠ 0, for which the roots in the interval [0, 1] are searched
        //output: a list of triples (c, k, h) representing isolating intervals of the form [c/2^k, (c+h)/2^k]
        debug_assert_ne!(self.evaluate(&Integer::zero()), Integer::zero());
        debug_assert_ne!(self.evaluate(&Integer::one()), Integer::zero());
        debug_assert_eq!(
            self.clone().primitive_squarefree_part().degree().unwrap(),
            self.degree().unwrap()
        );

        let mut l = vec![(Natural::from(0u8), 0, self.clone())];
        let mut isol = vec![];
        #[allow(clippy::manual_while_let_some)]
        while !l.is_empty() {
            let (c, k, mut q) = l.pop().unwrap();
            if q.evaluate(&Integer::from(0)) == Integer::from(0) {
                //q = q/x
                q = Self::div(&q, &Self::var()).unwrap();
                isol.push((c.clone(), k, false)); //rational root
            }
            let v = Self::compose(
                &q.reversed(),
                &Self::from_coeffs(vec![Integer::from(1), Integer::from(1)]),
            )
            .sign_variations();
            #[allow(clippy::comparison_chain)]
            if v == 1 {
                isol.push((c, k, true)); //root
            } else if v > 1 {
                //bisect
                //q_small(x) = 2^n q(x/2)
                let q_small = q.apply_map_with_powers(|(i, coeff)| {
                    coeff * Integer::from(Natural::TWO << (q.degree().unwrap() - i))
                });
                l.push((
                    (c.clone() << 1) + Natural::from(1u8),
                    k + 1,
                    Self::compose(
                        &q_small,
                        &Self::from_coeffs(vec![Integer::from(1), Integer::from(1)]),
                    ),
                ));
                l.push((c << 1, k + 1, q_small));
            }
        }
        isol
    }
}

impl Polynomial<Integer> {
    //isolate all real roots of a squarefree (no repeated roots) polynomial between a and b
    pub fn real_roots_squarefree(
        self,
        opt_a: Option<&Rational>,
        opt_b: Option<&Rational>,
        include_a: bool,
        include_b: bool,
    ) -> SquarefreePolyRealRoots {
        assert_ne!(self, Self::zero());
        //poly should be squarefree
        debug_assert_eq!(
            self.clone().primitive_squarefree_part().degree().unwrap(),
            self.degree().unwrap()
        );

        if let (Some(a), Some(b)) = (opt_a, opt_b) {
            assert!(a < b);
        }

        let d = self.degree().unwrap();
        if d == 0 {
            //constant polynomial has no roots
            SquarefreePolyRealRoots {
                poly_sqfr: self,
                intervals: vec![],
            }
        } else if d == 1 {
            //poly = a+bx
            //root = -a/b
            let root = unique_linear_root(&self);

            if {
                match opt_a {
                    Some(a) => match a.cmp(&root) {
                        std::cmp::Ordering::Less => true,
                        std::cmp::Ordering::Equal => include_a,
                        std::cmp::Ordering::Greater => false,
                    },
                    None => true,
                }
            } && {
                match opt_b {
                    Some(b) => match b.cmp(&root) {
                        std::cmp::Ordering::Greater => true,
                        std::cmp::Ordering::Equal => include_b,
                        std::cmp::Ordering::Less => false,
                    },
                    None => true,
                }
            } {
                SquarefreePolyRealRoots {
                    poly_sqfr: self,
                    intervals: vec![SquarefreePolyRealRootInterval::Rational(root)],
                }
            } else {
                SquarefreePolyRealRoots {
                    poly_sqfr: self,
                    intervals: vec![],
                }
            }
        } else {
            if opt_a.is_none() || opt_b.is_none() {
                //compute a bound M on the absolute value of any root
                //m = (Cauchy's bound + 1) https://captainblack.wordpress.com/2009/03/08/cauchys-upper-bound-for-the-roots-of-a-polynomial/
                let m = Rational::from(2)
                    + Rational::from_integers(
                        itertools::max((0..d).map(|i| self.coeff(i).abs().clone())).unwrap(),
                        self.coeff(d).abs().clone(),
                    );

                debug_assert!(m > Rational::ZERO);

                return match (opt_a, opt_b) {
                    (None, None) => {
                        let neg_m = -m.clone();
                        self.real_roots_squarefree(Some(&neg_m), Some(&m), include_a, include_b)
                    }
                    (None, Some(b_val)) => {
                        self.real_roots_squarefree(Some(&-m), Some(b_val), include_a, include_b)
                    }
                    (Some(a_val), None) => {
                        self.real_roots_squarefree(Some(a_val), Some(&m), include_a, include_b)
                    }
                    (Some(_a_val), Some(_b_val)) => {
                        panic!()
                    }
                };
            }
            let (a, b) = (opt_a.unwrap(), opt_b.unwrap());
            debug_assert!(a < b);

            //deal with end roots
            let mut poly_no_endroots = self.clone();
            let mut intervals = vec![];
            if evaluate_at_rational(&self, a) == Rational::from(0) {
                poly_no_endroots = Self::div(
                    &poly_no_endroots,
                    &Polynomial::from_coeffs(vec![-a.numerator(), a.denominator().into()]),
                )
                .unwrap();
                if include_a {
                    intervals.push(SquarefreePolyRealRootInterval::Rational(a.clone()));
                }
            }
            let mut do_add_b = false;
            if evaluate_at_rational(&self, b) == Rational::from(0) {
                poly_no_endroots = Self::div(
                    &poly_no_endroots,
                    &Polynomial::from_coeffs(vec![-b.numerator(), b.denominator().into()]),
                )
                .unwrap();
                if include_b {
                    do_add_b = true;
                }
            }

            debug_assert_ne!(
                evaluate_at_rational(&poly_no_endroots, a),
                Rational::from(0)
            );
            debug_assert_ne!(
                evaluate_at_rational(&poly_no_endroots, b),
                Rational::from(0)
            );

            //apply a transformation to p so that its roots in (a, b) are moved to roots in (0, 1)
            let (_, trans_poly) = Polynomial::compose(
                &poly_no_endroots.apply_map(|c| Rational::from(c)),
                &Polynomial::from_coeffs(vec![a.clone(), b.clone() - a.clone()]),
            )
            .factor_primitive_fof();

            'interval_loop: for (c, k, h) in trans_poly.isolate_real_roots_by_collin_akritas() {
                let d = Natural::from(1u8) << k;
                let mut interval_a = (b - a) * Rational::from_integers(c.clone(), d.clone()) + a;
                if h {
                    let mut interval_b =
                        (b - a) * Rational::from_integers(&c + Natural::from(1u8), d.clone()) + a;

                    //at the moment, interval_a and interval_b might be rational roots
                    //we need to strink them a little bit if so
                    if evaluate_at_rational(&self, &interval_a) == Rational::from(0)
                        || evaluate_at_rational(&self, &interval_b) == Rational::from(0)
                    {
                        let interval_m = (&interval_a + &interval_b) / Rational::from(2);
                        let mut shrunk_inerval_a = (&interval_a + &interval_m) / Rational::from(2);
                        let mut shrunk_inerval_b = (&interval_m + &interval_b) / Rational::from(2);
                        loop {
                            let shrunk_a_sign = evaluate_at_rational(&self, &shrunk_inerval_a)
                                .cmp(&Rational::from(0));
                            let shrunk_b_sign = evaluate_at_rational(&self, &shrunk_inerval_b)
                                .cmp(&Rational::from(0));
                            if shrunk_a_sign.is_eq() {
                                intervals.push(SquarefreePolyRealRootInterval::Rational(
                                    shrunk_inerval_a,
                                ));
                                continue 'interval_loop;
                            } else if shrunk_b_sign.is_eq() {
                                intervals.push(SquarefreePolyRealRootInterval::Rational(
                                    shrunk_inerval_b,
                                ));
                                continue 'interval_loop;
                            } else if shrunk_a_sign.is_ge() != shrunk_b_sign.is_ge() {
                                break;
                            }
                            shrunk_inerval_a =
                                (&interval_a + &shrunk_inerval_a) / Rational::from(2);
                            shrunk_inerval_b =
                                (&shrunk_inerval_b + &interval_b) / Rational::from(2);
                        }
                        interval_a = shrunk_inerval_a;
                        interval_b = shrunk_inerval_b;
                    }

                    let sign_a = evaluate_at_rational(&self, &interval_a) > Rational::from(0);
                    let sign_b = evaluate_at_rational(&self, &interval_b) > Rational::from(0);
                    debug_assert_ne!(evaluate_at_rational(&self, &interval_a), Rational::from(0));
                    debug_assert_ne!(evaluate_at_rational(&self, &interval_b), Rational::from(0));
                    debug_assert_ne!(sign_a, sign_b);
                    intervals.push(SquarefreePolyRealRootInterval::Real(
                        interval_a, interval_b, sign_b,
                    ));
                } else {
                    intervals.push(SquarefreePolyRealRootInterval::Rational(interval_a));
                }
            }

            if do_add_b {
                intervals.push(SquarefreePolyRealRootInterval::Rational(b.clone()));
            }

            let roots = SquarefreePolyRealRoots {
                poly_sqfr: self,
                intervals,
            };
            #[cfg(debug_assertions)]
            roots.check_invariants().unwrap();
            roots
        }
    }

    //isolate all real roots of the irreducible poly in the open interval (a, b)
    pub fn all_real_roots_squarefree(&self) -> SquarefreePolyRealRoots {
        self.clone().real_roots_squarefree(None, None, false, false)
    }

    //isolate all real roots of the irreducible poly in the open interval (a, b)
    pub fn real_roots_irreducible(
        &self,
        opt_a: Option<&Rational>,
        opt_b: Option<&Rational>,
        include_a: bool,
        include_b: bool,
    ) -> Vec<RealAlgebraic> {
        debug_assert!(self.is_irreducible());
        self.clone()
            .real_roots_squarefree(opt_a, opt_b, include_a, include_b)
            .to_real_roots()
    }

    //get the real roots with multiplicity of poly
    pub fn real_roots(
        &self,
        a: Option<&Rational>,
        b: Option<&Rational>,
        include_a: bool,
        include_b: bool,
    ) -> Vec<RealAlgebraic> {
        assert_ne!(self, &Self::zero());
        let factors = self.factor().unwrap();
        let mut roots = vec![];
        for (factor, k) in Polynomial::<Integer>::structure()
            .factorizations()
            .to_powers(&factors)
        {
            for root in factor.real_roots_irreducible(a, b, include_a, include_b) {
                let mut i = Natural::from(0u8);
                while &i < k {
                    roots.push(root.clone());
                    i += Natural::from(1u8);
                }
            }
        }
        roots
    }

    pub fn all_real_roots(&self) -> Vec<RealAlgebraic> {
        self.real_roots(None, None, false, false)
    }
}

impl Polynomial<Rational> {
    pub fn all_real_roots(&self) -> Vec<RealAlgebraic> {
        assert_ne!(self, &Self::zero());
        self.primitive_part_fof().all_real_roots()
    }
}

pub fn nth_root(x: &RealAlgebraic, n: usize) -> Result<RealAlgebraic, ()> {
    if n == 0 {
        panic!()
    } else if n == 1 {
        Ok(x.clone())
    } else {
        match x.cmp(&RealAlgebraic::zero()) {
            std::cmp::Ordering::Less => Err(()),
            std::cmp::Ordering::Equal => Ok(RealAlgebraic::zero()),
            std::cmp::Ordering::Greater => {
                let poly = match x {
                    RealAlgebraic::Rational(rat) => {
                        Polynomial::from_coeffs(vec![-rat.numerator(), rat.denominator().into()])
                    }
                    RealAlgebraic::Real(real) => real.poly.clone(),
                };
                let mut coeffs = vec![];
                for (i, c) in poly.into_coeffs().into_iter().enumerate() {
                    if i != 0 {
                        for _j in 0..(n - 1) {
                            coeffs.push(Integer::ZERO);
                        }
                    }
                    coeffs.push(c);
                }
                let nthroot_poly = Polynomial::from_coeffs(coeffs).primitive_squarefree_part();
                // println!("nthroot_poly = {:?}", nthroot_poly);
                let possible_nthroots = nthroot_poly.all_real_roots_squarefree();
                //search through the intervals starting from the largest.
                //if it is positive, check that the nth power of the bounds surrounds self
                //otherwise, the first one that does not have both bounds positive must be the nth root

                // println!("{:?}", possible_nthroots);

                for (idx, interval) in possible_nthroots.intervals.iter().enumerate().rev() {
                    match interval {
                        SquarefreePolyRealRootInterval::Rational(rat) => {
                            debug_assert!(rat > &Rational::ZERO);
                            match x {
                                RealAlgebraic::Rational(self_rat) => {
                                    debug_assert_eq!(&rat.nat_pow(&Natural::from(n)), self_rat);
                                }
                                RealAlgebraic::Real(_) => debug_assert!(false),
                            }
                            return Ok(possible_nthroots.get_real_root(idx));
                        }
                        SquarefreePolyRealRootInterval::Real(a, b, _dir) => {
                            debug_assert!(a < b);
                            if &Rational::ZERO < a {
                                let a_pow = a.nat_pow(&Natural::from(n));
                                let b_pow = b.nat_pow(&Natural::from(n));
                                if &RealAlgebraic::Rational(a_pow) <= x
                                    && x <= &RealAlgebraic::Rational(b_pow)
                                {
                                    return Ok(possible_nthroots.get_real_root(idx));
                                }
                            } else {
                                return Ok(possible_nthroots.get_real_root(idx));
                            }
                        }
                    }
                }

                unreachable!()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_squarefree_polynomial_real_root_isolation() {
        let f = Polynomial::product(vec![
            &Polynomial::from_coeffs(vec![
                Integer::from(-2),
                Integer::from(-4),
                Integer::from(-2),
            ]),
            &Polynomial::from_coeffs(vec![Integer::from(6), Integer::from(0), Integer::from(-3)]),
            &Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-3), Integer::from(1)]),
            &Polynomial::from_coeffs(vec![
                Integer::from(2),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ]),
            &Polynomial::from_coeffs(vec![
                Integer::from(1),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ]),
            &Polynomial::from_coeffs(vec![
                Integer::from(-1),
                Integer::from(12),
                Integer::from(-4),
                Integer::from(-15),
                Integer::from(5),
                Integer::from(3),
                Integer::from(-1),
            ]),
        ]);
        let f = f.primitive_squarefree_part();
        //f is a squarefree polynomial with lots of roots
        println!("f = {:?}", f);
        let intervals = Polynomial::real_roots_squarefree(f, None, None, false, false);
        println!("intervals = {:?}", &intervals);
        intervals.check_invariants().unwrap();

        let f =
            Polynomial::from_coeffs(vec![Integer::from(1), Integer::from(-3), Integer::from(1)]);
        println!("f = {:?}", f);
        let mut intervals = Polynomial::real_roots_squarefree(f, None, None, false, false);
        intervals.check_invariants().unwrap();
        intervals.clone().to_real_roots();
        for root in intervals.clone().to_real_roots() {
            println!("root = {:?}", root);
            root.check_invariants().unwrap();
        }
        println!("refine");
        for _i in 0..10 {
            intervals.refine_all();
        }
        for root in intervals.clone().to_real_roots() {
            println!("root = {:?}", root);
            root.check_invariants().unwrap();
        }

        let f = Polynomial::product(vec![
            &Polynomial::from_coeffs(vec![Integer::from(-1), Integer::from(1)]),
            &Polynomial::from_coeffs(vec![Integer::from(-2), Integer::from(1)]),
            &Polynomial::from_coeffs(vec![Integer::from(-3), Integer::from(1)]),
            &Polynomial::from_coeffs(vec![Integer::from(-4), Integer::from(1)]),
        ]);
        assert_eq!(
            f.clone()
                .real_roots_squarefree(
                    Some(&Rational::from(1)),
                    Some(&Rational::from(4)),
                    false,
                    false
                )
                .intervals
                .len(),
            2
        );
        assert_eq!(
            f.clone()
                .real_roots_squarefree(
                    Some(&Rational::from(1)),
                    Some(&Rational::from(4)),
                    true,
                    false
                )
                .intervals
                .len(),
            3
        );
        assert_eq!(
            f.clone()
                .real_roots_squarefree(
                    Some(&Rational::from(1)),
                    Some(&Rational::from(4)),
                    false,
                    true
                )
                .intervals
                .len(),
            3
        );
        assert_eq!(
            f.clone()
                .real_roots_squarefree(
                    Some(&Rational::from(1)),
                    Some(&Rational::from(4)),
                    true,
                    true
                )
                .intervals
                .len(),
            4
        );
    }

    #[test]
    fn test_real_roots_squarefree() {
        //a test where real_roots_squarefree find rational roots while refining the bounds on a real root such that no rational root lies on the end points
        let poly = Polynomial::from_coeffs(vec![
            Integer::from(0),
            Integer::from(9),
            Integer::from(0),
            Integer::from(-4),
        ]);
        let roots = poly.real_roots_squarefree(
            Some(&Rational::from(-2)),
            Some(&Rational::from(2)),
            false,
            true,
        );
        println!("{:?}", roots);
        roots.check_invariants().unwrap();
    }

    #[test]
    fn test_real_root_irreducible_count() {
        assert_eq!(
            Polynomial::from_coeffs(vec![
                Integer::from(3),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1)
            ])
            .real_roots_irreducible(None, None, false, false)
            .len(),
            1
        );
        assert_eq!(
            Polynomial::from_coeffs(vec![
                Integer::from(1),
                Integer::from(-3),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1)
            ])
            .real_roots_irreducible(None, None, false, false)
            .len(),
            3
        );
    }
}
