use std::fmt::Display;
use std::rc::Rc;

use itertools::Itertools;

use crate::linear::matrix::*;

use super::super::structure::structure::*;
use algebraeon_sets::structure::*;

use algebraeon_nzq::natural::*;

#[derive(Debug, Clone)]
pub struct Polynomial<Set> {
    //vec![c0, c1, c2, c3, ..., cn] represents the polynomial c0 + c1*x + c2*x^2 + c3*x^3 + ... + cn * x^n
    //if non-empty, the last item must not be zero
    coeffs: Vec<Set>,
}

impl<Set: std::hash::Hash> std::hash::Hash for Polynomial<Set> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.coeffs.hash(state);
    }
}

impl<Set> Polynomial<Set> {
    pub fn coeffs(&self) -> Vec<&Set> {
        self.coeffs.iter().collect()
    }

    pub fn into_coeffs(self) -> Vec<Set> {
        self.coeffs
    }

    pub fn from_coeffs(coeffs: Vec<impl Into<Set>>) -> Self {
        Self {
            coeffs: coeffs.into_iter().map(|x| x.into()).collect(),
        }
    }

    pub fn constant(x: Set) -> Self {
        Self::from_coeffs(vec![x])
    }

    pub fn apply_map<ImgSet>(&self, f: impl Fn(&Set) -> ImgSet) -> Polynomial<ImgSet> {
        Polynomial::from_coeffs(self.coeffs.iter().map(f).collect())
    }

    pub fn apply_map_with_powers<ImgSet>(
        &self,
        f: impl Fn((usize, &Set)) -> ImgSet,
    ) -> Polynomial<ImgSet> {
        Polynomial::from_coeffs(self.coeffs.iter().enumerate().map(f).collect())
    }
}

#[derive(Debug, Clone)]
pub struct PolynomialStructure<RS: SemiRingStructure> {
    coeff_ring_zero: RS::Set, //so that we can return a refernece to zero when getting polynomial coefficients out of range
    coeff_ring: Rc<RS>,
}

impl<RS: SemiRingStructure> PolynomialStructure<RS> {
    pub fn coeff_ring(&self) -> Rc<RS> {
        self.coeff_ring.clone()
    }
}

impl<RS: SemiRingStructure> Structure for PolynomialStructure<RS> {
    type Set = Polynomial<RS::Set>;
}

impl<RS: SemiRingStructure> PartialEq for PolynomialStructure<RS> {
    fn eq(&self, other: &Self) -> bool {
        self.coeff_ring == other.coeff_ring
    }
}

impl<RS: SemiRingStructure> Eq for PolynomialStructure<RS> {}

impl<RS: SemiRingStructure + ToStringStructure> ToStringStructure for PolynomialStructure<RS> {
    fn to_string(&self, elem: &Self::Set) -> String {
        if self.num_coeffs(elem) == 0 {
            "0".into()
        } else {
            let mut s = String::new();
            let mut first = true;
            for (k, c) in elem.coeffs.iter().enumerate() {
                if !self.coeff_ring.is_zero(c) {
                    if self.coeff_ring.equal(c, &self.coeff_ring.one()) {
                        if k == 0 {
                            s += "1";
                        } else {
                            s += "+";
                        }
                    } else {
                        if !first {
                            s += "+";
                        }
                        s += "(";
                        s += &self.coeff_ring.to_string(c);
                        s += ")";
                    }
                    if k == 0 {
                    } else if k == 1 {
                        s += "λ";
                    } else {
                        s += "λ";
                        s += "^";
                        s += &k.to_string();
                    }
                    first = false;
                }
            }
            s
        }
    }
}

impl<RS: SemiRingStructure> PartialEqStructure for PolynomialStructure<RS> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        for i in 0..std::cmp::max(a.coeffs.len(), b.coeffs.len()) {
            if !self.coeff_ring.equal(self.coeff(a, i), self.coeff(b, i)) {
                return false;
            }
        }
        true
    }
}

impl<RS: SemiRingStructure> EqStructure for PolynomialStructure<RS> {}

impl<RS: SemiRingStructure> SemiRingStructure for PolynomialStructure<RS> {
    fn zero(&self) -> Self::Set {
        Polynomial { coeffs: vec![] }
    }

    fn one(&self) -> Self::Set {
        Polynomial::from_coeffs(vec![self.coeff_ring.one()])
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.reduce_poly(Polynomial::from_coeffs(
            (0..std::cmp::max(a.coeffs.len(), b.coeffs.len()))
                .map(|i| self.coeff_ring.add(self.coeff(a, i), self.coeff(b, i)))
                .collect(),
        ))
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let mut coeffs = Vec::with_capacity(a.coeffs.len() + b.coeffs.len());
        for _k in 0..a.coeffs.len() + b.coeffs.len() {
            coeffs.push(self.coeff_ring.zero());
        }
        for i in 0..a.coeffs.len() {
            for j in 0..b.coeffs.len() {
                self.coeff_ring.add_mut(
                    &mut coeffs[i + j],
                    &self.coeff_ring.mul(&a.coeffs[i], &b.coeffs[j]),
                );
            }
        }
        self.reduce_poly(Polynomial::from_coeffs(coeffs))
    }
}

impl<RS: RingStructure> RingStructure for PolynomialStructure<RS> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        Polynomial::from_coeffs(
            a.coeffs()
                .into_iter()
                .map(|c| self.coeff_ring.neg(&c))
                .collect(),
        )
    }
}

impl<RS: SemiRingStructure> PolynomialStructure<RS> {
    pub fn new(coeff_ring: Rc<RS>) -> Self {
        Self {
            coeff_ring_zero: coeff_ring.zero(),
            coeff_ring,
        }
    }

    pub fn reduce_poly(&self, mut a: Polynomial<RS::Set>) -> Polynomial<RS::Set> {
        loop {
            if a.coeffs.len() == 0 {
                break;
            } else {
                if self
                    .coeff_ring
                    .equal(&a.coeffs.last().unwrap(), &self.coeff_ring.zero())
                {
                    a.coeffs.pop();
                } else {
                    break;
                }
            }
        }
        a
    }

    pub fn var(&self) -> Polynomial<RS::Set> {
        Polynomial::from_coeffs(vec![self.coeff_ring.zero(), self.coeff_ring.one()])
    }

    pub fn var_pow(&self, n: usize) -> Polynomial<RS::Set> {
        Polynomial::from_coeffs(
            (0..n + 1)
                .map(|i| {
                    if i < n {
                        self.coeff_ring.zero()
                    } else {
                        self.coeff_ring.one()
                    }
                })
                .collect(),
        )
    }

    fn constant_var_pow(&self, x: RS::Set, n: usize) -> Polynomial<RS::Set> {
        Polynomial::from_coeffs(
            (0..n + 1)
                .map(|i| {
                    if i < n {
                        self.coeff_ring.zero()
                    } else {
                        x.clone()
                    }
                })
                .collect(),
        )
    }

    pub fn coeff<'a>(&'a self, a: &'a Polynomial<RS::Set>, i: usize) -> &'a RS::Set {
        match a.coeffs.get(i) {
            Some(c) => c,
            None => &self.coeff_ring_zero,
        }
    }

    pub fn leading_coeff<'a>(&self, a: &'a Polynomial<RS::Set>) -> Option<&'a RS::Set> {
        Some(a.coeffs.get(self.degree(a)?).unwrap())
    }

    pub fn evaluate(&self, p: &Polynomial<RS::Set>, x: &RS::Set) -> RS::Set {
        // f(x) = a + bx + cx^2 + dx^3
        // evaluate as f(x) = a + x(b + x(c + x(d)))
        let mut y = self.coeff_ring.zero();
        for c in p.coeffs().into_iter().rev() {
            self.coeff_ring.mul_mut(&mut y, x);
            self.coeff_ring.add_mut(&mut y, c)
        }
        y
    }

    //find p(q(x))
    pub fn compose(&self, p: &Polynomial<RS::Set>, q: &Polynomial<RS::Set>) -> Polynomial<RS::Set> {
        PolynomialStructure::new(Rc::new(self.clone()))
            .evaluate(&p.apply_map(|c| Polynomial::constant(c.clone())), q)
    }

    //if n = deg(p)
    //return x^n * p(1/x)
    pub fn reversed(&self, p: &Polynomial<RS::Set>) -> Polynomial<RS::Set> {
        Polynomial::from_coeffs(p.coeffs.clone().into_iter().rev().collect())
    }

    pub fn mul_var_pow(&self, p: &Polynomial<RS::Set>, n: usize) -> Polynomial<RS::Set> {
        let mut coeffs = vec![];
        for _i in 0..n {
            coeffs.push(self.coeff_ring.zero());
        }
        for c in &p.coeffs {
            coeffs.push(c.clone());
        }
        Polynomial { coeffs }
    }

    pub fn mul_scalar(&self, p: &Polynomial<RS::Set>, x: &RS::Set) -> Polynomial<RS::Set> {
        self.reduce_poly(Polynomial::from_coeffs(
            p.coeffs.iter().map(|c| self.coeff_ring.mul(c, x)).collect(),
        ))
    }

    pub fn num_coeffs(&self, p: &Polynomial<RS::Set>) -> usize {
        match self.degree(p) {
            Some(n) => n + 1,
            None => 0,
        }
    }

    //zero -> None
    //const -> 0
    //linear -> 1
    //quadratic -> 2
    //...
    pub fn degree(&self, p: &Polynomial<RS::Set>) -> Option<usize> {
        //the polynomial representation might not be reduced
        for i in (0..p.coeffs.len()).rev() {
            if !self.coeff_ring.is_zero(&p.coeffs[i]) {
                return Some(i);
            }
        }
        None
    }

    pub fn as_constant(&self, p: &Polynomial<RS::Set>) -> Option<RS::Set> {
        if self.num_coeffs(p) == 0 {
            Some(self.coeff_ring.zero())
        } else if self.num_coeffs(p) == 1 {
            Some(p.coeffs[0].clone())
        } else {
            None
        }
    }

    pub fn is_monic(&self, p: &Polynomial<RS::Set>) -> bool {
        match self.leading_coeff(p) {
            Some(lc) => self.coeff_ring.equal(lc, &self.coeff_ring.one()),
            None => false,
        }
    }

    pub fn derivative(&self, mut p: Polynomial<RS::Set>) -> Polynomial<RS::Set> {
        p = self.reduce_poly(p);
        if self.num_coeffs(&p) > 0 {
            for i in 0..self.num_coeffs(&p) - 1 {
                p.coeffs[i] = p.coeffs[i + 1].clone();
                self.coeff_ring.mul_mut(
                    &mut p.coeffs[i],
                    &self.coeff_ring.from_nat(Natural::from(i + 1)),
                );
            }
            p.coeffs.pop();
        }
        p
    }
}

// #[derive(Debug)]
// pub struct PseudoRemainderSubresultantSequence<RS: IntegralDomainStructure> {
//     poly_ring: PolynomialStructure<RS>,
//     n: usize,
//     rs: Vec<Polynomial<RS::Set>>,
//     subresultants: Vec<Polynomial<RS::Set>>,
//     smallest_nonzero_subresultant_idx: usize,
// }

// impl<RS: IntegralDomainStructure> PseudoRemainderSubresultantSequence<RS> {
//     fn new(poly_ring: PolynomialStructure<RS>, n: usize) -> Self {
//         let subresultants = (0..n).map(|_| poly_ring.zero()).collect();
//         Self {
//             poly_ring,
//             n,
//             rs: vec![],
//             subresultants,
//             smallest_nonzero_subresultant_idx: n,
//         }
//     }

//     fn push(&mut self, r: Polynomial<RS::Set>) {
//         let deg = self.poly_ring.degree(&r).unwrap();
//         match self.rs.last() {
//             Some(prev_r) => {
//                 let prev_deg = self.poly_ring.degree(prev_r).unwrap();
//                 debug_assert!(deg < prev_deg);
//                 self.subresultants[prev_deg - 1] = r.clone();
//                 self.smallest_nonzero_subresultant_idx = prev_deg - 1;
//                 if deg < prev_deg - 1 {
//                     self.subresultants[deg] = self.poly_ring.mul(
//                         &self
//                             .poly_ring
//                             .nat_pow(&r, &Natural::from(prev_deg - deg - 1)),
//                         &r,
//                     );
//                     self.smallest_nonzero_subresultant_idx = deg;
//                 }
//             }
//             None => {
//                 debug_assert_eq!(self.poly_ring.degree(&r).unwrap(), self.n);
//             }
//         };
//         self.rs.push(r);
//     }

//     pub fn resultant(&self) -> RS::Set {
//         debug_assert_ne!(self.rs.len(), 0); // Should have added something by the time it comes to using it
//         self.poly_ring.as_constant(&self.subresultants[0]).unwrap()
//     }

//     pub fn subresultant_gcd(&self) -> &Polynomial<RS::Set> {
//         debug_assert_ne!(self.rs.len(), 0); // Should have added something by the time it comes to using it
//         &self.subresultants[self.smallest_nonzero_subresultant_idx]
//     }

//     pub fn into_subresultant_gcd(self) -> Polynomial<RS::Set> {
//         debug_assert_ne!(self.rs.len(), 0); // Should have added something by the time it comes to using it
//         self.subresultants
//             .into_iter()
//             .nth(self.smallest_nonzero_subresultant_idx)
//             .unwrap()
//     }
// }

impl<RS: IntegralDomainStructure> PolynomialStructure<RS> {
    pub fn try_quorem(
        &self,
        a: &Polynomial<RS::Set>,
        b: &Polynomial<RS::Set>,
    ) -> Result<(Polynomial<RS::Set>, Polynomial<RS::Set>), RingDivisionError> {
        //try to find q such that q*b == a
        // a0 + a1*x + a2*x^2 + ... + am*x^m = (q0 + q1*x + q2*x^2 + ... + qk*x^k) * (b0 + b1*x + b2*x^2 + ... + bn*x^n)
        // 1 + x + x^2 + x^3 + x^4 + x^5 = (?1 + ?x + ?x^2) * (1 + x + x^2 + x^3)      m=6 k=3 n=4

        let mut a = a.clone();

        let m = self.num_coeffs(&a);
        let n = self.num_coeffs(b);
        if n == 0 {
            Err(RingDivisionError::DivideByZero)
        } else if m < n {
            Ok((self.zero(), a))
        } else {
            let k = m - n + 1;
            let mut q_coeffs = (0..k).map(|_i| self.coeff_ring.zero()).collect_vec();
            for i in (0..k).rev() {
                //a[i+n-1] = q[i] * b[n-1]
                match self
                    .coeff_ring
                    .div(self.coeff(&a, i + n - 1), &self.coeff(b, n - 1))
                {
                    Ok(qc) => {
                        //a -= qc*x^i*b
                        self.add_mut(
                            &mut a,
                            &self.neg(&self.mul_var_pow(&self.mul_scalar(&b, &qc), i)),
                        );
                        q_coeffs[i] = qc;
                    }
                    Err(RingDivisionError::NotDivisible) => {
                        return Err(RingDivisionError::NotDivisible);
                    }
                    Err(RingDivisionError::DivideByZero) => panic!(),
                }
            }
            Ok((Polynomial::from_coeffs(q_coeffs), a))
        }
    }

    pub fn div_impl(
        &self,
        a: &Polynomial<RS::Set>,
        b: &Polynomial<RS::Set>,
    ) -> Result<Polynomial<RS::Set>, RingDivisionError> {
        match self.try_quorem(&a, &b) {
            Ok((q, r)) => {
                debug_assert!(self.equal(&self.add(&self.mul(&q, &b), &r), &a));
                if self.is_zero(&r) {
                    Ok(q)
                } else {
                    Err(RingDivisionError::NotDivisible)
                }
            }
            Err(RingDivisionError::NotDivisible) => Err(RingDivisionError::NotDivisible),
            Err(RingDivisionError::DivideByZero) => Err(RingDivisionError::DivideByZero),
        }
    }

    //None if b = 0
    //error if deg(a) < deg(b)
    pub fn pseudorem(
        &self,
        mut a: Polynomial<RS::Set>,
        b: &Polynomial<RS::Set>,
    ) -> Option<Result<Polynomial<RS::Set>, &'static str>> {
        let m = self.num_coeffs(&a);
        let n = self.num_coeffs(b);

        if n == 0 {
            None
        } else if m < n {
            Some(Err("Should have deg(a) >= deg(b) for pseudo remainder"))
        } else {
            self.mul_mut(
                &mut a,
                &Polynomial::constant(
                    self.coeff_ring
                        .nat_pow(&self.coeff(b, n - 1), &Natural::from(m - n + 1)),
                ),
            );

            match self.try_quorem(&a, b) {
                Ok((_q, r)) => Some(Ok(r)),
                Err(_) => panic!(),
            }
        }
    }

    /// Output: (r, s) where
    ///     r is the pseudo remainder sequence
    ///     s is the scalar subresultants
    pub fn pseudo_remainder_subresultant_sequence(
        &self,
        mut a: Polynomial<RS::Set>,
        mut b: Polynomial<RS::Set>,
    ) -> (Vec<Polynomial<RS::Set>>, Vec<RS::Set>) {
        match (self.degree(&a), self.degree(&b)) {
            (None, None) => (vec![], vec![]),
            (None, Some(_)) => (vec![b], vec![self.coeff_ring().one()]),
            (Some(_), None) => (vec![a], vec![self.coeff_ring().one()]),
            (Some(mut a_deg), Some(mut b_deg)) => {
                if a_deg < b_deg {
                    (a, b, a_deg, b_deg) = (b, a, b_deg, a_deg);
                }
                debug_assert!(a_deg >= b_deg);
                let mut prs = vec![a.clone(), b.clone()];
                let mut diff_deg = a_deg - b_deg;
                let mut beta = self.coeff_ring().nat_pow(
                    &self.coeff_ring().neg(&self.coeff_ring().one()),
                    &Natural::from(diff_deg + 1),
                );
                let mut r = self.mul_scalar(&self.pseudorem(a, &b).unwrap().unwrap(), &beta);
                let mut lc_b = self.leading_coeff(&b).unwrap().clone();
                let mut gamma = self.coeff_ring().nat_pow(&lc_b, &Natural::from(diff_deg));
                let mut ssres = vec![self.coeff_ring().one(), gamma.clone()];
                gamma = self.coeff_ring.neg(&gamma);
                loop {
                    match self.degree(&r) {
                        Some(r_deg) => {
                            prs.push(r.clone());
                            (a, b, b_deg, diff_deg) = (b, r, r_deg, b_deg - r_deg);
                            beta = self.coeff_ring.mul(
                                &self.coeff_ring.neg(&lc_b),
                                &self.coeff_ring().nat_pow(&gamma, &Natural::from(diff_deg)),
                            );
                            r = self
                                .div(
                                    &self.pseudorem(a, &b).unwrap().unwrap(),
                                    &Polynomial::constant(beta),
                                )
                                .unwrap();
                            lc_b = self.leading_coeff(&b).unwrap().clone();
                            gamma = if diff_deg > 1 {
                                self.coeff_ring()
                                    .div(
                                        &self.coeff_ring().nat_pow(
                                            &self.coeff_ring().neg(&lc_b),
                                            &Natural::from(diff_deg),
                                        ),
                                        &self
                                            .coeff_ring()
                                            .nat_pow(&gamma, &Natural::from(diff_deg - 1)),
                                    )
                                    .unwrap()
                            } else {
                                self.coeff_ring().neg(&lc_b)
                            };
                            ssres.push(self.coeff_ring().neg(&gamma));
                        }
                        None => {
                            debug_assert!(self.is_zero(&r));
                            break;
                        }
                    }
                }
                (prs, ssres)
            }
        }
    }

    // efficiently compute the gcd of a and b up to scalar multipication using pseudo-remainder subresultant sequence
    pub fn subresultant_gcd(
        &self,
        a: Polynomial<RS::Set>,
        b: Polynomial<RS::Set>,
    ) -> Polynomial<RS::Set> {
        // The GCD (up to scalar mul) is the last subresultant pseudoremainder
        let (mut prs, _) = self.pseudo_remainder_subresultant_sequence(a, b);
        prs.pop().unwrap()
    }

    pub fn resultant(&self, a: Polynomial<RS::Set>, b: Polynomial<RS::Set>) -> RS::Set {
        if self.is_zero(&a) || self.is_zero(&b) {
            self.coeff_ring().zero()
        } else {
            let (mut prs, mut ssres) = self.pseudo_remainder_subresultant_sequence(a, b);
            if self.degree(&prs.pop().unwrap()).unwrap() > 0 {
                self.coeff_ring().zero()
            } else {
                ssres.pop().unwrap()
            }
        }
    }

    pub fn is_squarefree(&self, p: &Polynomial<RS::Set>) -> bool {
        let dp = self.derivative(p.clone());
        self.degree(&self.subresultant_gcd(p.clone(), dp)).unwrap() == 0
    }

    pub fn discriminant(&self, p: Polynomial<RS::Set>) -> Result<RS::Set, &'static str> {
        match self.degree(&p) {
            Some(n) => {
                if n == 0 {
                    Err("Discriminant of a constant polynomial is undefined.")
                } else {
                    let an = self.coeff(&p, n).clone(); // leading coeff
                    let dp = self.derivative(p.clone());
                    let disc = self.coeff_ring().div(&self.resultant(p, dp), &an).unwrap();
                    // multiply by (-1)^{n(n+1)/2}
                    match n % 4 {
                        0 | 1 => Ok(disc),
                        2 | 3 => Ok(self.coeff_ring().neg(&disc)),
                        _ => unreachable!(),
                    }
                }
            }
            None => Err("Discriminant of zero polynomial is undefined."),
        }
    }
}

impl<RS: IntegralDomainStructure> IntegralDomainStructure for PolynomialStructure<RS> {
    fn div(&self, a: &Self::Set, b: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        self.div_impl(a, b)
    }
}

// impl<R: UniqueFactorizationDomain> UniqueFactorizationDomain for Polynomial<R> {}

impl<RS: GreatestCommonDivisorStructure> PolynomialStructure<RS> {
    pub fn factor_primitive(
        &self,
        mut p: Polynomial<RS::Set>,
    ) -> Option<(RS::Set, Polynomial<RS::Set>)> {
        if self.is_zero(&p) {
            None
        } else {
            let g = self.coeff_ring.gcd_list(p.coeffs.iter().collect());
            for i in 0..p.coeffs.len() {
                p.coeffs[i] = self.coeff_ring.div(&p.coeffs[i], &g).unwrap()
            }
            Some((g, p))
        }
    }

    pub fn is_primitive(&self, p: Polynomial<RS::Set>) -> bool {
        match self.factor_primitive(p) {
            Some((unit, _)) => self.coeff_ring().is_unit(&unit),
            None => false,
        }
    }

    pub fn primitive_part(&self, p: Polynomial<RS::Set>) -> Option<Polynomial<RS::Set>> {
        match self.factor_primitive(p) {
            Some((_unit, prim)) => Some(prim),
            None => None,
        }
    }

    pub fn gcd_by_primitive_subresultant(
        &self,
        a: Polynomial<RS::Set>,
        b: Polynomial<RS::Set>,
    ) -> Polynomial<RS::Set> {
        if self.is_zero(&a) {
            b
        } else if self.is_zero(&b) {
            a
        } else {
            let (a_unit, a_prim) = self.factor_primitive(a).unwrap();
            let (b_unit, b_prim) = self.factor_primitive(b).unwrap();
            let g_unit = self.coeff_ring.gcd(&a_unit, &b_unit);
            let g_prim = self
                .factor_primitive(self.subresultant_gcd(a_prim, b_prim))
                .unwrap()
                .1;
            let g = self.mul(&Polynomial::constant(g_unit), &g_prim);
            self.primitive_part(g).unwrap()
        }
    }
}

impl<FS: FieldStructure> GreatestCommonDivisorStructure for PolynomialStructure<FS> {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        self.euclidean_gcd(x.clone(), y.clone())
    }
}

impl<FS: FieldStructure> BezoutDomainStructure for PolynomialStructure<FS> {
    fn xgcd(&self, x: &Self::Set, y: &Self::Set) -> (Self::Set, Self::Set, Self::Set) {
        self.euclidean_xgcd(x.clone(), y.clone())
    }
}

impl<RS: GreatestCommonDivisorStructure + CharZeroStructure> PolynomialStructure<RS> {
    pub fn primitive_squarefree_part(&self, f: Polynomial<RS::Set>) -> Polynomial<RS::Set> {
        if self.is_zero(&f) {
            f
        } else {
            let g = self.subresultant_gcd(f.clone(), self.derivative(f.clone()));
            let (_c, g_prim) = self.factor_primitive(g).unwrap();
            let (_c, f_prim) = self.factor_primitive(f).unwrap();
            let f_prim_sqfree = self.div(&f_prim, &g_prim).unwrap();
            f_prim_sqfree
        }
    }
}

impl<RS: FavoriteAssociateStructure + IntegralDomainStructure> FavoriteAssociateStructure
    for PolynomialStructure<RS>
{
    fn factor_fav_assoc(
        &self,
        a: &Polynomial<RS::Set>,
    ) -> (Polynomial<RS::Set>, Polynomial<RS::Set>) {
        if self.is_zero(&a) {
            (self.one(), self.zero())
        } else {
            let mut a = a.clone();
            let (u, _c) = self
                .coeff_ring
                .factor_fav_assoc(&a.coeffs[self.num_coeffs(&a) - 1]);
            for i in 0..a.coeffs.len() {
                a.coeffs[i] = self.coeff_ring.div(&a.coeffs[i], &u).unwrap()
            }
            (Polynomial::constant(u), a.clone())
        }
    }
}

impl<RS: CharZeroStructure> CharZeroStructure for PolynomialStructure<RS> {}

impl<RS: IntegralDomainStructure + FiniteUnitsStructure> FiniteUnitsStructure
    for PolynomialStructure<RS>
{
    fn all_units(&self) -> Vec<Self::Set> {
        self.coeff_ring
            .all_units()
            .into_iter()
            .map(|u| Polynomial::constant(u))
            .collect()
    }
}

// pub trait InterpolatablePolynomials: ComRS {
//     fn interpolate(points: &Vec<(Self::ElemT, Self::ElemT)>) -> Option<Polynomial<Self>>;
// }

impl<FS: FieldStructure> EuclideanDivisionStructure for PolynomialStructure<FS> {
    fn norm(&self, elem: &Self::Set) -> Option<Natural> {
        Some(Natural::from(self.degree(elem)?))
    }

    fn quorem(&self, a: &Self::Set, b: &Self::Set) -> Option<(Self::Set, Self::Set)> {
        match self.try_quorem(a, b) {
            Ok((q, r)) => Some((q, r)),
            Err(RingDivisionError::NotDivisible) => panic!(),
            Err(RingDivisionError::DivideByZero) => None,
        }
    }
}

impl<RS: IntegralDomainStructure> PolynomialStructure<RS> {
    pub fn interpolate_by_lagrange_basis(
        &self,
        points: &Vec<(RS::Set, RS::Set)>,
    ) -> Option<Polynomial<RS::Set>> {
        /*
        points should be a list of pairs (xi, yi) where the xi are distinct

        find f such that f(x1) = y1, f(x2) = y2, f(x3) = y3
        find f1(x), f2(x), f3(x) such that fi(xi)=1 and fi(xj)=0 if i!=j

            (x-x2) (x-x3)         (x-x1) (x-x3)         (x-x1) (x-x2)
        f1= --------------    f2= --------------    f3= --------------
            (x1-x2)(x1-x3)        (x2-x1)(x2-x3)        (x3-x1)(x3-x2)

        so f(x) = y1f1(x) + y2f2(x) + y3f3(x)
                = (   y1 (x-x2)(x-x3) (x2-x3)
                    + y2 (x1-x)(x-x3) (x1-x3)
                    + y3 (x1-x)(x2-x) (x1-x2) )
                  / (x1-x2)(x1-x3)(x2-x3)

         */
        let mut numerator = self.zero();
        for i in 0..points.len() {
            let (_xi, yi) = &points[i];
            let mut term = Polynomial::constant(yi.clone());

            for j in i + 1..points.len() {
                // (x - xj) for j<i
                let (xj, _yj) = &points[j];
                self.mul_mut(
                    &mut term,
                    &Polynomial::from_coeffs(vec![self.coeff_ring.neg(xj), self.coeff_ring.one()]),
                );
            }
            for j in 0..i {
                // (xj - x) for i<j
                let (xj, _yj) = &points[j];
                self.mul_mut(
                    &mut term,
                    &Polynomial::from_coeffs(vec![
                        xj.clone(),
                        self.coeff_ring.neg(&self.coeff_ring.one()),
                    ]),
                );
            }

            for j in 0..points.len() {
                for k in j + 1..points.len() {
                    if i != j && i != k {
                        let (xj, _yj) = &points[j];
                        let (xk, _yk) = &points[k];
                        self.mul_mut(
                            &mut term,
                            &Polynomial::constant(
                                self.coeff_ring.add(&self.coeff_ring.neg(xk), xj),
                            ),
                        );
                    }
                }
            }
            self.add_mut(&mut numerator, &term);
        }

        let mut denominator = self.one();
        for i in 0..points.len() {
            let (xi, _yi) = &points[i];
            for j in i + 1..points.len() {
                let (xj, _yj) = &points[j];
                self.mul_mut(
                    &mut denominator,
                    &Polynomial::constant(self.coeff_ring.add(&self.coeff_ring.neg(xj), xi)),
                );
            }
        }

        match self.div(&numerator, &denominator) {
            Ok(interp_poly) => Some(interp_poly),
            Err(RingDivisionError::NotDivisible) => {
                //no such polynomial exists
                None
            }
            Err(RingDivisionError::DivideByZero) => {
                panic!("are the input points distinct?");
            }
        }
    }
}

impl<RS: BezoutDomainStructure> PolynomialStructure<RS> {
    pub fn interpolate_by_linear_system(
        &self,
        points: &Vec<(RS::Set, RS::Set)>,
    ) -> Option<Polynomial<RS::Set>> {
        /*
        e.g. finding a degree 2 polynomial f(x)=a+bx+cx^2 such that
        f(1)=3
        f(2)=1
        f(3)=-2
        is the same as solving the linear system for a, b, c
        / 1 1 1 \ / a \   / 3  \
        | 1 2 4 | | b | = | 1  |
        \ 1 3 9 / \ c /   \ -2 /
        */

        let matrix_structure = MatrixStructure::new(self.coeff_ring());

        let n = points.len();
        let mut mat = matrix_structure.zero(n, n);
        for r in 0..n {
            let (x, _y) = &points[r];
            let mut x_pow = self.coeff_ring().one();
            for c in 0..n {
                *mat.at_mut(r, c).unwrap() = x_pow.clone();
                self.coeff_ring().mul_mut(&mut x_pow, x);
            }
        }

        let mut output_vec = matrix_structure.zero(n, 1);
        for r in 0..n {
            let (_x, y) = &points[r];
            *output_vec.at_mut(r, 0).unwrap() = y.clone();
        }

        match matrix_structure.col_solve(&mat, output_vec) {
            Some(coeff_vec) => Some(Polynomial::from_coeffs(
                (0..n)
                    .map(|i| coeff_vec.at(i, 0).unwrap().clone())
                    .collect(),
            )),
            None => None,
        }
    }
}

impl<FS: FieldOfFractionsStructure> PolynomialStructure<FS>
where
    FS::RS: GreatestCommonDivisorStructure,
{
    pub fn factor_primitive_fof(
        &self,
        p: &Polynomial<FS::Set>,
    ) -> (FS::Set, Polynomial<<FS::RS as Structure>::Set>) {
        let div = self.coeff_ring.base_ring_structure().lcm_list(
            p.coeffs()
                .into_iter()
                .map(|c| self.coeff_ring.denominator(&c))
                .collect(),
        );

        let (mul, prim) = PolynomialStructure::new(self.coeff_ring.base_ring_structure())
            .factor_primitive(p.apply_map(|c| {
                self.coeff_ring
                    .as_base_ring(
                        self.coeff_ring
                            .mul(&self.coeff_ring.from_base_ring(div.clone()), c),
                    )
                    .unwrap()
            }))
            .unwrap();

        (
            self.coeff_ring
                .div(
                    &self.coeff_ring.from_base_ring(mul),
                    &self.coeff_ring.from_base_ring(div),
                )
                .unwrap(),
            prim,
        )
    }

    pub fn primitive_part_fof(
        &self,
        p: &Polynomial<FS::Set>,
    ) -> Polynomial<<FS::RS as Structure>::Set> {
        self.factor_primitive_fof(p).1
    }
}

impl<R: MetaType> MetaType for Polynomial<R>
where
    R::Structure: SemiRingStructure,
{
    type Structure = PolynomialStructure<R::Structure>;

    fn structure() -> Rc<Self::Structure> {
        PolynomialStructure::new(R::structure()).into()
    }
}

impl<R: MetaType> Polynomial<R>
where
    R::Structure: SemiRingStructure<Set = R>,
{
    fn reduce(self) -> Self {
        Self::structure().reduce_poly(self)
    }
}

impl<R: MetaType> Display for Polynomial<R>
where
    R::Structure: SemiRingStructure + ToStringStructure,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", Self::structure().to_string(self))
    }
}

impl<R: MetaType> PartialEq for Polynomial<R>
where
    R::Structure: SemiRingStructure,
{
    fn eq(&self, other: &Self) -> bool {
        Self::structure().equal(self, other)
    }
}

impl<R: MetaType> Eq for Polynomial<R> where R::Structure: SemiRingStructure {}

impl<R: MetaType> Polynomial<R>
where
    R::Structure: SemiRingStructure<Set = R>,
{
    pub fn var() -> Self {
        Self::structure().var()
    }

    pub fn var_pow(n: usize) -> Self {
        Self::structure().var_pow(n)
    }

    pub fn coeff(&self, i: usize) -> R {
        Self::structure().coeff(self, i).clone()
    }

    pub fn leading_coeff(&self) -> Option<R> {
        Self::structure().leading_coeff(self).cloned()
    }

    pub fn evaluate(&self, x: &R) -> R {
        Self::structure().evaluate(self, x)
    }

    //find p(q(x))
    pub fn compose(p: &Self, q: &Self) -> Self {
        Self::structure().compose(p, q)
    }

    pub fn num_coeffs(&self) -> usize {
        Self::structure().num_coeffs(self)
    }

    //if n = deg(p)
    //return x^n * p(1/x)
    pub fn reversed(&self) -> Self {
        Self::structure().reversed(self)
    }

    pub fn degree(&self) -> Option<usize> {
        Self::structure().degree(self)
    }

    pub fn as_constant(&self) -> Option<R> {
        Self::structure().as_constant(self)
    }

    pub fn is_monic(&self) -> bool {
        Self::structure().is_monic(self)
    }

    pub fn derivative(self) -> Self {
        Self::structure().derivative(self)
    }
}

impl<R: MetaType> Polynomial<R>
where
    R::Structure: IntegralDomainStructure,
{
    pub fn try_quorem(a: &Self, b: &Self) -> Result<(Self, Self), RingDivisionError> {
        Self::structure().try_quorem(a, b)
    }

    pub fn pseudorem(a: &Self, b: &Self) -> Option<Result<Polynomial<R>, &'static str>> {
        Self::structure().pseudorem(a.clone(), b)
    }

    pub fn pseudo_remainder_subresultant_sequence(a: Self, b: Self) -> (Vec<Self>, Vec<R>) {
        Self::structure().pseudo_remainder_subresultant_sequence(a, b)
    }

    pub fn subresultant_gcd(a: &Self, b: &Self) -> Self {
        Self::structure().subresultant_gcd(a.clone(), b.clone())
    }

    pub fn resultant(a: &Self, b: &Self) -> R {
        Self::structure().resultant(a.clone(), b.clone())
    }

    pub fn is_squarefree(&self) -> bool {
        Self::structure().is_squarefree(self)
    }

    pub fn discriminant(self) -> Result<R, &'static str> {
        Self::structure().discriminant(self)
    }

    pub fn interpolate_by_lagrange_basis(points: &Vec<(R, R)>) -> Option<Self> {
        Self::structure().interpolate_by_lagrange_basis(points)
    }
}

impl<R: MetaType> Polynomial<R>
where
    R::Structure: BezoutDomainStructure,
{
    pub fn interpolate_by_linear_system(points: &Vec<(R, R)>) -> Option<Self> {
        Self::structure().interpolate_by_linear_system(points)
    }
}

impl<R: MetaType> Polynomial<R>
where
    R::Structure: GreatestCommonDivisorStructure,
{
    pub fn factor_primitive(self) -> Option<(R, Polynomial<R>)> {
        Self::structure().factor_primitive(self)
    }

    pub fn primitive_part(self) -> Option<Polynomial<R>> {
        Self::structure().primitive_part(self)
    }

    pub fn gcd_by_primitive_subresultant(a: Polynomial<R>, b: Polynomial<R>) -> Polynomial<R> {
        Self::structure().gcd_by_primitive_subresultant(a, b)
    }
}

impl<R: MetaType> Polynomial<R>
where
    R::Structure: GreatestCommonDivisorStructure + CharZeroStructure,
{
    pub fn primitive_squarefree_part(&self) -> Self {
        Self::structure().primitive_squarefree_part(self.clone())
    }
}

impl<F: MetaType> Polynomial<F>
where
    F::Structure: FieldOfFractionsStructure,
    <F::Structure as FieldOfFractionsStructure>::RS: GreatestCommonDivisorStructure,
{
    pub fn factor_primitive_fof(
        &self,
    ) -> (
        F,
        Polynomial<<<F::Structure as FieldOfFractionsStructure>::RS as Structure>::Set>,
    ) {
        Self::structure().factor_primitive_fof(self)
    }

    pub fn primitive_part_fof(
        &self,
    ) -> Polynomial<<<F::Structure as FieldOfFractionsStructure>::RS as Structure>::Set> {
        Self::structure().primitive_part_fof(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::number::finite_fields::quaternary_field::*;
    use crate::structure::elements::IntoErgonomic;

    use super::super::super::structure::structure::*;

    use algebraeon_nzq::integer::*;
    use algebraeon_nzq::rational::*;

    use super::*;

    #[test]
    fn test_constant_var_pow() {
        let ring = Polynomial::<Integer>::structure();
        let p = ring.constant_var_pow(Integer::from(2), 7);
        let q = ring.mul(&ring.from_int(Integer::from(2)), &ring.var_pow(7));
        assert_eq!(p, q);
    }

    #[test]
    fn test_display_poly_over_display_cannonical_ring() {
        let f = Polynomial {
            coeffs: vec![
                Integer::from(-2),
                Integer::from(1),
                Integer::from(2),
                Integer::from(4),
            ],
        };

        //test that this compliles
        println!("{}", f);
        println!("{}", f.into_ergonomic());
        //    Integer : Display
        // => CannonicalRS<Integer> : DisplayableRSStructure
        // => PolynomialStructure<CannonicalRS<Integer>> : DisplayableRSStructure
        // => CannonicalRS<Polynomial<Integer>> : DisplayableRSStructure
        // => Polynomial<Integer> : Display AND RSElement<CannonicalRS<Polynomial<Integer>>> : Display
    }

    #[test]
    fn invariant_reduction() {
        let unreduced = Polynomial::<Integer>::from_coeffs(vec![
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
        ]);
        let reduced = Polynomial::from_coeffs(vec![Integer::from(0), Integer::from(1)]);
        assert_eq!(unreduced, reduced);

        let unreduced = Polynomial::from_coeffs(vec![
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
        ]);
        let reduced = Polynomial::<Integer>::zero();
        assert_eq!(unreduced, reduced);
    }

    #[test]
    fn divisibility_over_integers() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        match Polynomial::div(a.ref_set(), b.ref_set()) {
            Ok(c) => {
                println!("{:?} {:?} {:?}", a, b, c);
                assert_eq!(a, b * c.into_ergonomic())
            }
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) + 1;
        match Polynomial::div(a.ref_set(), b.ref_set()) {
            Ok(_c) => panic!(),
            Err(RingDivisionError::NotDivisible) => {}
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        match Polynomial::div(a.ref_set(), b.ref_set()) {
            Ok(_c) => panic!(),
            Err(RingDivisionError::NotDivisible) => {}
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        let b = 0 * x;
        match Polynomial::div(a.ref_set(), b.ref_set()) {
            Ok(_c) => panic!(),
            Err(RingDivisionError::DivideByZero) => {}
            Err(_) => panic!(),
        }

        let a = 0 * x;
        let b = (x - x) + 5;
        match Polynomial::div(a.ref_set(), b.ref_set()) {
            Ok(c) => {
                assert_eq!(c, Polynomial::zero())
            }
            Err(RingDivisionError::DivideByZero) => panic!(),
            Err(_) => panic!(),
        }

        let a = 3087 * x - 8805 * x.pow(2) + 607 * x.pow(3) + x.pow(4);
        let b = (x - x) + 1;
        match Polynomial::div(a.ref_set(), b.ref_set()) {
            Ok(c) => {
                assert_eq!(c.into_ergonomic(), a)
            }
            Err(RingDivisionError::DivideByZero) => panic!(),
            Err(_) => panic!(),
        }
    }

    #[test]
    fn divisibility_over_f4() {
        let x = &Polynomial::<QuaternaryField>::var().into_ergonomic();

        let a = 1 + x + x.pow(2);
        let b = Polynomial::constant(QuaternaryField::Alpha).into_ergonomic() + x;
        match Polynomial::div(a.ref_set(), b.ref_set()) {
            Ok(c) => {
                println!("{:?} {:?} {:?}", a, b, c);
                assert_eq!(a, b * c.into_ergonomic())
            }
            Err(e) => panic!("{:?}", e),
        }
    }

    #[test]
    fn euclidean() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();

        let a = 1 + x + 3 * x.pow(2) + x.pow(3) + 7 * x.pow(4) + x.pow(5);
        let b = 1 + x + 3 * x.pow(2) + 2 * x.pow(3);
        let (q, r) = Polynomial::quorem(&a.ref_set(), &b.ref_set()).unwrap();
        let (q, r) = (q.into_ergonomic(), r.into_ergonomic());
        println!("{:?} = {:?} * {:?} + {:?}", a, b, q, r);
        assert_eq!(a, &b * &q + &r);

        let a = 3 * x;
        let b = 2 * x;
        let (q, r) = Polynomial::quorem(&a.ref_set(), &b.ref_set()).unwrap();
        let (q, r) = (q.into_ergonomic(), r.into_ergonomic());
        println!("{:?} = {:?} * {:?} + {:?}", a, b, q, r);
        assert_eq!(a, &b * &q + &r);

        let a = 3 * x + 5;
        let b = 2 * x + 1;
        let c = 1 + x + x.pow(2);
        let x = &a * &b;
        let y = &b * &c;

        let g = Polynomial::gcd(x.ref_set(), y.ref_set());

        println!("gcd({:?} , {:?}) = {:?}", x, y, g);
        Polynomial::div(&g, &b.ref_set()).unwrap();
        Polynomial::div(&b.ref_set(), &g).unwrap();
    }

    #[test]
    fn test_pseudo_remainder() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();
        {
            let f = (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5)
                .into_verbose();
            let g = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).into_verbose();

            println!("f = {}", f.to_string());
            println!("g = {}", g.to_string());

            let r1 = Polynomial::pseudorem(&f, &g).unwrap().unwrap();
            println!("r1 = {}", r1.to_string());
            assert_eq!(
                r1.clone().into_ergonomic(),
                -15 * x.pow(4) + 3 * x.pow(2) - 9
            );

            let r2 = Polynomial::pseudorem(&g, &r1).unwrap().unwrap();
            println!("r2 = {}", r2.to_string());
            assert_eq!(r2.into_ergonomic(), 15795 * x.pow(2) + 30375 * x - 59535);
        }
        println!();
        {
            let f = (4 * x.pow(3) + 2 * x - 7).into_verbose();
            let g = Polynomial::zero();

            println!("f = {}", f.to_string());
            println!("g = {}", g.to_string());

            if let None = Polynomial::pseudorem(&f, &g) {
            } else {
                assert!(false);
            }
        }
        println!();
        {
            let f = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).into_verbose();
            let g = (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5)
                .into_verbose();

            println!("f = {}", f.to_string());
            println!("g = {}", g.to_string());

            if let Err(_msg) = Polynomial::pseudorem(&f, &g).unwrap() {
            } else {
                assert!(false);
            }
        }
    }

    #[test]
    fn integer_primitive_and_assoc() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();
        let p1 = (-2 - 4 * x.pow(2)).into_verbose();
        let (g, p2) = p1.factor_primitive().unwrap();
        assert_eq!(g, Integer::from(2));
        let (u, p3) = p2.factor_fav_assoc();
        assert_eq!(u.coeffs[0], Integer::from(-1));
        assert_eq!(p3.into_ergonomic(), 1 + 2 * x.pow(2));
    }

    #[test]
    fn test_evaluate() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();
        let f = (1 + x + 3 * x.pow(2) + x.pow(3) + 7 * x.pow(4) + x.pow(5)).into_verbose();
        assert_eq!(f.evaluate(&Integer::from(3)), Integer::from(868));

        let f = Polynomial::zero();
        assert_eq!(f.evaluate(&Integer::from(3)), Integer::from(0));
    }

    #[test]
    fn test_interpolate_by_lagrange_basis() {
        for points in vec![
            vec![
                (Rational::from(-2), Rational::from(-5)),
                (Rational::from(7), Rational::from(4)),
                (Rational::from(-1), Rational::from(-3)),
                (Rational::from(4), Rational::from(1)),
            ],
            vec![(Rational::from(0), Rational::from(0))],
            vec![(Rational::from(0), Rational::from(1))],
            vec![],
            vec![
                (Rational::from(0), Rational::from(0)),
                (Rational::from(1), Rational::from(1)),
                (Rational::from(2), Rational::from(2)),
            ],
        ] {
            let f = Polynomial::interpolate_by_lagrange_basis(&points).unwrap();
            for (inp, out) in &points {
                assert_eq!(&f.evaluate(&inp), out);
            }
        }

        //f(x)=2x
        match Polynomial::interpolate_by_lagrange_basis(&vec![
            (Integer::from(0), Integer::from(0)),
            (Integer::from(1), Integer::from(2)),
        ]) {
            Some(f) => {
                assert_eq!(
                    f,
                    Polynomial::from_coeffs(vec![Integer::from(0), Integer::from(2)])
                )
            }
            None => panic!(),
        }

        //f(x)=1/2x does not have integer coefficients
        match Polynomial::interpolate_by_lagrange_basis(&vec![
            (Integer::from(0), Integer::from(0)),
            (Integer::from(2), Integer::from(1)),
        ]) {
            Some(_f) => panic!(),
            None => {}
        }
    }

    #[test]
    fn test_interpolate_by_linear_system() {
        for points in vec![
            vec![
                (Rational::from(-2), Rational::from(-5)),
                (Rational::from(7), Rational::from(4)),
                (Rational::from(-1), Rational::from(-3)),
                (Rational::from(4), Rational::from(1)),
            ],
            vec![(Rational::from(0), Rational::from(0))],
            vec![(Rational::from(0), Rational::from(1))],
            vec![],
            vec![
                (Rational::from(0), Rational::from(0)),
                (Rational::from(1), Rational::from(1)),
                (Rational::from(2), Rational::from(2)),
            ],
        ] {
            let f = Polynomial::interpolate_by_linear_system(&points).unwrap();
            for (inp, out) in &points {
                assert_eq!(&f.evaluate(&inp), out);
            }
        }

        //f(x)=2x
        match Polynomial::interpolate_by_linear_system(&vec![
            (Integer::from(0), Integer::from(0)),
            (Integer::from(1), Integer::from(2)),
        ]) {
            Some(f) => {
                assert_eq!(
                    f,
                    Polynomial::from_coeffs(vec![Integer::from(0), Integer::from(2)])
                )
            }
            None => panic!(),
        }

        //f(x)=1/2x does not have integer coefficients
        match Polynomial::interpolate_by_linear_system(&vec![
            (Integer::from(0), Integer::from(0)),
            (Integer::from(2), Integer::from(1)),
        ]) {
            Some(_f) => panic!(),
            None => {}
        }
    }

    #[test]
    fn test_derivative() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();
        let f = (2 + 3 * x - x.pow(2) + 7 * x.pow(3)).into_verbose();
        let g = (3 - 2 * x + 21 * x.pow(2)).into_verbose();
        assert_eq!(f.derivative(), g);

        let f = Polynomial::<Integer>::zero();
        let g = Polynomial::<Integer>::zero();
        assert_eq!(f.derivative(), g);

        let f = Polynomial::<Integer>::one();
        let g = Polynomial::<Integer>::zero();
        assert_eq!(f.derivative(), g);
    }

    #[test]
    fn test_monic() {
        assert!(!Polynomial::<Integer>::zero().is_monic());
        assert!(Polynomial::<Integer>::one().is_monic());
        let x = &Polynomial::<Integer>::var().into_ergonomic();
        let f = (2 + 3 * x - x.pow(2) + 7 * x.pow(3) + x.pow(4)).into_verbose();
        let g = (3 - 2 * x + 21 * x.pow(2)).into_verbose();
        assert!(f.is_monic());
        assert!(!g.is_monic());
    }

    #[test]
    fn test_subresultant_gcd() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();

        let f = (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5)
            .into_verbose();
        let g = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).into_verbose();
        assert_eq!(
            Polynomial::subresultant_gcd(&f, &g),
            Polynomial::constant(Integer::from(260708))
        );

        let f = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).into_verbose();
        let g = (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5)
            .into_verbose();
        assert_eq!(
            Polynomial::subresultant_gcd(&f, &g),
            Polynomial::constant(Integer::from(260708))
        );

        let f = ((x + 2).pow(2) * (2 * x - 3).pow(2)).into_verbose();
        let g = ((3 * x - 1) * (2 * x - 3).pow(2)).into_verbose();
        assert_eq!(
            Polynomial::subresultant_gcd(&f, &g).into_ergonomic(),
            7056 - 9408 * x + 3136 * x.pow(2)
        );

        let f = (x.pow(4) + 1).into_verbose();
        let g = (3 * x.pow(2)).into_verbose();
        println!(
            "{:#?}",
            Polynomial::pseudo_remainder_subresultant_sequence(f.clone(), g.clone())
        );
        println!("{:#?}", Polynomial::resultant(&f, &g));
    }

    // #[test]
    // fn test_squarefree_part_by_yuns() {
    //     let x = &Ergonomic::new(Polynomial::<Integer>::var());
    //     let f = ((x + 1).pow(3) * (2 * x + 3).pow(2)).elem();
    //     let g = ((x + 1) * (2 * x + 3)).elem();
    //     assert_eq!(squarefree_part_by_yuns(&f), g);
    // }

    #[test]
    fn test_discriminant() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();
        // Constant -> undefined
        debug_assert!((0 * x.pow(0)).into_verbose().discriminant().is_err());
        debug_assert!((1 * x.pow(0)).into_verbose().discriminant().is_err());
        debug_assert!((2 * x.pow(0)).into_verbose().discriminant().is_err());

        // Linear f(x) = ax+b -> disc(f) = 1
        debug_assert_eq!(
            (3 * x.pow(1) + 2).into_verbose().discriminant().unwrap(),
            Integer::from(1)
        );

        // Quadratic f(x) = ax^2 + bx + c  ->  disc(f) = b^2-4ac
        debug_assert_eq!(
            (x.pow(2) + 1).into_verbose().discriminant().unwrap(),
            Integer::from(-4)
        );
        debug_assert_eq!(
            (3 * x.pow(2) + 1).into_verbose().discriminant().unwrap(),
            Integer::from(-12)
        );
        debug_assert_eq!(
            (x.pow(2) + 3).into_verbose().discriminant().unwrap(),
            Integer::from(-12)
        );
        debug_assert_eq!(
            (3 * x.pow(2) + 3).into_verbose().discriminant().unwrap(),
            Integer::from(-36)
        );
        debug_assert_eq!(
            (x.pow(2) + x.pow(1) + 1)
                .into_verbose()
                .discriminant()
                .unwrap(),
            Integer::from(1 - 4)
        );
        debug_assert_eq!(
            (x.pow(2) + 2 * x.pow(1) + 1)
                .into_verbose()
                .discriminant()
                .unwrap(),
            Integer::from(4 - 4)
        );
        debug_assert_eq!(
            (x.pow(2) + 3 * x.pow(1) + 1)
                .into_verbose()
                .discriminant()
                .unwrap(),
            Integer::from(9 - 4)
        );

        //Cubic f(x) = ax^3 + bx^2 + cx + d  ->  disc(f) = b^2c^2 - 4ac^3 - 4b^3d - 27a^2d^2 + 18abcd
        for (a, b, c, d) in [
            (1, 0, 0, 0),
            (1, 0, 1, 0),
            (1, 0, 0, 1),
            (1, 1, 1, 1),
            (3, 1, 1, 1),
            (3, 3, 3, 3),
            (2, 3, 5, 7),
            (7, 5, 3, 2),
        ] {
            println!(
                "{}",
                (a * x.pow(3) + b * x.pow(2) + c * x.pow(1) + d).into_verbose()
            );
            println!(
                "{:?}",
                (a * x.pow(3) + b * x.pow(2) + c * x.pow(1) + d)
                    .into_verbose()
                    .discriminant()
            );
            debug_assert_eq!(
                (a * x.pow(3) + b * x.pow(2) + c * x.pow(1) + d)
                    .into_verbose()
                    .discriminant()
                    .unwrap(),
                Integer::from(
                    b * b * c * c - 4 * a * c * c * c - 4 * b * b * b * d - 27 * a * a * d * d
                        + 18 * a * b * c * d
                )
            );
        }
    }

    #[test]
    fn test_factor_primitive_fof() {
        for (f, exp) in vec![
            (
                Polynomial::from_coeffs(vec![
                    Rational::from_integers(1, 2),
                    Rational::from_integers(1, 3),
                ]),
                Polynomial::from_coeffs(vec![Integer::from(3), Integer::from(2)]),
            ),
            (
                Polynomial::from_coeffs(vec![
                    Rational::from_integers(4, 1),
                    Rational::from_integers(6, 1),
                ]),
                Polynomial::from_coeffs(vec![Integer::from(2), Integer::from(3)]),
            ),
        ] {
            let (mul, ans) = f.factor_primitive_fof();
            assert!(Polynomial::are_associate(&ans, &exp));
            assert_eq!(
                Polynomial::mul(
                    &ans.apply_map(|c| Rational::from(c)),
                    &Polynomial::constant(mul)
                ),
                f
            );
        }
    }
}
