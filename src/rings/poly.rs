#![allow(dead_code)]

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;

use std::collections::HashMap;
use std::hash::Hash;

use super::matrix::*;
use super::nzq::*;
use super::ring::*;

pub const ZZ_POLY: PolynomialRing<IntegerRing> = PolynomialRing { ring: &ZZ };
pub const QQ_POLY: PolynomialRing<RationalField> = PolynomialRing { ring: &QQ };

#[derive(Debug, Clone)]
pub struct Polynomial<ElemT: Clone> {
    //vec![c0, c1, c2, c3, ..., cn] represents the polynomial c0 + c1*x + c2*x^2 + c3*x^3 + ... + cn * x^n
    //if non-empty, the last item must not be zero
    coeffs: Vec<ElemT>,
}

impl<ElemT: Clone> Polynomial<ElemT> {
    pub fn coeffs(&self) -> Vec<ElemT> {
        self.coeffs.clone()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolynomialRing<'a, R: ComRing> {
    ring: &'a R,
}

impl<'a, R: ComRing> PolynomialRing<'a, R> {
    pub fn new(ring: &'a R) -> Self {
        Self { ring }
    }
}

impl<'a, R: ComRing> ComRing for PolynomialRing<'a, R> {
    type ElemT = Polynomial<R::ElemT>;

    fn to_string(&self, poly: &Self::ElemT) -> String {
        if poly.coeffs.len() == 0 {
            String::from("0")
        } else {
            let mut s = String::new();
            let mut first = true;
            for (k, c) in poly.coeffs.iter().enumerate() {
                if !self.ring.equal(c, &self.ring.zero()) {
                    if self.ring.equal(c, &self.ring.one()) {
                        if k == 0 {
                            s += "1";
                        } else {
                            s += "+"
                        }
                    // } else if self.ring.equal(c, &self.ring.neg(self.ring.one())) {
                    //     if k == 0 {
                    //         s += "-1";
                    //     } else {
                    //         s += "-"
                    //     }
                    } else {
                        if !first {
                            s += "+";
                        }
                        s += "(";
                        s += &self.ring.to_string(c);
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

    fn equal(&self, a: &Self::ElemT, b: &Self::ElemT) -> bool {
        let n = a.coeffs.len();
        if n != b.coeffs.len() {
            false
        } else {
            (0..n).all(|i| self.ring.equal(&a.coeffs[i], &b.coeffs[i]))
        }
    }

    fn zero(&self) -> <Self as ComRing>::ElemT {
        Polynomial { coeffs: vec![] }
    }

    fn one(&self) -> <Self as ComRing>::ElemT {
        Polynomial {
            coeffs: vec![self.ring.one()],
        }
    }

    fn neg_mut(&self, poly: &mut <Self as ComRing>::ElemT) {
        for coeff in &mut poly.coeffs {
            self.ring.neg_mut(coeff);
        }
    }

    fn neg_ref(&self, poly: &<Self as ComRing>::ElemT) -> Self::ElemT {
        self.neg(poly.clone())
    }

    fn neg(&self, mut poly: Self::ElemT) -> Self::ElemT {
        self.neg_mut(&mut poly);
        poly
    }

    fn add_mut(&self, mut poly: &mut Self::ElemT, x: &Self::ElemT) {
        for i in 0..x.coeffs.len() {
            if i < poly.coeffs.len() {
                self.ring.add_mut(&mut poly.coeffs[i], &x.coeffs[i]);
            } else {
                poly.coeffs.push(x.coeffs[i].clone());
            }
        }
        self.reduce(&mut poly);
    }

    fn mul_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Self::ElemT {
        let mut coeffs = Vec::with_capacity(a.coeffs.len() + b.coeffs.len());
        for _k in 0..a.coeffs.len() + b.coeffs.len() {
            coeffs.push(self.ring.zero());
        }
        for i in 0..a.coeffs.len() {
            for j in 0..b.coeffs.len() {
                self.ring.add_mut(
                    &mut coeffs[i + j],
                    &self.ring.mul_refs(&a.coeffs[i], &b.coeffs[j]),
                );
            }
        }
        let mut ans = Self::ElemT { coeffs };
        self.reduce(&mut ans); //TODO: dont have to do this over an integral domain
        ans
    }

    fn mul_mut(&self, poly: &mut Self::ElemT, x: &Self::ElemT) {
        poly.clone_from(&self.mul_refs(poly, x));
    }

    fn div(&self, a: Self::ElemT, b: Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        self.div_rref(a, &b)
    }

    fn div_lref(&self, a: &Self::ElemT, b: Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        self.div_refs(a, &b)
    }

    fn div_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Result<Self::ElemT, RingDivisionError> {
        let q_res = self.div_rref(a.clone(), b);
        match &q_res {
            Ok(q) => debug_assert!(self.equal(&self.mul_refs(q, b), a)),
            Err(_) => {}
        };
        q_res
    }

    fn div_rref(
        &self,
        mut a: Self::ElemT,
        b: &Self::ElemT,
    ) -> Result<Self::ElemT, RingDivisionError> {
        //try to find q such that q*b == a
        // a0 + a1*x + a2*x^2 + ... + am*x^m = (q0 + q1*x + q2*x^2 + ... + qk*x^k) * (b0 + b1*x + b2*x^2 + ... + bn*x^n)
        // 1 + x + x^2 + x^3 + x^4 + x^5 = (?1 + ?x + ?x^2) * (1 + x + x^2 + x^3)      m=6 k=3 n=4

        let m = a.coeffs.len();
        let n = b.coeffs.len();
        if n == 0 {
            Err(RingDivisionError::DivideByZero)
        } else if m == 0 {
            Ok(self.zero())
        } else if m < n {
            Err(RingDivisionError::NotDivisible)
        } else {
            let k = m - n + 1;
            let mut q = Self::ElemT {
                coeffs: (0..k).map(|_i| self.ring.zero()).collect(),
            };
            for i in (0..k).rev() {
                //a[i+n-1] = q[i] * b[n-1]
                match self
                    .ring
                    .div_refs(&self.coeff(&a, i + n - 1), &self.coeff(b, n - 1))
                {
                    Ok(qc) => {
                        //a -= qc*x^i*b
                        self.add_mut(
                            &mut a,
                            &self.neg(self.mul_var_pow(&self.mul_scalar(&b, &qc), i)),
                        );
                        q.coeffs[i] = qc;
                    }
                    Err(RingDivisionError::NotDivisible) => {
                        return Err(RingDivisionError::NotDivisible);
                    }
                    Err(_) => panic!(),
                }
            }
            if !self.equal(&a, &self.zero()) {
                return Err(RingDivisionError::NotDivisible);
            }
            Ok(q)
        }
    }
}

impl<'a, R: ComRing> PolynomialRing<'a, R> {
    pub fn ring(&self) -> &R {
        self.ring
    }

    fn check_invariants(&self, poly: <Self as ComRing>::ElemT) -> Result<(), &'static str> {
        match poly.coeffs.len() {
            0 => {}
            n => {
                if self.ring.equal(&poly.coeffs[n - 1], &self.ring.zero()) {
                    return Err("polynomial coefficients must not end with a zero");
                }
            }
        };
        Ok(())
    }

    fn reduce(&self, poly: &mut <Self as ComRing>::ElemT) {
        loop {
            if poly.coeffs.len() == 0 {
                return;
            } else {
                if self
                    .ring
                    .equal(&poly.coeffs[poly.coeffs.len() - 1], &self.ring.zero())
                {
                    poly.coeffs.pop();
                } else {
                    return;
                }
            }
        }
    }

    pub fn apply_map<S: ComRing>(
        &self,
        new_ring: &S,
        poly: &<Self as ComRing>::ElemT,
        f: impl Fn(&R::ElemT) -> S::ElemT,
    ) -> Polynomial<S::ElemT> {
        PolynomialRing::new(new_ring).from_coeffs(poly.coeffs.iter().map(f).collect())
    }

    pub fn apply_map_with_powers<S: ComRing>(
        &self,
        new_ring: &S,
        poly: &<Self as ComRing>::ElemT,
        f: impl Fn((usize, &R::ElemT)) -> S::ElemT,
    ) -> Polynomial<S::ElemT> {
        PolynomialRing::new(new_ring).from_coeffs(poly.coeffs.iter().enumerate().map(f).collect())
    }

    //find p(q(x))
    pub fn compose(
        &self,
        p: &<Self as ComRing>::ElemT,
        q: &<Self as ComRing>::ElemT,
    ) -> <Self as ComRing>::ElemT {
        let pp_ring = PolynomialRing::new(self);
        pp_ring.evaluate(&self.apply_map(self, p, |c| self.constant(c.clone())), q)
    }

    //if n = deg(p)
    //return x^n * p(1/x)
    pub fn reversed(&self, poly: &<Self as ComRing>::ElemT) -> <Self as ComRing>::ElemT {
        self.from_coeffs(poly.coeffs.clone().into_iter().rev().collect())
    }

    pub fn from_coeffs(&self, coeffs: Vec<R::ElemT>) -> <Self as ComRing>::ElemT {
        let mut p = Polynomial { coeffs };
        self.reduce(&mut p);
        p
    }

    pub fn constant(&self, x: R::ElemT) -> <Self as ComRing>::ElemT {
        if self.ring.equal(&x, &self.ring.zero()) {
            Polynomial { coeffs: vec![] }
        } else {
            Polynomial { coeffs: vec![x] }
        }
    }

    pub fn constant_var_pow(&self, x: R::ElemT, n: usize) -> <Self as ComRing>::ElemT {
        Polynomial {
            coeffs: (0..n + 1)
                .map(|i| if i < n { self.ring.zero() } else { x.clone() })
                .collect(),
        }
    }

    pub fn var(&self) -> <Self as ComRing>::ElemT {
        Polynomial {
            coeffs: vec![self.ring.zero(), self.ring.one()],
        }
    }

    pub fn var_pow(&self, n: usize) -> <Self as ComRing>::ElemT {
        Polynomial {
            coeffs: (0..n + 1)
                .map(|i| {
                    if i < n {
                        self.ring.zero()
                    } else {
                        self.ring.one()
                    }
                })
                .collect(),
        }
    }

    pub fn mul_var_pow(
        &self,
        poly: &<Self as ComRing>::ElemT,
        n: usize,
    ) -> <Self as ComRing>::ElemT {
        let mut coeffs = vec![];
        for _i in 0..n {
            coeffs.push(self.ring.zero());
        }
        for c in &poly.coeffs {
            coeffs.push(c.clone());
        }
        Polynomial { coeffs }
    }

    pub fn mul_scalar(
        &self,
        poly: &<Self as ComRing>::ElemT,
        x: &R::ElemT,
    ) -> <Self as ComRing>::ElemT {
        let mut ans = Polynomial {
            coeffs: poly
                .coeffs
                .iter()
                .map(|c| self.ring.mul_refs(c, x))
                .collect(),
        };
        self.reduce(&mut ans);
        ans
    }

    //getting stuff
    pub fn coeff(&self, poly: &<Self as ComRing>::ElemT, i: usize) -> R::ElemT {
        if i < poly.coeffs.len() {
            poly.coeffs[i].clone()
        } else {
            self.ring.zero()
        }
    }

    //zero -> None
    //const -> 0
    //linear -> 1
    //quadratic -> 2
    //etc.
    pub fn degree(&self, poly: &<Self as ComRing>::ElemT) -> Option<usize> {
        if poly.coeffs.len() == 0 {
            None
        } else {
            Some(poly.coeffs.len() - 1)
        }
    }

    pub fn as_constant(&self, poly: &<Self as ComRing>::ElemT) -> Option<R::ElemT> {
        if poly.coeffs.len() == 0 {
            Some(self.ring.zero())
        } else if poly.coeffs.len() == 1 {
            Some(poly.coeffs[0].clone())
        } else {
            None
        }
    }

    pub fn evaluate(&self, poly: &<Self as ComRing>::ElemT, x: &R::ElemT) -> R::ElemT {
        // f(x) = a + bx + cx^2 + dx^3
        // evaluate as f(x) = a + x(b + x(c + x(d)))
        let mut y = self.ring.zero();
        for c in poly.coeffs.iter().rev() {
            self.ring.mul_mut(&mut y, x);
            self.ring.add_mut(&mut y, c)
        }
        y
    }

    pub fn derivative(&self, mut poly: <Self as ComRing>::ElemT) -> <Self as ComRing>::ElemT {
        if poly.coeffs.len() > 0 {
            for i in 0..poly.coeffs.len() - 1 {
                poly.coeffs[i] = poly.coeffs[i + 1].clone();
                self.ring.mul_mut(
                    &mut poly.coeffs[i],
                    &self.ring.from_int(&Integer::from(i + 1)),
                );
            }
            poly.coeffs.pop();
        }
        poly
    }
}

impl<'a, R: IntegralDomain> PolynomialRing<'a, R> {
    pub fn pseudorem(
        &self,
        a: <Self as ComRing>::ElemT,
        b: <Self as ComRing>::ElemT,
    ) -> Option<Result<<Self as ComRing>::ElemT, &'static str>> {
        self.pseudorem_rref(a, &b)
    }

    pub fn pseudorem_lref(
        &self,
        a: &<Self as ComRing>::ElemT,
        b: <Self as ComRing>::ElemT,
    ) -> Option<Result<<Self as ComRing>::ElemT, &'static str>> {
        self.pseudorem_refs(a, &b)
    }

    pub fn pseudorem_refs(
        &self,
        a: &<Self as ComRing>::ElemT,
        b: &<Self as ComRing>::ElemT,
    ) -> Option<Result<<Self as ComRing>::ElemT, &'static str>> {
        self.pseudorem_rref(a.clone(), b)
    }

    //None if b = 0
    //error if deg(a) < deg(b)
    pub fn pseudorem_rref(
        &self,
        mut a: <Self as ComRing>::ElemT,
        b: &<Self as ComRing>::ElemT,
    ) -> Option<Result<<Self as ComRing>::ElemT, &'static str>> {
        let m = a.coeffs.len();
        let n = b.coeffs.len();

        if n == 0 {
            None
        } else if m < n {
            Some(Err("Should have deg(a) >= deg(b) for pseudo remainder"))
        } else {
            self.mul_mut(
                &mut a,
                &self.constant(
                    self.ring
                        .nat_pow(&self.coeff(b, n - 1), &Natural::from(m - n + 1)),
                ),
            );

            let k = m - n + 1;
            let mut q = Polynomial {
                coeffs: (0..k).map(|_i| self.ring.zero()).collect(),
            };
            for i in (0..k).rev() {
                //a[i+n-1] = q[i] * b[n-1]
                match self
                    .ring
                    .div_rref(self.coeff(&a, i + n - 1), &b.coeffs[n - 1])
                {
                    Ok(qc) => {
                        //a -= qc*x^i*b
                        self.add_mut(
                            &mut a,
                            &self.neg(self.mul_var_pow(&self.mul_scalar(&b, &qc), i)),
                        );
                        q.coeffs[i] = qc;
                    }
                    Err(_) => panic!(),
                }
            }
            Some(Ok(a))
        }
    }

    //efficiently compute the gcd of a and b up to scalar multipication using pseudo remainder subresultant sequence
    //the returned polynomial should the smallest non-zero subresultant polynomial
    //https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Trivial_pseudo-remainder_sequence
    pub fn subresultant_gcd(
        &self,
        mut a: Polynomial<R::ElemT>,
        mut b: Polynomial<R::ElemT>,
    ) -> Polynomial<R::ElemT> {
        match self.degree(&a) {
            None => b,
            Some(mut a_deg) => match self.degree(&b) {
                None => a,
                Some(mut b_deg) => {
                    if a_deg < b_deg {
                        (a, b) = (b, a);
                        (a_deg, b_deg) = (b_deg, a_deg);
                    }
                    let mut beta = {
                        if (a_deg - b_deg) % 2 == 0 {
                            self.ring.from_int(&Integer::from(-1))
                        } else {
                            self.ring.from_int(&Integer::from(1))
                        }
                    };
                    let mut psi = self.ring.from_int(&Integer::from(-1));
                    loop {
                        let d = self.degree(&a).unwrap() - self.degree(&b).unwrap();
                        let gamma = self.coeff(&b, self.degree(&b).unwrap());
                        let r = self
                            .div(
                                self.pseudorem_rref(a, &b).unwrap().unwrap(),
                                self.constant(beta),
                            )
                            .unwrap();
                        (a, b) = (b, r);

                        if self.equal(&b, &self.zero()) {
                            break;
                        }

                        if d == 0 {
                            //can only happen in the first loop
                            debug_assert!(self.ring.equal(&psi, &self.ring.neg(self.ring.one())));
                            psi = self.ring.one();
                        } else {
                            psi = self
                                .ring
                                .div(
                                    self.ring
                                        .nat_pow(&self.ring.neg_ref(&gamma), &Natural::from(d)),
                                    self.ring.nat_pow(&psi, &Natural::from(d - 1)),
                                )
                                .unwrap();
                        }
                        beta = self.ring.mul(
                            self.ring.neg(gamma),
                            self.ring.nat_pow(
                                &psi,
                                &Natural::from(self.degree(&a).unwrap() - self.degree(&b).unwrap()),
                            ),
                        );
                    }
                    a
                }
            },
        }
    }

    pub fn resultant(&self, a: Polynomial<R::ElemT>, b: Polynomial<R::ElemT>) -> R::ElemT {
        let subresultant_gcd = self.subresultant_gcd(a, b);
        match self.as_constant(&subresultant_gcd) {
            Some(res) => res,
            None => self.ring.zero(),
        }
    }
}

impl<'a, R: IntegralDomain> IntegralDomain for PolynomialRing<'a, R> {}

// impl<R: UniqueFactorizationDomain> UniqueFactorizationDomain for PolynomialRing<R> {}

impl<'a, R: GreatestCommonDivisorDomain> PolynomialRing<'a, R> {
    pub fn factor_primitive(
        &self,
        mut poly: Polynomial<R::ElemT>,
    ) -> Option<(R::ElemT, Polynomial<R::ElemT>)> {
        if self.equal(&poly, &self.zero()) {
            None
        } else {
            let g = self.ring.gcd_list(poly.coeffs.iter().collect());
            for i in 0..poly.coeffs.len() {
                poly.coeffs[i] = self.ring.div_refs(&poly.coeffs[i], &g).unwrap()
            }
            Some((g, poly))
        }
    }

    pub fn primitive_part(&self, poly: Polynomial<R::ElemT>) -> Option<Polynomial<R::ElemT>> {
        match self.factor_primitive(poly) {
            Some((unit, prim)) => Some(prim),
            None => None,
        }
    }
}

impl<'a, R: GreatestCommonDivisorDomain + CharacteristicZero> PolynomialRing<'a, R> {
    pub fn primitive_squarefree_part(&self, f: Polynomial<R::ElemT>) -> Polynomial<R::ElemT> {
        if self.equal(&f, &self.zero()) {
            f
        } else {
            let g = self.subresultant_gcd(f.clone(), self.derivative(f.clone()));
            let (_c, g_prim) = self.factor_primitive(g).unwrap();
            let (_c, f_prim) = self.factor_primitive(f).unwrap();
            let f_prim_sqfree = self.div(f_prim, g_prim).unwrap();
            f_prim_sqfree
        }
    }
}

impl<'a, R: FavoriteAssociate + IntegralDomain> FavoriteAssociate for PolynomialRing<'a, R> {
    fn factor_fav_assoc(
        &self,
        mut elem: Polynomial<R::ElemT>,
    ) -> (Polynomial<R::ElemT>, Polynomial<R::ElemT>) {
        if self.equal(&elem, &self.zero()) {
            (self.one(), self.zero())
        } else {
            let (u, _c) = self
                .ring
                .factor_fav_assoc(elem.coeffs[elem.coeffs.len() - 1].clone());
            for i in 0..elem.coeffs.len() {
                elem.coeffs[i] = self.ring.div_refs(&elem.coeffs[i], &u).unwrap()
            }
            (self.constant(u), elem)
        }
    }
}

impl<'a, R: CharacteristicZero> CharacteristicZero for PolynomialRing<'a, R> {}

impl<'a, R: IntegralDomain + FiniteUnits> FiniteUnits for PolynomialRing<'a, R> {
    fn all_units(&self) -> Vec<Polynomial<R::ElemT>> {
        self.ring
            .all_units()
            .into_iter()
            .map(|u| self.constant(u))
            .collect()
    }
}

// pub trait InterpolatablePolynomials: ComRing {
//     fn interpolate(points: &Vec<(Self::ElemT, Self::ElemT)>) -> Option<Polynomial<Self>>;
// }

impl<'a, F: Field> EuclideanDomain for PolynomialRing<'a, F> {
    fn norm(&self, elem: &Self::ElemT) -> Option<Natural> {
        if self.equal(elem, &self.zero()) {
            None
        } else {
            Some(Natural::from(elem.coeffs.len() - 1))
        }
    }

    fn quorem(&self, a: Self::ElemT, b: Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        self.quorem_rref(a, &b)
    }

    fn quorem_lref(&self, a: &Self::ElemT, b: Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        self.quorem_refs(a, &b)
    }

    fn quorem_refs(&self, a: &Self::ElemT, b: &Self::ElemT) -> Option<(Self::ElemT, Self::ElemT)> {
        let res = self.quorem_rref(a.clone(), b);
        match &res {
            Some((q, r)) => debug_assert!(self.equal(&self.add_ref(self.mul_refs(q, b), r), a)),
            None => {}
        };
        res
    }

    fn quorem_rref(
        &self,
        mut a: Self::ElemT,
        b: &Self::ElemT,
    ) -> Option<(Self::ElemT, Self::ElemT)> {
        //try to find q such that q*b == a
        // a0 + a1*x + a2*x^2 + ... + am*x^m = (q0 + q1*x + q2*x^2 + ... + qk*x^k) * (b0 + b1*x + b2*x^2 + ... + bn*x^n)
        // 1 + x + x^2 + x^3 + x^4 + x^5 = (?1 + ?x + ?x^2) * (1 + x + x^2 + x^3)      m=6 k=3 n=4
        let m = a.coeffs.len();
        let n = b.coeffs.len();

        if n == 0 {
            None
        } else if m < n {
            Some((self.zero(), a))
        } else {
            let k = m - n + 1;
            let mut q = Self::ElemT {
                coeffs: (0..k).map(|_i| self.ring.zero()).collect(),
            };
            for i in (0..k).rev() {
                //a[i+n-1] = q[i] * b[n-1]
                match self
                    .ring
                    .div_rref(self.coeff(&a, i + n - 1), &b.coeffs[n - 1])
                {
                    Ok(qc) => {
                        //a -= qc*x^i*b
                        self.add_mut(
                            &mut a,
                            &self.neg(self.mul_var_pow(&self.mul_scalar(b, &qc), i)),
                        );
                        q.coeffs[i] = qc;
                    }
                    Err(_) => panic!(),
                }
            }
            Some((q, a))
        }
    }
}

impl<'a, R: IntegralDomain> PolynomialRing<'a, R> {
    pub fn interpolate_by_lagrange_basis(
        &self,
        points: &Vec<(R::ElemT, R::ElemT)>,
    ) -> Option<Polynomial<R::ElemT>> {
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
            let mut term = self.constant(yi.clone());

            for j in i + 1..points.len() {
                // (x - xj) for j<i
                let (xj, _yj) = &points[j];
                self.mul_mut(
                    &mut term,
                    &self.from_coeffs(vec![self.ring.neg_ref(xj), self.ring.one()]),
                );
            }
            for j in 0..i {
                // (xj - x) for i<j
                let (xj, _yj) = &points[j];
                self.mul_mut(
                    &mut term,
                    &self.from_coeffs(vec![xj.clone(), self.ring.neg(self.ring.one())]),
                );
            }

            for j in 0..points.len() {
                for k in j + 1..points.len() {
                    if i != j && i != k {
                        let (xj, _yj) = &points[j];
                        let (xk, _yk) = &points[k];
                        self.mul_mut(
                            &mut term,
                            &self.constant(self.ring.add_ref(self.ring.neg_ref(xk), xj)),
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
                    &self.constant(self.ring.add_ref(self.ring.neg_ref(xj), xi)),
                );
            }
        }

        match self.div(numerator, denominator) {
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

impl<'a, R: PrincipalIdealDomain> PolynomialRing<'a, R> {
    pub fn interpolate_by_linear_system(
        &self,
        points: &Vec<(R::ElemT, R::ElemT)>,
    ) -> Option<Polynomial<R::ElemT>> {
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
        let n = points.len();
        let mut mat = MatrixStructure::new(self.ring).zero(n, n);
        for r in 0..n {
            let (x, _y) = &points[r];
            let mut x_pow = self.ring.one();
            for c in 0..n {
                *mat.at_mut(r, c).unwrap() = x_pow.clone();
                self.ring.mul_mut(&mut x_pow, x);
            }
        }

        let mut output_vec = MatrixStructure::new(self.ring).zero(n, 1);
        for r in 0..n {
            let (_x, y) = &points[r];
            *output_vec.at_mut(r, 0).unwrap() = y.clone();
        }

        match MatrixStructure::new(self.ring).col_solve(&mat, output_vec) {
            Some(coeff_vec) => Some(
                self.from_coeffs(
                    (0..n)
                        .map(|i| coeff_vec.at(i, 0).unwrap().clone())
                        .collect(),
                ),
            ),
            None => None,
        }
    }
}

impl<'a, R: UniqueFactorizationDomain + GreatestCommonDivisorDomain + FiniteUnits>
    PolynomialRing<'a, R>
{
    fn factor_primitive_linear_part(
        &self,
        mut f: Polynomial<R::ElemT>,
    ) -> (Factored<Polynomial<R::ElemT>>, Polynomial<R::ElemT>) {
        //f should be a primitive polynomial over R
        //linear factor (a+bx) of c0 + c1*x + c2*x^2 + ... + cn*x^n
        //must be such that a divides c0 and b divides cn
        //so just factor c0 and cn and check all divisors

        let mut linear_factors = Factored::new_unchecked(self.one(), vec![]);
        'seek_linear_factor: while self.degree(&f).unwrap() > 0 {
            let c0 = self.coeff(&f, 0);
            if self.ring.equal(&c0, &self.ring.zero()) {
                //linear factor of x
                f = self.div(f, self.var()).unwrap();
                linear_factors = Factored::mul(
                    self,
                    linear_factors,
                    Factored::factored_irreducible_unchecked(self, self.var()),
                );
                continue 'seek_linear_factor;
            } else {
                //look for linear factors of the form (a+bx)
                let c0fs = self.ring.factor(&self.coeff(&f, 0)).unwrap();
                let cnfs = self
                    .ring
                    .factor(&self.coeff(&f, self.degree(&f).unwrap()))
                    .unwrap();
                for a_assoc in self.ring.divisors(&c0fs) {
                    for u in self.ring.all_units() {
                        let a = self.ring.mul_ref(u, &a_assoc);
                        for b in self.ring.divisors(&cnfs) {
                            //a ranges over all divisors of c0
                            //b ranges over all divisors factors of cn up to associates
                            //try the linear factor (a+bx)
                            let lin = self.from_coeffs(vec![a.clone(), b]);
                            match self.div_refs(&f, &lin) {
                                Ok(new_f) => {
                                    f = new_f;
                                    linear_factors = Factored::mul(
                                        self,
                                        linear_factors,
                                        Factored::factored_irreducible_unchecked(self, lin),
                                    );
                                    continue 'seek_linear_factor;
                                }
                                Err(RingDivisionError::NotDivisible) => {}
                                Err(RingDivisionError::DivideByZero) => panic!(),
                            }
                        }
                    }
                }
            }

            break;
        }
        (linear_factors, f)
    }
}

impl<
        'a,
        R: UniqueFactorizationDomain + GreatestCommonDivisorDomain + CharacteristicZero + FiniteUnits,
    > PolynomialRing<'a, R>
{
    fn factorize_primitive_polynomial_by_kroneckers_method(
        &self,
        f: Polynomial<R::ElemT>,
    ) -> Factored<Polynomial<R::ElemT>> {
        debug_assert!(!self.equal(&f, &self.zero()));
        fn partial_factor<
            'a,
            R: UniqueFactorizationDomain
                + GreatestCommonDivisorDomain
                + CharacteristicZero
                + FiniteUnits,
        >(
            poly_ring: &PolynomialRing<'a, R>,
            f: Polynomial<R::ElemT>,
        ) -> Option<(Polynomial<R::ElemT>, Polynomial<R::ElemT>)> {
            /*
            Suppose we want to factor f(x) = 2 + x + x^2 + x^4 + x^5
            Assume it has a proper factor g(x). wlog g(x) has degree <= 2
            g(x) is determined by its value at 3 points, say at x=0, x=1, x=-1
            f(0)=2, f(1)=6, f(-1)=2     if one of these was zero, then we would have found a linear factor
            g(0) divides 2, g(1) divides 6, g(-1) divides 2
            there are finitely many possible values of g(0), g(1) and g(-1) which satisfy these
            infact there are 4*8*4=128 possible triples
            however, only 64 need to be checked as the other half are their negatives
            more abstractly, some possibilities can be avoided because we only care about g up to multiplication by a unit
             */
            let f_deg = poly_ring.degree(&f).unwrap();
            if f_deg == 1 {
                //linear factor is irreducible
                None
            } else {
                let max_factor_degree = f_deg / 2;
                let mut f_points = vec![];
                let mut elem_gen = poly_ring.ring.generate_distinct_elements();
                //take more samples than necessary, then take the subset with the smallest number of divisors
                while f_points.len() < 3 * (max_factor_degree + 1) {
                    //loop terminates because polynomial over integral domain has finitely many roots
                    let x = elem_gen.next().unwrap();
                    let y = poly_ring.evaluate(&f, &x);
                    if !poly_ring.ring.equal(&y, &poly_ring.ring.zero()) {
                        f_points.push((x, poly_ring.ring.factor(&y).unwrap()));
                    }
                }

                //compute all factors of each y value. choose the y with the most divisors to only factor up to units
                f_points.sort_by_cached_key(|(_x, yf)| poly_ring.ring.count_divisors(yf));
                let _ = f_points.split_off(max_factor_degree + 1);
                //possible_g_points is (x, possible_y_values)
                let all_possible_g_points: Vec<(R::ElemT, Vec<R::ElemT>)> = f_points
                    .into_iter()
                    .rev()
                    .enumerate()
                    .map(|(i, (x, yf))| {
                        let mut y_divs = vec![];
                        for d in poly_ring.ring.divisors(&yf) {
                            if i == 0 {
                                //take divisors up to associates for one, because we only care about g up to associates
                                y_divs.push(d);
                            } else {
                                //take _all_ divisors for the rest
                                for u in poly_ring.ring.all_units() {
                                    y_divs.push(poly_ring.ring.mul_ref(u, &d));
                                }
                            }
                        }
                        (x, y_divs)
                    })
                    .collect();

                for possible_g_points in itertools::Itertools::multi_cartesian_product(
                    all_possible_g_points
                        .into_iter()
                        .map(|(x, y_divs)| y_divs.into_iter().map(move |y_div| (x.clone(), y_div))),
                ) {
                    // println!("{:?}", possible_g_points);
                    match poly_ring.interpolate_by_lagrange_basis(&possible_g_points) {
                        Some(g) => {
                            if poly_ring.degree(&g).unwrap() >= 1 {
                                //g is a possible proper divisor of f
                                match poly_ring.div_refs(&f, &g) {
                                    Ok(h) => {
                                        //g really is a proper divisor of f
                                        return Some((g, h));
                                    }
                                    Err(RingDivisionError::NotDivisible) => {}
                                    Err(RingDivisionError::DivideByZero) => panic!(),
                                }
                            }
                        }
                        None => {}
                    }
                }
                //f is irreducible
                None
            }
        }
        full_factor_using_partial_factor(self, f, &partial_factor::<R>)
    }

    pub fn factorize_by_kroneckers_method(
        &self,
        f: &Polynomial<R::ElemT>,
    ) -> Option<Factored<Polynomial<R::ElemT>>> {
        if self.equal(f, &self.zero()) {
            return None;
        }

        // println!();
        // println!("factorize_by_kroneckers_method");
        // println!("f = {:?}", f);
        let (scalar_part, f) = self.factor_primitive(f.clone()).unwrap();
        // println!("scalar_part = {:?}", scalar_part);
        // println!("primitive_part = {:?}", f);
        let factored_scalar_part = self.ring.factor(&scalar_part).unwrap();
        // println!("factored_scalar_part = {:?}", factored_scalar_part);
        let factored_scalar_part_poly = Factored::new_unchecked(
            self.constant(factored_scalar_part.unit().clone()),
            factored_scalar_part
                .factors()
                .iter()
                .map(|(p, k)| (self.constant(p.clone()), k.clone()))
                .collect(),
        );
        let (linear_factors, f) = self.factor_primitive_linear_part(f);

        let f_sqfree = self.primitive_squarefree_part(f.clone());
        let f_sqfree_factors = self.factorize_primitive_polynomial_by_kroneckers_method(f_sqfree);

        let mut f = f;
        let mut f_factors = Factored::new_unchecked(self.one(), vec![]);
        for (sqfree_factor, k) in f_sqfree_factors.factors() {
            debug_assert_eq!(k, &Natural::from(1u8)); //the factorization should be squarefree
            loop {
                match self.div_refs(&f, sqfree_factor) {
                    Ok(new_f) => {
                        f = new_f;
                        f_factors = Factored::mul(
                            self,
                            f_factors,
                            Factored::new_unchecked(
                                self.one(),
                                vec![(sqfree_factor.clone(), Natural::from(1u8))],
                            ),
                        )
                    }
                    Err(RingDivisionError::NotDivisible) => {
                        break;
                    }
                    Err(RingDivisionError::DivideByZero) => panic!(),
                }
            }
        }
        debug_assert_eq!(self.degree(&f).unwrap(), 0);
        debug_assert!(self.is_unit(f.clone()));
        f_factors = Factored::mul(self, f_factors, Factored::factored_unit_unchecked(self, f));

        return Some(Factored::mul(
            self,
            Factored::mul(self, linear_factors, factored_scalar_part_poly),
            f_factors,
        ));
    }
}

impl<'a, R: Field + FiniteUnits> PolynomialRing<'a, R> {
    pub fn factorize_by_trying_all_factors(
        &self,
        f: Polynomial<R::ElemT>,
    ) -> Factored<Polynomial<R::ElemT>> {
        debug_assert!(!self.equal(&f, &self.zero()));
        fn partial_factor<'a, R: Field + FiniteUnits>(
            poly_ring: &PolynomialRing<'a, R>,
            f: Polynomial<R::ElemT>,
        ) -> Option<(Polynomial<R::ElemT>, Polynomial<R::ElemT>)> {
            let f_deg = poly_ring.degree(&f).unwrap();
            let max_factor_degree = f_deg / 2;
            for d in 0..max_factor_degree {
                for mut coeffs in
                    itertools::Itertools::multi_cartesian_product((0..d + 1).into_iter().map(|d| {
                        let mut all_elems = vec![poly_ring.ring().zero()];
                        all_elems.append(&mut poly_ring.ring.all_units());
                        all_elems
                    }))
                {
                    coeffs.push(poly_ring.ring().one());
                    let g = poly_ring.from_coeffs(coeffs);
                    // println!("{}", self.to_string(&g));
                    match poly_ring.div_refs(&f, &g) {
                        Ok(h) => {
                            return Some((g, h));
                        }
                        Err(RingDivisionError::NotDivisible) => {}
                        Err(RingDivisionError::DivideByZero) => panic!(),
                    }
                }
            }
            None
        }
        full_factor_using_partial_factor(self, f, &partial_factor::<R>)
    }
}

impl<'a, F: FieldOfFractions> PolynomialRing<'a, F>
where
    F::R: EuclideanDomain + FavoriteAssociate,
{
    pub fn factor_primitive_fof(
        &self,
        poly: &Polynomial<F::ElemT>,
    ) -> (
        F::ElemT,
        Polynomial<<<F as FieldOfFractions>::R as ComRing>::ElemT>,
    ) {
        let div = self.ring.base_ring().lcm_list(
            poly.coeffs
                .iter()
                .map(|c| self.ring.denominator(&c))
                .collect(),
        );
        let (mul, prim) = PolynomialRing::new(self.ring.base_ring())
            .factor_primitive(self.apply_map(self.ring.base_ring(), poly, |c| {
                self.ring
                    .as_base_ring(self.ring.mul_ref(self.ring.from_base_ring(div.clone()), c))
                    .unwrap()
            }))
            .unwrap();
        (
            self.ring
                .div(self.ring.from_base_ring(mul), self.ring.from_base_ring(div))
                .unwrap(),
            prim,
        )
    }
}

// fn subsylvester_matrix<R: ComRing>(
//     f_deg: usize,
//     g_deg: usize,
//     f: &Polynomial<R>,
//     g: &Polynomial<R>,
//     k: usize,
// ) -> Matrix<R::ElemT> {
//     match f.degree() {
//         Some(d) => assert!(d <= f_deg),
//         None => {}
//     }
//     match g.degree() {
//         Some(d) => assert!(d <= g_deg),
//         None => {}
//     }
//     assert!(k <= f_deg);
//     assert!(k <= g_deg);

//     let mut smat = Matrix::zero(f_deg + g_deg - k, f_deg + g_deg - 2 * k);
//     for j in 0..g_deg - k {
//         for i in 0..f_deg + 1 {
//             *smat.at_mut(j + i, j).unwrap() = f.coeff(f_deg - i);
//         }
//     }
//     for i in 0..f_deg - k {
//         for j in 0..g_deg + 1 {
//             *smat.at_mut(j + i, g_deg - k + i).unwrap() = g.coeff(g_deg - j);
//         }
//     }
//     smat
// }

// pub fn sylvester_matrix<R: ComRing>(
//     f_deg: usize,
//     g_deg: usize,
//     f: &Polynomial<R>,
//     g: &Polynomial<R>,
// ) -> Matrix<R> {
//     subsylvester_matrix(f_deg, g_deg, f, g, 0)
// }

// pub fn resultant_naive<R: ComRing>(
//     f_deg: usize,
//     g_deg: usize,
//     f: &Polynomial<R>,
//     g: &Polynomial<R>,
// ) -> R {
//     sylvester_matrix(f_deg, g_deg, f, g).det_naive().unwrap()
// }

// //determinant of this is the subresultant polynomial
// pub fn subresultant_matrix<R: ComRing>(
//     f_deg: usize,
//     g_deg: usize,
//     f: &Polynomial<R>,
//     g: &Polynomial<R>,
//     k: usize,
// ) -> Matrix<PolynomialRing<R>> {
//     assert!(k <= f_deg);
//     assert!(k <= g_deg);
//     let tmat = subsylvester_matrix(f_deg, g_deg, f, g, k).apply_map(|a| Polynomial::from(a));
//     let mut vmat = Matrix::<Polynomial<R>>::zero(f_deg + g_deg - 2 * k, f_deg + g_deg - k);
//     for i in 0..f_deg + g_deg - 2 * k - 1 {
//         *vmat.at_mut(i, i).unwrap() = Polynomial::one();
//     }
//     for i in 0..k + 1 {
//         *vmat
//             .at_mut(f_deg + g_deg - 2 * k - 1, f_deg + g_deg - 2 * k - 1 + i)
//             .unwrap() = Polynomial::var_pow(k - i);
//     }
//     let prod_mat = Matrix::mul_refs(&vmat, &tmat).unwrap();
//     prod_mat
// }

// //bad way to compute subresultant polynomials
// pub fn subresultant_naive<R: ComRing>(
//     f_deg: usize,
//     g_deg: usize,
//     f: &Polynomial<R>,
//     g: &Polynomial<R>,
//     k: usize,
// ) -> Polynomial<R> {
//     subresultant_matrix(f_deg, g_deg, f, g, k)
//         .det_naive()
//         .unwrap()
// }

fn hensel_quadratic_lift<const IS_FIELD: bool, R: EuclideanDomain + UniqueFactorizationDomain>(
    ring: R,
    old_ring: EuclideanQuotient<IS_FIELD, R>,
    f: Polynomial<R::ElemT>,
    g: Polynomial<R::ElemT>,
    h: Polynomial<R::ElemT>,
    s: Polynomial<R::ElemT>,
    t: Polynomial<R::ElemT>,
) -> (
    Polynomial<R::ElemT>,
    Polynomial<R::ElemT>,
    Polynomial<R::ElemT>,
    Polynomial<R::ElemT>,
) {
    /*
       Lifts a factorization f=gh (mod m) and Bezout coefficients (s, t) for (g, h) in Z/m, (g, h)
       being coprime modulo m, to a factorization f=g*h* (mod m^2) and Bezout coefficients (s*,t*).

       Input: polynomials f, g, h, s, t in Z[x] and a natural number m such that
              f = gh (mod m) and sg + th = 1 (mod m). We also assume that lc(f)
              is invertible modulo m, h is monic and deg(s) < deg(h) and
              deg(t) < deg(g).

       Output: a list [g*, h*, s*, t*] of polynomials  in Z[x] such that
               f = g*h* (mod m^2) and s*h* + t*h* = 1 (mod m^2). We also have
               g* = g (mod m), h* = h (mod m), s* = s (mod m) and
               t* = t (mod m), with g*, h*, s*, t* also satisfying the
               inequalities on the degrees.
    */
    let old_poly_ring = PolynomialRing::new(&old_ring);

    let new_ring = EuclideanQuotient::new_ring(
        ring.clone(),
        ring.mul_refs(&old_ring.get_n(), &old_ring.get_n()),
    );
    let new_poly_ring = PolynomialRing::new(&new_ring);

    // f = g*h
    debug_assert!(old_poly_ring.equal(&f, &old_poly_ring.mul_refs(&g, &h)));

    // g*s + h*t = 1
    debug_assert!(old_poly_ring.equal(
        &old_poly_ring.one(),
        &old_poly_ring.add(
            old_poly_ring.mul_refs(&g, &s),
            old_poly_ring.mul_refs(&h, &t)
        )
    ));

    // // Lifting the factorization.
    // let e = new_poly_ring.add(f, new_poly_ring.neg(new_poly_ring.mul(g, h))); //(f - g*h)
    // let q, r = (s*e).quo_rem(h)  // Division of se by h in Z/m^2.
    // g_, h_ = g + t*e + q*g, h + r
    
    // //Lifting the Bezout coefficients.
    // b = s*g_ + t*h_ - 1
    // c, d = (s*b).quo_rem(h_)  // Division of sb by h_ in Z/m^2.
    // s_, t_ = s - d, t - t*b - c*g_
    
    // // Return the polynomials as embedded in Z[x].
    // return [g_.change_ring(ZZ), h_.change_ring(ZZ), 
    //         s_.change_ring(ZZ), t_.change_ring(ZZ)]

    todo!();
}

impl PartialEq for Polynomial<Integer> {
    fn eq(&self, other: &Self) -> bool {
        self.coeffs == other.coeffs
    }
}

impl Eq for Polynomial<Integer> {}

impl Hash for Polynomial<Integer> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.coeffs.hash(state);
    }
}

#[cfg(test)]
mod tests {
    use core::panic;

    use super::super::ergonomic::*;
    use super::*;
    use malachite_nz::integer::Integer;
    use malachite_q::Rational;

    #[test]
    fn invariant_reduction() {
        let mut unreduced = Polynomial {
            coeffs: vec![
                Integer::from(0),
                Integer::from(1),
                Integer::from(0),
                Integer::from(0),
            ],
        };
        let reduced = Polynomial {
            coeffs: vec![Integer::from(0), Integer::from(1)],
        };
        ZZ_POLY.reduce(&mut unreduced);
        assert!(ZZ_POLY.equal(&unreduced, &reduced));

        let mut unreduced = Polynomial {
            coeffs: vec![
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
            ],
        };
        let reduced = Polynomial::<Integer> { coeffs: vec![] };
        ZZ_POLY.reduce(&mut unreduced);
        assert!(ZZ_POLY.equal(&unreduced, &reduced));
    }

    #[test]
    fn divisibility() {
        let x = &Ergonomic::new(&ZZ_POLY, ZZ_POLY.var());

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        match ZZ_POLY.div(a.elem(), b.elem()) {
            Ok(c) => {
                println!("{:?} {:?} {:?}", a, b, c);
                assert_eq!(a, b * Ergonomic::new(&ZZ_POLY, c))
            }
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) + 1;
        match ZZ_POLY.div(a.elem(), b.elem()) {
            Ok(_c) => panic!(),
            Err(RingDivisionError::NotDivisible) => {}
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        match ZZ_POLY.div(a.elem(), b.elem()) {
            Ok(_c) => panic!(),
            Err(RingDivisionError::NotDivisible) => {}
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        let b = 0 * x;
        match ZZ_POLY.div(a.elem(), b.elem()) {
            Ok(_c) => panic!(),
            Err(RingDivisionError::DivideByZero) => {}
            Err(_) => panic!(),
        }

        let a = 0 * x;
        let b = (x - x) + 5;
        match ZZ_POLY.div(a.elem(), b.elem()) {
            Ok(c) => {
                assert!(ZZ_POLY.equal(&c, &ZZ_POLY.zero()))
            }
            Err(RingDivisionError::DivideByZero) => panic!(),
            Err(_) => panic!(),
        }

        let a = 3087 * x - 8805 * x.pow(2) + 607 * x.pow(3) + x.pow(4);
        let b = (x - x) + 1;
        match ZZ_POLY.div(a.elem(), b.elem()) {
            Ok(c) => {
                assert!(ZZ_POLY.equal(&c, &a.elem()))
            }
            Err(RingDivisionError::DivideByZero) => panic!(),
            Err(_) => panic!(),
        }
    }

    #[test]
    fn euclidean() {
        let x = &Ergonomic::new(&QQ_POLY, QQ_POLY.var());

        let a = 1 + x + 3 * x.pow(2) + x.pow(3) + 7 * x.pow(4) + x.pow(5);
        let b = 1 + x + 3 * x.pow(2) + 2 * x.pow(3);
        let (q, r) = QQ_POLY.quorem_refs(&a.elem(), &b.elem()).unwrap();
        let (q, r) = (Ergonomic::new(&QQ_POLY, q), Ergonomic::new(&QQ_POLY, r));
        println!("{:?} = {:?} * {:?} + {:?}", a, b, q, r);
        assert_eq!(a, &b * &q + &r);

        let a = 3 * x;
        let b = 2 * x;
        let (q, r) = QQ_POLY.quorem_refs(&a.elem(), &b.elem()).unwrap();
        let (q, r) = (Ergonomic::new(&QQ_POLY, q), Ergonomic::new(&QQ_POLY, r));
        println!("{:?} = {:?} * {:?} + {:?}", a, b, q, r);
        assert_eq!(a, &b * &q + &r);

        let a = 3 * x + 5;
        let b = 2 * x + 1;
        let c = 1 + x + x.pow(2);
        let x = &a * &b;
        let y = &b * &c;
        let g = QQ_POLY.gcd(x.elem(), y.elem());

        println!("gcd({:?} , {:?}) = {:?}", x, y, g);
        QQ_POLY.div_refs(&g, &b.elem()).unwrap();
        QQ_POLY.div_refs(&b.elem(), &g).unwrap();
    }

    #[test]
    fn test_pseudo_remainder() {
        let x = &Ergonomic::new(&ZZ_POLY, ZZ_POLY.var());
        {
            let f = (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5)
                .elem();
            let g = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).elem();

            println!("f = {}", ZZ_POLY.to_string(&f));
            println!("g = {}", ZZ_POLY.to_string(&g));

            let r1 = ZZ_POLY.pseudorem_refs(&f, &g).unwrap().unwrap();
            println!("r1 = {}", ZZ_POLY.to_string(&r1));
            assert!(ZZ_POLY.equal(&r1, &(-15 * x.pow(4) + 3 * x.pow(2) - 9).elem()));

            let r2 = ZZ_POLY.pseudorem_refs(&g, &r1).unwrap().unwrap();
            println!("r2 = {}", ZZ_POLY.to_string(&r2));
            assert!(ZZ_POLY.equal(&r2, &(15795 * x.pow(2) + 30375 * x - 59535).elem()));
        }
        println!();
        {
            let f = (4 * x.pow(3) + 2 * x - 7).elem();
            let g = ZZ_POLY.zero();

            println!("f = {}", ZZ_POLY.to_string(&f));
            println!("g = {}", ZZ_POLY.to_string(&g));

            if let None = ZZ_POLY.pseudorem_refs(&f, &g) {
            } else {
                assert!(false);
            }
        }
        println!();
        {
            let f = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).elem();
            let g = (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5)
                .elem();

            println!("f = {}", ZZ_POLY.to_string(&f));
            println!("g = {}", ZZ_POLY.to_string(&g));

            if let Err(_msg) = ZZ_POLY.pseudorem_refs(&f, &g).unwrap() {
            } else {
                assert!(false);
            }
        }
    }

    #[test]
    fn integer_primitive_and_assoc() {
        let x = &Ergonomic::new(&ZZ_POLY, ZZ_POLY.var());
        let p1 = (-2 - 4 * x.pow(2)).elem();
        let (g, p2) = ZZ_POLY.factor_primitive(p1).unwrap();
        assert_eq!(g, Integer::from(2));
        let (u, p3) = ZZ_POLY.factor_fav_assoc(p2);
        assert_eq!(u.coeffs[0], Integer::from(-1));
        assert_eq!(Ergonomic::new(&ZZ_POLY, p3), 1 + 2 * x.pow(2));
    }

    #[test]
    fn test_evaluate() {
        let x = &Ergonomic::new(&ZZ_POLY, ZZ_POLY.var());
        let f = (1 + x + 3 * x.pow(2) + x.pow(3) + 7 * x.pow(4) + x.pow(5)).elem();
        assert_eq!(ZZ_POLY.evaluate(&f, &Integer::from(3)), Integer::from(868));

        let f = ZZ_POLY.zero();
        assert_eq!(ZZ_POLY.evaluate(&f, &Integer::from(3)), Integer::from(0));
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
            let f = QQ_POLY.interpolate_by_lagrange_basis(&points).unwrap();
            for (inp, out) in &points {
                assert_eq!(&QQ_POLY.evaluate(&f, &inp), out);
            }
        }

        //f(x)=2x
        match ZZ_POLY.interpolate_by_lagrange_basis(&vec![
            (Integer::from(0), Integer::from(0)),
            (Integer::from(1), Integer::from(2)),
        ]) {
            Some(f) => {
                assert!(ZZ_POLY.equal(
                    &f,
                    &ZZ_POLY.from_coeffs(vec![Integer::from(0), Integer::from(2)])
                ))
            }
            None => panic!(),
        }

        //f(x)=1/2x does not have integer coefficients
        match ZZ_POLY.interpolate_by_lagrange_basis(&vec![
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
            let f = QQ_POLY.interpolate_by_linear_system(&points).unwrap();
            for (inp, out) in &points {
                assert_eq!(&QQ_POLY.evaluate(&f, &inp), out);
            }
        }

        //f(x)=2x
        match ZZ_POLY.interpolate_by_linear_system(&vec![
            (Integer::from(0), Integer::from(0)),
            (Integer::from(1), Integer::from(2)),
        ]) {
            Some(f) => {
                assert!(ZZ_POLY.equal(
                    &f,
                    &ZZ_POLY.from_coeffs(vec![Integer::from(0), Integer::from(2)])
                ))
            }
            None => panic!(),
        }

        //f(x)=1/2x does not have integer coefficients
        match ZZ_POLY.interpolate_by_linear_system(&vec![
            (Integer::from(0), Integer::from(0)),
            (Integer::from(2), Integer::from(1)),
        ]) {
            Some(_f) => panic!(),
            None => {}
        }
    }

    #[test]
    fn test_factor_by_kroneckers_method_over_integers() {
        let x = &Ergonomic::<PolynomialRing<IntegerRing>>::new(&ZZ_POLY, ZZ_POLY.var());

        //primitive cases
        let f = ((1 + x).pow(2)).elem();
        assert!(Factored::equal(
            &ZZ_POLY,
            &ZZ_POLY.factorize_by_kroneckers_method(&f).unwrap(),
            &Factored::new_unchecked(ZZ_POLY.one(), vec![((1 + x).elem(), Natural::from(2u8))])
        ));

        let f = (-1 - 2 * x).elem();
        let fs1 = ZZ_POLY.factorize_by_kroneckers_method(&f).unwrap();
        let fs2 = &Factored::new_unchecked(
            ZZ_POLY.neg(ZZ_POLY.one()),
            vec![((1 + 2 * x).elem(), Natural::from(1u8))],
        );
        println!("fs1={:?} fs2={:?}", fs1, fs2);
        assert!(Factored::equal(&ZZ_POLY, &fs1, &fs2));

        let f = (x.pow(5) + x.pow(4) + x.pow(2) + x + 2).elem();
        assert!(Factored::equal(
            &ZZ_POLY,
            &ZZ_POLY.factorize_by_kroneckers_method(&f).unwrap(),
            &Factored::new_unchecked(
                ZZ_POLY.one(),
                vec![
                    ((1 + x + x.pow(2)).elem(), Natural::from(1u8)),
                    ((2 - x + x.pow(3)).elem(), Natural::from(1u8))
                ]
            )
        ));

        let f = (1 + x + x.pow(2)).pow(2).elem();
        assert!(Factored::equal(
            &ZZ_POLY,
            &ZZ_POLY.factorize_by_kroneckers_method(&f).unwrap(),
            &Factored::new_unchecked(
                ZZ_POLY.one(),
                vec![((1 + x + x.pow(2)).elem(), Natural::from(2u8))]
            )
        ));

        //non-primitive cases
        let f = (2 + 2 * x).elem();
        assert!(Factored::equal(
            &ZZ_POLY,
            &ZZ_POLY.factorize_by_kroneckers_method(&f).unwrap(),
            &Factored::new_unchecked(
                ZZ_POLY.one(),
                vec![
                    (ZZ_POLY.from_int(&Integer::from(2)), Natural::from(1u8)),
                    ((1 + x).elem(), Natural::from(1u8))
                ]
            )
        ));

        let f = (12 * (2 + 3 * x) * (x - 1).pow(2)).elem();
        assert!(Factored::equal(
            &ZZ_POLY,
            &ZZ_POLY.factorize_by_kroneckers_method(&f).unwrap(),
            &Factored::new_unchecked(
                ZZ_POLY.one(),
                vec![
                    (ZZ_POLY.from_int(&Integer::from(2)), Natural::from(2u8)),
                    (ZZ_POLY.from_int(&Integer::from(3)), Natural::from(1u8)),
                    ((2 + 3 * x).elem(), Natural::from(1u8)),
                    ((x - 1).elem(), Natural::from(2u8))
                ]
            )
        ));

        let f = ZZ_POLY.one();
        assert!(Factored::equal(
            &ZZ_POLY,
            &ZZ_POLY.factorize_by_kroneckers_method(&f).unwrap(),
            &Factored::new_unchecked(ZZ_POLY.one(), vec![])
        ));

        let f = ((x.pow(4) + x + 1) * (x.pow(3) + x + 1)).elem();
        assert!(Factored::equal(
            &ZZ_POLY,
            &ZZ_POLY.factorize_by_kroneckers_method(&f).unwrap(),
            &Factored::new_unchecked(
                ZZ_POLY.one(),
                vec![
                    ((x.pow(4) + x + 1).elem(), Natural::from(1u8)),
                    ((x.pow(3) + x + 1).elem(), Natural::from(1u8))
                ]
            )
        ));
    }

    #[test]
    fn test_derivative() {
        let x = &Ergonomic::new(&ZZ_POLY, ZZ_POLY.var());
        let f = (2 + 3 * x - x.pow(2) + 7 * x.pow(3)).elem();
        let g = (3 - 2 * x + 21 * x.pow(2)).elem();
        assert!(ZZ_POLY.equal(&ZZ_POLY.derivative(f), &g));

        let f = ZZ_POLY.zero();
        let g = ZZ_POLY.zero();
        assert!(ZZ_POLY.equal(&ZZ_POLY.derivative(f), &g));

        let f = ZZ_POLY.one();
        let g = ZZ_POLY.zero();
        assert!(ZZ_POLY.equal(&ZZ_POLY.derivative(f), &g));
    }

    // #[test]
    // fn test_sylvester_matrix() {
    //     let f = Polynomial::new(vec![
    //         Integer::from(1),
    //         Integer::from(2),
    //         Integer::from(3),
    //         Integer::from(4),
    //     ]);
    //     let g = Polynomial::new(vec![Integer::from(5), Integer::from(6), Integer::from(7)]);
    //     assert_eq!(
    //         sylvester_matrix(3, 2, &f, &g),
    //         Matrix::from_rows(vec![
    //             vec![
    //                 Integer::from(4),
    //                 Integer::from(0),
    //                 Integer::from(7),
    //                 Integer::from(0),
    //                 Integer::from(0)
    //             ],
    //             vec![
    //                 Integer::from(3),
    //                 Integer::from(4),
    //                 Integer::from(6),
    //                 Integer::from(7),
    //                 Integer::from(0)
    //             ],
    //             vec![
    //                 Integer::from(2),
    //                 Integer::from(3),
    //                 Integer::from(5),
    //                 Integer::from(6),
    //                 Integer::from(7)
    //             ],
    //             vec![
    //                 Integer::from(1),
    //                 Integer::from(2),
    //                 Integer::from(0),
    //                 Integer::from(5),
    //                 Integer::from(6)
    //             ],
    //             vec![
    //                 Integer::from(0),
    //                 Integer::from(1),
    //                 Integer::from(0),
    //                 Integer::from(0),
    //                 Integer::from(5)
    //             ]
    //         ])
    //     );
    //     assert_eq!(
    //         subsylvester_matrix(3, 2, &f, &g, 1),
    //         Matrix::from_rows(vec![
    //             vec![Integer::from(4), Integer::from(7), Integer::from(0)],
    //             vec![Integer::from(3), Integer::from(6), Integer::from(7)],
    //             vec![Integer::from(2), Integer::from(5), Integer::from(6)],
    //             vec![Integer::from(1), Integer::from(0), Integer::from(5)],
    //         ])
    //     );

    //     let f = Polynomial::new(vec![Integer::from(1)]);
    //     let g = Polynomial::new(vec![Integer::from(2)]);
    //     let mat = sylvester_matrix(0, 0, &f, &g);
    //     assert_eq!(mat, Matrix::zero(0, 0));

    //     let f = Polynomial::new(vec![Integer::from(1)]);
    //     let g = Polynomial::new(vec![Integer::from(2), Integer::from(3)]);
    //     let mat = sylvester_matrix(0, 1, &f, &g);
    //     assert_eq!(mat, Matrix::from_rows(vec![vec![Integer::from(1)],]));
    // }

    // #[test]
    // fn test_resultant() {
    //     //naive algorithm (naive det of sylvester matrix)
    //     let f = Polynomial::<ZZ>::new(vec![
    //         Integer::from(1),
    //         Integer::from(2),
    //         Integer::from(3),
    //         Integer::from(4),
    //     ]);
    //     let g = Polynomial::<ZZ>::new(vec![Integer::from(5), Integer::from(6), Integer::from(7)]);
    //     assert_eq!(resultant_naive(3, 2, &f, &g), Integer::from(832));

    //     let f = Polynomial::new(vec![
    //         Integer::from(1),
    //         Integer::from(2),
    //         Integer::from(3),
    //         Integer::from(4),
    //     ]);
    //     let g = Polynomial::<ZZ>::new(vec![Integer::from(5), Integer::from(6), Integer::from(7)]);
    //     assert_eq!(resultant_naive(4, 3, &f, &g), Integer::from(0));

    // //subresultant
    // let f = Polynomial::new(vec![
    //     Integer::from(1),
    //     Integer::from(2),
    //     Integer::from(3),
    //     Integer::from(4),
    // ]);
    // let g = Polynomial::<ZZ>::new(vec![Integer::from(5), Integer::from(6), Integer::from(7)]);
    // let subres = subresultant_naive(3, 2, &f, &g, 1);
    // assert_eq!(
    //     subres,
    //     Polynomial::new(vec![Integer::from(64), Integer::from(-24)])
    // );
    // }

    #[test]
    fn test_subresultant_gcd() {
        let x = &Ergonomic::new(&ZZ_POLY, ZZ_POLY.var());

        let f =
            (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5).elem();
        let g = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).elem();
        assert!(ZZ_POLY.equal(
            &ZZ_POLY.subresultant_gcd(f, g),
            &ZZ_POLY.constant(Integer::from(260708))
        ));

        let f = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).elem();
        let g =
            (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5).elem();
        assert!(ZZ_POLY.equal(
            &ZZ_POLY.subresultant_gcd(f, g),
            &ZZ_POLY.constant(Integer::from(260708))
        ));

        let f = ((x + 2).pow(2) * (2 * x - 3).pow(2)).elem();
        let g = ((3 * x - 1) * (2 * x - 3).pow(2)).elem();
        assert!(ZZ_POLY.equal(
            &ZZ_POLY.subresultant_gcd(f, g),
            &(7056 - 9408 * x + 3136 * x.pow(2)).elem()
        ));
    }

    // #[test]
    // fn test_squarefree_part_by_yuns() {
    //     let x = &Ergonomic::new(Polynomial::<Integer>::var());
    //     let f = ((x + 1).pow(3) * (2 * x + 3).pow(2)).elem();
    //     let g = ((x + 1) * (2 * x + 3)).elem();
    //     assert_eq!(squarefree_part_by_yuns(&f), g);
    // }

    #[test]
    fn test_factor_primitive_fof() {
        for (f, exp) in vec![
            (
                QQ_POLY.from_coeffs(vec![
                    Rational::from_signeds(1, 2),
                    Rational::from_signeds(1, 3),
                ]),
                ZZ_POLY.from_coeffs(vec![Integer::from(3), Integer::from(2)]),
            ),
            (
                QQ_POLY.from_coeffs(vec![
                    Rational::from_signeds(4, 1),
                    Rational::from_signeds(6, 1),
                ]),
                ZZ_POLY.from_coeffs(vec![Integer::from(2), Integer::from(3)]),
            ),
        ] {
            let (mul, ans) = QQ_POLY.factor_primitive_fof(&f);
            assert!(ZZ_POLY.are_associate(&ans, &exp));
            assert!(QQ_POLY.equal(
                &QQ_POLY.mul_scalar(&ZZ_POLY.apply_map(&QQ, &ans, |c| Rational::from(c)), &mul),
                &f
            ));
        }
    }
}
