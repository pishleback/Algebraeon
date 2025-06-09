use super::{super::structure::*, Polynomial};
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use std::borrow::Borrow;

#[derive(Debug, Clone)]
pub struct PolynomialSemiRingStructure<RS: SemiRingSignature, RSB: BorrowedStructure<RS>> {
    coeff_ring_zero: RS::Set, //so that we can return a reference to zero when getting polynomial coefficients out of range
    coeff_ring: RSB,
}

impl<RS: SemiRingSignature, RSB: BorrowedStructure<RS>> PolynomialSemiRingStructure<RS, RSB> {
    fn new(coeff_ring: RSB) -> Self {
        Self {
            coeff_ring_zero: coeff_ring.borrow().zero(),
            coeff_ring,
        }
    }

    pub fn coeff_ring(&self) -> &RS {
        self.coeff_ring.borrow()
    }
}

pub trait SemiRingToPolynomialSemiRingSignature: SemiRingSignature {
    fn polynomial_semiring<'a>(&'a self) -> PolynomialSemiRingStructure<Self, &'a Self> {
        PolynomialSemiRingStructure::new(self)
    }

    fn into_polynomial_semiring(self) -> PolynomialSemiRingStructure<Self, Self> {
        PolynomialSemiRingStructure::new(self)
    }
}

impl<RS: SemiRingSignature> SemiRingToPolynomialSemiRingSignature for RS {}

impl<RS: SemiRingSignature, RSB: BorrowedStructure<RS>> Signature
    for PolynomialSemiRingStructure<RS, RSB>
{
}

impl<RS: SemiRingSignature, RSB: BorrowedStructure<RS>> SetSignature
    for PolynomialSemiRingStructure<RS, RSB>
{
    type Set = Polynomial<RS::Set>;

    fn is_element(&self, _x: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl<RS: SemiRingSignature, RSB: BorrowedStructure<RS>> PartialEq
    for PolynomialSemiRingStructure<RS, RSB>
{
    fn eq(&self, other: &Self) -> bool {
        self.coeff_ring == other.coeff_ring
    }
}

impl<RS: SemiRingSignature, RSB: BorrowedStructure<RS>> Eq
    for PolynomialSemiRingStructure<RS, RSB>
{
}

impl<RS: SemiRingSignature + ToStringSignature, RSB: BorrowedStructure<RS>> ToStringSignature
    for PolynomialSemiRingStructure<RS, RSB>
{
    fn to_string(&self, elem: &Self::Set) -> String {
        if self.num_coeffs(elem) == 0 {
            "0".into()
        } else {
            let mut s = String::new();
            let mut first = true;
            for (k, c) in elem.coeffs.iter().enumerate() {
                if !self.coeff_ring().is_zero(c) {
                    if self.coeff_ring().equal(c, &self.coeff_ring().one()) {
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
                        s += &self.coeff_ring().to_string(c);
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

impl<RS: SemiRingSignature, RSB: BorrowedStructure<RS>> EqSignature
    for PolynomialSemiRingStructure<RS, RSB>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        for i in 0..std::cmp::max(a.coeffs.len(), b.coeffs.len()) {
            if !self.coeff_ring().equal(self.coeff(a, i), self.coeff(b, i)) {
                return false;
            }
        }
        true
    }
}

impl<RS: SemiRingSignature, RSB: BorrowedStructure<RS>> PolynomialSemiRingStructure<RS, RSB> {
    pub fn add_impl<'a, C: Borrow<RS::Set>>(
        &self,
        mut a: &'a Polynomial<C>,
        mut b: &'a Polynomial<C>,
    ) -> Polynomial<RS::Set> {
        if a.coeffs.len() > b.coeffs.len() {
            (a, b) = (b, a);
        }
        let x = a.coeffs.len();
        let y = b.coeffs.len();

        let mut coeffs = vec![];
        let mut i = 0;
        while i < x {
            coeffs.push(
                self.coeff_ring()
                    .add(a.coeffs[i].borrow(), b.coeffs[i].borrow()),
            );
            i += 1;
        }
        while i < y {
            coeffs.push(b.coeffs[i].borrow().to_owned());
            i += 1;
        }

        Polynomial::from_coeffs(coeffs)
    }

    pub fn mul_naive(
        &self,
        a: &Polynomial<impl Borrow<RS::Set>>,
        b: &Polynomial<impl Borrow<RS::Set>>,
    ) -> Polynomial<RS::Set> {
        let mut coeffs = Vec::with_capacity(a.coeffs.len() + b.coeffs.len());
        for _k in 0..a.coeffs.len() + b.coeffs.len() {
            coeffs.push(self.coeff_ring().zero());
        }
        for i in 0..a.coeffs.len() {
            for j in 0..b.coeffs.len() {
                self.coeff_ring().add_mut(
                    &mut coeffs[i + j],
                    &self
                        .coeff_ring()
                        .mul(a.coeffs[i].borrow(), b.coeffs[j].borrow()),
                );
            }
        }
        Polynomial::from_coeffs(coeffs)
    }
}

impl<RS: SemiRingSignature, RSB: BorrowedStructure<RS>> AdditiveMonoidSignature
    for PolynomialSemiRingStructure<RS, RSB>
{
    fn zero(&self) -> Self::Set {
        Polynomial { coeffs: vec![] }
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.reduce_poly(self.add_impl(a, b))
    }
}

impl<RS: SemiRingSignature, RSB: BorrowedStructure<RS>> SemiRingSignature
    for PolynomialSemiRingStructure<RS, RSB>
{
    fn one(&self) -> Self::Set {
        Polynomial::from_coeffs(vec![self.coeff_ring().one()])
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.reduce_poly(self.mul_naive(a, b))
    }
}

impl<RS: CharacteristicSignature, RSB: BorrowedStructure<RS>> CharacteristicSignature
    for PolynomialSemiRingStructure<RS, RSB>
{
    fn characteristic(&self) -> Natural {
        self.coeff_ring().characteristic()
    }
}

impl<RS: SemiRingSignature, RSB: BorrowedStructure<RS>> PolynomialSemiRingStructure<RS, RSB> {
    pub fn reduce_poly(&self, mut a: Polynomial<RS::Set>) -> Polynomial<RS::Set> {
        loop {
            if a.coeffs.is_empty() {
                break;
            }
            if self
                .coeff_ring()
                .equal(a.coeffs.last().unwrap(), &self.coeff_ring().zero())
            {
                a.coeffs.pop();
            } else {
                break;
            }
        }
        a
    }

    pub fn var(&self) -> Polynomial<RS::Set> {
        Polynomial::from_coeffs(vec![self.coeff_ring().zero(), self.coeff_ring().one()])
    }

    pub fn var_pow(&self, n: usize) -> Polynomial<RS::Set> {
        Polynomial::from_coeffs(
            (0..=n)
                .map(|i| {
                    if i < n {
                        self.coeff_ring().zero()
                    } else {
                        self.coeff_ring().one()
                    }
                })
                .collect(),
        )
    }

    pub fn constant_var_pow(&self, x: RS::Set, n: usize) -> Polynomial<RS::Set> {
        Polynomial::from_coeffs(
            (0..=n)
                .map(|i| {
                    if i < n {
                        self.coeff_ring().zero()
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
        let mut y = self.coeff_ring().zero();
        for c in p.coeffs().into_iter().rev() {
            self.coeff_ring().mul_mut(&mut y, x);
            self.coeff_ring().add_mut(&mut y, c);
        }
        y
    }

    /// evaluate p(x^k)
    pub fn evaluate_at_var_pow(&self, p: Polynomial<RS::Set>, k: usize) -> Polynomial<RS::Set> {
        if k == 0 {
            Polynomial::constant(self.coeff_ring().sum(p.coeffs()))
        } else {
            let mut coeffs = vec![];
            for (i, c) in p.into_coeffs().into_iter().enumerate() {
                if i != 0 {
                    for _j in 0..(k - 1) {
                        coeffs.push(self.coeff_ring().zero());
                    }
                }
                coeffs.push(c);
            }
            Polynomial::from_coeffs(coeffs)
        }
    }

    //find p(q(x))
    pub fn compose(&self, p: &Polynomial<RS::Set>, q: &Polynomial<RS::Set>) -> Polynomial<RS::Set> {
        PolynomialSemiRingStructure::new(self.clone())
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
            coeffs.push(self.coeff_ring().zero());
        }
        for c in &p.coeffs {
            coeffs.push(c.clone());
        }
        Polynomial { coeffs }
    }

    pub fn eval_var_pow(&self, p: &Polynomial<RS::Set>, n: usize) -> Polynomial<RS::Set> {
        if n == 0 {
            Polynomial::constant(self.coeff_ring().sum(p.coeffs()))
        } else {
            let gap = n - 1;
            let mut coeffs = vec![];
            for (i, coeff) in p.coeffs.iter().enumerate() {
                if i != 0 {
                    for _ in 0..gap {
                        coeffs.push(self.coeff_ring().zero());
                    }
                }
                coeffs.push(coeff.clone());
            }
            Polynomial { coeffs }
        }
    }

    pub fn mul_scalar(&self, p: &Polynomial<RS::Set>, x: &RS::Set) -> Polynomial<RS::Set> {
        self.reduce_poly(Polynomial::from_coeffs(
            p.coeffs
                .iter()
                .map(|c| self.coeff_ring().mul(c, x))
                .collect(),
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
        #[allow(clippy::manual_find)]
        for i in (0..p.coeffs.len()).rev() {
            if !self.coeff_ring().is_zero(&p.coeffs[i]) {
                return Some(i);
            }
        }
        None
    }

    pub fn as_constant(&self, p: &Polynomial<RS::Set>) -> Option<RS::Set> {
        if self.num_coeffs(p) == 0 {
            Some(self.coeff_ring().zero())
        } else if self.num_coeffs(p) == 1 {
            Some(p.coeffs[0].clone())
        } else {
            None
        }
    }

    pub fn is_monic(&self, p: &Polynomial<RS::Set>) -> bool {
        match self.leading_coeff(p) {
            Some(lc) => self.coeff_ring().equal(lc, &self.coeff_ring().one()),
            None => false,
        }
    }

    pub fn derivative(&self, mut p: Polynomial<RS::Set>) -> Polynomial<RS::Set> {
        p = self.reduce_poly(p);
        if self.num_coeffs(&p) > 0 {
            for i in 0..self.num_coeffs(&p) - 1 {
                p.coeffs[i] = p.coeffs[i + 1].clone();
                self.coeff_ring().mul_mut(
                    &mut p.coeffs[i],
                    &self.coeff_ring().from_nat(Natural::from(i + 1)),
                );
            }
            p.coeffs.pop();
        }
        p
    }
}

// impl Display for Polynomial<Natural> {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         write!(f, "{}", Self::structure().to_string(self))
//     }
// }
