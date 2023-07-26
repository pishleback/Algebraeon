#![allow(dead_code)]

use malachite_nz::natural::Natural;

use std::{hash::Hash, marker::PhantomData};

use super::ring::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polynomial<R: ComRing> {
    //vec![c0, c1, c2, c3, ..., cn] represents the polynomial c0 + c1*x + c2*x^2 + c3*x^3 + ... + cn * x^n
    //if non-empty, the last item must not be zero
    coeffs: Vec<R>,
}

impl<R: ComRing + ToString> ToString for Polynomial<R> {
    fn to_string(&self) -> String {
        if self == &Self::zero() {
            String::from("0")
        } else {
            let mut s = String::new();
            let mut first = true;
            for (k, c) in self.coeffs.iter().enumerate() {
                if c != &R::zero() {
                    if !first {
                        s += "+";
                    }
                    s += "(";
                    s += &c.to_string();
                    s += ")";
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

impl<R: ComRing> Polynomial<R> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        match self.coeffs.len() {
            0 => {}
            n => {
                if self.coeffs[n - 1] == R::zero() {
                    return Err("polynomial coefficients must not end with a zero");
                }
            }
        };
        Ok(())
    }

    fn reduce(&mut self) {
        loop {
            if self.coeffs.len() == 0 {
                return;
            } else {
                if self.coeffs[self.coeffs.len() - 1] == R::zero() {
                    self.coeffs.pop();
                } else {
                    return;
                }
            }
        }
    }

    pub fn new(coeffs: Vec<R>) -> Self {
        let mut p = Self { coeffs };
        p.reduce();
        p
    }

    pub fn var() -> Self {
        Self {
            coeffs: vec![R::zero(), R::one()],
        }
    }

    pub fn var_pow(n: usize) -> Self {
        Self {
            coeffs: (0..n + 1)
                .map(|i| if i < n { R::zero() } else { R::one() })
                .collect(),
        }
    }

    pub fn mul_var_pow(&self, n: usize) -> Self {
        let mut coeffs = vec![];
        for _i in 0..n {
            coeffs.push(R::zero());
        }
        for c in &self.coeffs {
            coeffs.push(c.clone());
        }
        Self { coeffs }
    }

    pub fn mul_scalar(&self, x: &R) -> Self {
        let mut ans = Self {
            coeffs: self.coeffs.iter().map(|c| R::mul_refs(c, x)).collect(),
        };
        ans.reduce();
        ans
    }

    pub fn coeff(&self, i: usize) -> R {
        if i < self.coeffs.len() {
            self.coeffs[i].clone()
        } else {
            R::zero()
        }
    }

    //zero -> None
    //const -> 0
    //linear -> 1
    //quadratic -> 2
    //etc.
    pub fn degree(&self) -> Option<usize> {
        if self.coeffs.len() == 0 {
            None
        } else {
            Some(self.coeffs.len() - 1)
        }
    }
}

impl<R: ComRing> From<R> for Polynomial<R> {
    fn from(x: R) -> Self {
        if x == R::zero() {
            Self { coeffs: vec![] }
        } else {
            Self { coeffs: vec![x] }
        }
    }
}

impl<R: ComRing> ComRing for Polynomial<R> {
    fn zero() -> Self {
        Self { coeffs: vec![] }
    }

    fn one() -> Self {
        Self {
            coeffs: vec![R::one()],
        }
    }

    fn neg_mut(&mut self) {
        for coeff in &mut self.coeffs {
            coeff.neg_mut();
        }
    }

    fn add_mut(&mut self, x: &Self) {
        for i in 0..x.coeffs.len() {
            if i < self.coeffs.len() {
                self.coeffs[i].add_mut(&x.coeffs[i]);
            } else {
                self.coeffs.push(x.coeffs[i].clone());
            }
        }
        self.reduce();
    }

    fn mul_refs(a: &Self, b: &Self) -> Self {
        let mut coeffs = Vec::<R>::with_capacity(a.coeffs.len() + b.coeffs.len());
        for _k in 0..a.coeffs.len() + b.coeffs.len() {
            coeffs.push(R::zero());
        }
        for i in 0..a.coeffs.len() {
            for j in 0..b.coeffs.len() {
                coeffs[i + j].add_mut(&R::mul_refs(&a.coeffs[i], &b.coeffs[j]));
            }
        }
        let mut ans = Self { coeffs };
        ans.reduce(); //TODO: dont have to do this over an integral domain
        ans
    }

    fn mul_mut(&mut self, x: &Self) {
        self.clone_from(&Self::mul_refs(self, x));
    }

    fn div(a: Self, b: Self) -> Result<Self, RingOppErr> {
        Self::div_rref(a, &b)
    }

    fn div_lref(a: &Self, b: Self) -> Result<Self, RingOppErr> {
        Self::div_refs(a, &b)
    }

    fn div_refs(a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        let q_res = Self::div_rref(a.clone(), b);
        match &q_res {
            Ok(q) => debug_assert_eq!(&Self::mul_refs(q, b), a),
            Err(_) => {}
        };
        q_res
    }

    fn div_rref(mut a: Self, b: &Self) -> Result<Self, RingOppErr> {
        //try to find q such that q*b == a
        // a0 + a1*x + a2*x^2 + ... + am*x^m = (q0 + q1*x + q2*x^2 + ... + qk*x^k) * (b0 + b1*x + b2*x^2 + ... + bn*x^n)
        // 1 + x + x^2 + x^3 + x^4 + x^5 = (?1 + ?x + ?x^2) * (1 + x + x^2 + x^3)      m=6 k=3 n=4
        let m = a.coeffs.len();
        let n = b.coeffs.len();
        if m == 0 {
            Ok(Self::zero())
        } else if m < n {
            Err(RingOppErr::NotDivisible)
        } else if n == 0 {
            Err(RingOppErr::DivideByZero)
        } else {
            let k = m - n + 1;
            let mut q = Self {
                coeffs: (0..k).map(|_i| R::zero()).collect(),
            };
            for i in (0..k).rev() {
                //a[i+n-1] = q[i] * b[n-1]
                match R::div_refs(&a.coeff(i + n - 1), &b.coeff(n - 1)) {
                    Ok(qc) => {
                        //a -= qc*x^i*b
                        a.add_mut(&b.mul_scalar(&qc).mul_var_pow(i).neg());
                        q.coeffs[i] = qc;
                    }
                    Err(RingOppErr::NotDivisible) => {
                        return Err(RingOppErr::NotDivisible);
                    }
                    Err(_) => panic!(),
                }
            }
            if a != Self::zero() {
                return Err(RingOppErr::NotDivisible);
            }
            Ok(q)
        }
    }
}

impl<R: IntegralDomain> IntegralDomain for Polynomial<R> {}

impl<R: PrincipalIdealDomain> Polynomial<R> {
    fn factor_primitive(mut self) -> Option<(R, Self)> {
        if self == Self::zero() {
            None
        } else {
            let g = R::gcd_list(self.coeffs.iter().collect());
            for i in 0..self.coeffs.len() {
                self.coeffs[i] = R::div_refs(&self.coeffs[i], &g).unwrap()
            }

            Some((g, self))
        }
    }
}

impl<R: FavoriteAssociate> FavoriteAssociate for Polynomial<R> {
    fn factor_fav_assoc(mut self) -> (Self, Self) {
        if self == Self::zero() {
            (Self::one(), Self::zero())
        } else {
            let (u, _c) = self.coeffs[self.coeffs.len() - 1]
                .clone()
                .factor_fav_assoc();
            for i in 0..self.coeffs.len() {
                self.coeffs[i] = R::div_refs(&self.coeffs[i], &u).unwrap()
            }
            (Self::from(u), self)
        }
    }
}

impl<R: CharacteristicZero> CharacteristicZero for Polynomial<R> {}

impl<R: FiniteUnits> FiniteUnits for Polynomial<R> {
    fn all_units() -> Vec<Self> {
        R::all_units().into_iter().map(|u| Self::from(u)).collect()
    }
}

impl<R: PrincipalIdealDomain> Polynomial<R> {
    //find a polynomial of degree no more than the passed degree
    pub fn interpolate(degree: usize, points: Vec<(R, R)>) -> Option<Self> {
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
        todo!();
    }
}

impl<R: ComRing + Hash> Hash for Polynomial<R> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.coeffs.hash(state);
    }
}

impl<R: UniqueFactorizationDomain> UniqueFactorizationDomain for Polynomial<R> {}

//Kronecker's method for factoring polynomials over infinite ufds with finitely many units (eg the integers)
pub struct KroneckerFactorizer();
impl<R: UniqueFactorizationDomain + FavoriteAssociate + Hash + InfiniteRing + FiniteUnits>
    UniqueFactorizer<Polynomial<R>> for KroneckerFactorizer
{
    fn factor(&mut self, a: &Polynomial<R>) -> Option<UniqueFactorization<Polynomial<R>>> {
        /*
        Suppose we want to factor f(x) = 2 + x + x^2 + x^4 + x^5
        Assume it has a proper factor g(x). wlog g(x) has degree <= 2
        g(x) is determined by its value at 3 points, say at x=0, x=1, x=-1
        f(0)=2, f(1)=6, f(-1)=2     if one of these was zero, then we would have found a linear factor
        g(0) divides 2, g(1) divides 6, g(-1) divides 2
        there are finitely many possible values of g(0), g(1) and g(-1) which satisfy these
        infact there are 4*8*4=128 possible tripples
        however, only 64 need to be checked as the other half are their negatives
        more abstractly, some possibilities can be avoided because we only care about g up to multiplication by a unit
         */

        println!("{}", a.to_string());
        todo!()
    }
}

pub struct FractionFieldPolynomialFactorizer<
    F: FieldOfFractions + Hash,
    RingPolyFactorizerT: UniqueFactorizer<Polynomial<F::R>>,
> where
    F::R: UniqueFactorizationDomain + FavoriteAssociate + Hash,
{
    _f: PhantomData<F>,
    ring_poly_factorizer: RingPolyFactorizerT,
}

impl<F: FieldOfFractions + Hash, RingPolyFactorizerT: UniqueFactorizer<Polynomial<F::R>>>
    FractionFieldPolynomialFactorizer<F, RingPolyFactorizerT>
where
    F::R: UniqueFactorizationDomain + FavoriteAssociate + Hash,
{
    pub fn new(ring_poly_factorizer: RingPolyFactorizerT) -> Self {
        Self {
            _f: PhantomData,
            ring_poly_factorizer,
        }
    }
}

impl<F: FieldOfFractions + Hash, RingPolyFactorizerT: UniqueFactorizer<Polynomial<F::R>>>
    UniqueFactorizer<Polynomial<F>> for FractionFieldPolynomialFactorizer<F, RingPolyFactorizerT>
where
    F::R: UniqueFactorizationDomain + FavoriteAssociate + Hash,
{
    fn factor(&mut self, a: &Polynomial<F>) -> Option<UniqueFactorization<Polynomial<F>>> {
        // F::R::one();
        println!(
            "{:?}",
            self.ring_poly_factorizer
                .factor(&Polynomial::from(F::R::from_int(
                    &malachite_nz::integer::Integer::from(12)
                )))
        );
        todo!()
    }
}

impl<F: Field> EuclideanDomain for Polynomial<F> {
    fn norm(&self) -> Option<Natural> {
        if self == &Self::zero() {
            None
        } else {
            Some(Natural::from(self.coeffs.len() - 1))
        }
    }

    fn quorem(a: Self, b: Self) -> Result<(Self, Self), RingOppErr> {
        Self::quorem_rref(a, &b)
    }

    fn quorem_lref(a: &Self, b: Self) -> Result<(Self, Self), RingOppErr> {
        Self::quorem_refs(a, &b)
    }

    fn quorem_refs(a: &Self, b: &Self) -> Result<(Self, Self), RingOppErr> {
        let res = Self::quorem_rref(a.clone(), b);
        match &res {
            Ok((q, r)) => debug_assert_eq!(&Self::add_ref(Self::mul_refs(q, b), r), a),
            Err(_) => {}
        };
        res
    }

    fn quorem_rref(mut a: Self, b: &Self) -> Result<(Self, Self), RingOppErr> {
        //try to find q such that q*b == a
        // a0 + a1*x + a2*x^2 + ... + am*x^m = (q0 + q1*x + q2*x^2 + ... + qk*x^k) * (b0 + b1*x + b2*x^2 + ... + bn*x^n)
        // 1 + x + x^2 + x^3 + x^4 + x^5 = (?1 + ?x + ?x^2) * (1 + x + x^2 + x^3)      m=6 k=3 n=4
        let m = a.coeffs.len();
        let n = b.coeffs.len();
        if m < n {
            Ok((Self::zero(), a))
        } else if n == 0 {
            Err(RingOppErr::DivideByZero)
        } else {
            let k = m - n + 1;
            let mut q = Self {
                coeffs: (0..k).map(|_i| F::zero()).collect(),
            };
            for i in (0..k).rev() {
                //a[i+n-1] = q[i] * b[n-1]
                match F::div_rref(a.coeff(i + n - 1), &b.coeffs[n - 1]) {
                    Ok(qc) => {
                        //a -= qc*x^i*b
                        a.add_mut(&b.mul_scalar(&qc).mul_var_pow(i).neg());
                        q.coeffs[i] = qc;
                    }
                    Err(_) => panic!(),
                }
            }
            Ok((q, a))
        }
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
        let mut unreduced = Polynomial::<Integer> {
            coeffs: vec![
                Integer::zero(),
                Integer::one(),
                Integer::zero(),
                Integer::zero(),
            ],
        };
        let reduced = Polynomial::<Integer> {
            coeffs: vec![Integer::zero(), Integer::one()],
        };
        unreduced.reduce();
        assert_eq!(unreduced, reduced);

        let mut unreduced = Polynomial::<Integer> {
            coeffs: vec![
                Integer::zero(),
                Integer::zero(),
                Integer::zero(),
                Integer::zero(),
            ],
        };
        let reduced = Polynomial::<Integer> { coeffs: vec![] };
        unreduced.reduce();
        assert_eq!(unreduced, reduced);
    }

    #[test]
    fn divisibility() {
        let x = &Ergonomic::new(Polynomial::<Integer>::var());

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        match Polynomial::<Integer>::div(a.elem(), b.elem()) {
            Ok(c) => {
                println!("{:?} {:?} {:?}", a, b, c);
                assert_eq!(a, b * Ergonomic::new(c))
            }
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) + 1;
        match Polynomial::<Integer>::div(a.elem(), b.elem()) {
            Ok(_c) => panic!(),
            Err(RingOppErr::NotDivisible) => {}
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        match Polynomial::<Integer>::div(a.elem(), b.elem()) {
            Ok(_c) => panic!(),
            Err(RingOppErr::NotDivisible) => {}
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        let b = 0 * x;
        match Polynomial::<Integer>::div(a.elem(), b.elem()) {
            Ok(_c) => panic!(),
            Err(RingOppErr::DivideByZero) => {}
            Err(_) => panic!(),
        }

        let a = 0 * x;
        let b = (x - x) + 5;
        match Polynomial::<Integer>::div(a.elem(), b.elem()) {
            Ok(c) => {
                assert_eq!(c, Polynomial::zero())
            }
            Err(RingOppErr::DivideByZero) => panic!(),
            Err(_) => panic!(),
        }

        let a = 3087 * x - 8805 * x.pow(2) + 607 * x.pow(3) + x.pow(4);
        let b = (x - x) + 1;
        match Polynomial::<Integer>::div(a.elem(), b.elem()) {
            Ok(c) => {
                assert_eq!(c, a.elem())
            }
            Err(RingOppErr::DivideByZero) => panic!(),
            Err(_) => panic!(),
        }
    }

    #[test]
    fn euclidean() {
        let x = &Ergonomic::new(Polynomial::<Rational>::var());

        let a = 1 + x + 3 * x.pow(2) + x.pow(3) + 7 * x.pow(4) + x.pow(5);
        let b = 1 + x + 3 * x.pow(2) + 2 * x.pow(3);
        let (q, r) = Polynomial::<Rational>::quorem_refs(&a.elem(), &b.elem()).unwrap();
        let (q, r) = (Ergonomic::new(q), Ergonomic::new(r));
        println!("{:?} = {:?} * {:?} + {:?}", a, b, q, r);
        assert_eq!(a, &b * &q + &r);

        let a = 3 * x;
        let b = 2 * x;
        let (q, r) = Polynomial::<Rational>::quorem_refs(&a.elem(), &b.elem()).unwrap();
        let (q, r) = (Ergonomic::new(q), Ergonomic::new(r));
        println!("{:?} = {:?} * {:?} + {:?}", a, b, q, r);
        assert_eq!(a, &b * &q + &r);

        let a = 3 * x + 5;
        let b = 2 * x + 1;
        let c = 1 + x + x.pow(2);
        let x = &a * &b;
        let y = &b * &c;
        let g = Polynomial::<Rational>::gcd(x.elem(), y.elem());

        println!("gcd({:?} , {:?}) = {:?}", x, y, g);
        Polynomial::<Rational>::div_refs(&g, &b.elem()).unwrap();
        Polynomial::<Rational>::div_refs(&b.elem(), &g).unwrap();
    }

    #[test]
    fn integer_primitive_and_assoc() {
        let x = &Ergonomic::new(Polynomial::<Integer>::var());
        let p1 = (-2 - 4 * x.pow(2)).elem();
        let (g, p2) = p1.factor_primitive().unwrap();
        assert_eq!(g, Integer::from(2));
        let (u, p3) = p2.factor_fav_assoc();
        assert_eq!(u.coeffs[0], Integer::from(-1));
        assert_eq!(Ergonomic::new(p3), 1 + 2 * x.pow(2));
    }
}
