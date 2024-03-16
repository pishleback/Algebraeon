use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;

use std::fmt::Display;
use std::hash::Hash;

use super::matrix::*;
use super::nzq::*;
use super::ring::*;

#[derive(Debug, Clone)]
pub struct Polynomial<Ring: ComRing> {
    //vec![c0, c1, c2, c3, ..., cn] represents the polynomial c0 + c1*x + c2*x^2 + c3*x^3 + ... + cn * x^n
    //if non-empty, the last item must not be zero
    coeffs: Vec<Ring>,
}

impl<Ring: ComRing> Polynomial<Ring> {
    pub fn coeffs(&self) -> Vec<Ring> {
        self.coeffs.clone()
    }
}

impl<Ring: ComRing> PartialEq for Polynomial<Ring> {
    fn eq(&self, other: &Self) -> bool {
        let n = self.coeffs.len();
        if n != other.coeffs.len() {
            false
        } else {
            (0..n).all(|i| self.coeffs[i] == other.coeffs[i])
        }
    }
}

impl<Ring: ComRing> Eq for Polynomial<Ring> {}

impl<Ring: ComRing + Display> Display for Polynomial<Ring> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.coeffs.len() == 0 {
            write!(f, "0")
        } else {
            let mut first = true;
            for (k, c) in self.coeffs.iter().enumerate() {
                if c != &Ring::zero() {
                    if c == &Ring::one() {
                        if k == 0 {
                            write!(f, "1");
                        } else {
                            write!(f, "+");
                        }
                    } else {
                        if !first {
                            write!(f, "+");
                        }
                        write!(f, "(");
                        write!(f, "{}", c);
                        write!(f, ")");
                    }
                    if k == 0 {
                    } else if k == 1 {
                        write!(f, "λ");
                    } else {
                        write!(f, "λ");
                        write!(f, "^");
                        write!(f, "{}", k);
                    }
                    first = false;
                }
            }
            Ok(())
        }
    }
}

impl<Ring: ComRing> ComRing for Polynomial<Ring> {
    fn zero() -> Self {
        Polynomial { coeffs: vec![] }
    }

    fn one() -> Self {
        Polynomial {
            coeffs: vec![Ring::one()],
        }
    }

    fn neg_mut(&mut self) {
        for coeff in &mut self.coeffs {
            coeff.neg_mut();
        }
    }

    // fn neg_ref(&self) -> Self {
    //     self.clone().neg()
    // }

    // fn neg(mut poly: Self) -> Self {
    //     self.neg_mut(&mut poly);
    //     poly
    // }

    fn add_mut(&mut self, x: &Self) {
        for i in 0..x.coeffs.len() {
            if i < self.coeffs.len() {
                Ring::add_mut(&mut self.coeffs[i], &x.coeffs[i]);
            } else {
                self.coeffs.push(x.coeffs[i].clone());
            }
        }
        self.reduce();
    }

    fn mul_refs(a: &Self, b: &Self) -> Self {
        let mut coeffs = Vec::with_capacity(a.coeffs.len() + b.coeffs.len());
        for _k in 0..a.coeffs.len() + b.coeffs.len() {
            coeffs.push(Ring::zero());
        }
        for i in 0..a.coeffs.len() {
            for j in 0..b.coeffs.len() {
                Ring::add_mut(
                    &mut coeffs[i + j],
                    &Ring::mul_refs(&a.coeffs[i], &b.coeffs[j]),
                );
            }
        }
        let mut ans = Self { coeffs };
        Self::reduce(&mut ans); //TODO: dont have to do this over an integral domain
        ans
    }

    fn mul_mut(&mut self, x: &Self) {
        self.clone_from(&Self::mul_refs(self, x));
    }

    fn div(a: Self, b: Self) -> Result<Self, RingDivisionError> {
        Self::div_rref(a, &b)
    }

    fn div_lref(a: &Self, b: Self) -> Result<Self, RingDivisionError> {
        Self::div_refs(a, &b)
    }

    fn div_refs(a: &Self, b: &Self) -> Result<Self, RingDivisionError> {
        Self::div_rref(a.clone(), b)
    }

    fn div_rref(a: Self, b: &Self) -> Result<Self, RingDivisionError> {
        match Self::try_quorem_rref(a, b) {
            Ok((q, r)) => {
                if r == Self::zero() {
                    Ok(q)
                } else {
                    Err(RingDivisionError::NotDivisible)
                }
            }
            Err(RingDivisionError::NotDivisible) => Err(RingDivisionError::NotDivisible),
            Err(RingDivisionError::DivideByZero) => Err(RingDivisionError::DivideByZero),
        }
    }
}

impl<Ring: ComRing> Polynomial<Ring> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        match self.coeffs.len() {
            0 => {}
            n => {
                if self.coeffs[n - 1] == Ring::zero() {
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
                if self.coeffs[self.coeffs.len() - 1] == Ring::zero() {
                    self.coeffs.pop();
                } else {
                    return;
                }
            }
        }
    }

    pub fn try_quorem(a: Self, b: Self) -> Result<(Self, Self), RingDivisionError> {
        Self::try_quorem_rref(a, &b)
    }

    pub fn try_quorem_lref(a: &Self, b: Self) -> Result<(Self, Self), RingDivisionError> {
        Self::try_quorem_refs(a, &b)
    }

    pub fn try_quorem_refs(a: &Self, b: &Self) -> Result<(Self, Self), RingDivisionError> {
        Self::try_quorem_rref(a.clone(), b)
    }

    pub fn try_quorem_rref(mut a: Self, b: &Self) -> Result<(Self, Self), RingDivisionError> {
        //try to find q such that q*b == a
        // a0 + a1*x + a2*x^2 + ... + am*x^m = (q0 + q1*x + q2*x^2 + ... + qk*x^k) * (b0 + b1*x + b2*x^2 + ... + bn*x^n)
        // 1 + x + x^2 + x^3 + x^4 + x^5 = (?1 + ?x + ?x^2) * (1 + x + x^2 + x^3)      m=6 k=3 n=4

        let m = a.coeffs.len();
        let n = b.coeffs.len();
        if n == 0 {
            Err(RingDivisionError::DivideByZero)
        } else if m < n {
            Ok((Self::zero(), a))
        } else {
            let k = m - n + 1;
            let mut q_coeffs: Vec<_> = (0..k).map(|_i| Ring::zero()).collect();
            for i in (0..k).rev() {
                //a[i+n-1] = q[i] * b[n-1]
                match Ring::div_refs(&a.coeff(i + n - 1), &b.coeff(n - 1)) {
                    Ok(qc) => {
                        //a -= qc*x^i*b
                        Self::add_mut(
                            &mut a,
                            &Self::mul_var_pow(&Self::mul_scalar(&b, &qc), i).neg(),
                        );
                        q_coeffs[i] = qc;
                    }
                    Err(RingDivisionError::NotDivisible) => {
                        return Err(RingDivisionError::NotDivisible);
                    }
                    Err(RingDivisionError::DivideByZero) => panic!(),
                }
            }
            Ok((Self::from_coeffs(q_coeffs), a))
        }
    }

    pub fn apply_map<ImgRing: ComRing>(&self, f: impl Fn(&Ring) -> ImgRing) -> Polynomial<ImgRing> {
        Polynomial::from_coeffs(self.coeffs.iter().map(f).collect())
    }

    pub fn apply_map_with_powers<ImgRing: ComRing>(
        &self,
        f: impl Fn((usize, &Ring)) -> ImgRing,
    ) -> Polynomial<ImgRing> {
        Polynomial::from_coeffs(self.coeffs.iter().enumerate().map(f).collect())
    }

    //find p(q(x))
    pub fn compose(p: &Self, q: &Self) -> Self {
        p.apply_map(|c| Self::constant(c.clone())).evaluate(q)
    }

    //if n = deg(p)
    //return x^n * p(1/x)
    pub fn reversed(&self) -> Self {
        Self::from_coeffs(self.coeffs.clone().into_iter().rev().collect())
    }

    pub fn from_coeffs(coeffs: Vec<Ring>) -> Self {
        let mut p = Polynomial { coeffs };
        p.reduce();
        p
    }

    pub fn constant(x: Ring) -> Self {
        if x == Ring::zero() {
            Polynomial { coeffs: vec![] }
        } else {
            Polynomial { coeffs: vec![x] }
        }
    }

    pub fn constant_var_pow(x: Ring, n: usize) -> Self {
        Polynomial {
            coeffs: (0..n + 1)
                .map(|i| if i < n { Ring::zero() } else { x.clone() })
                .collect(),
        }
    }

    pub fn var() -> Self {
        Polynomial {
            coeffs: vec![Ring::zero(), Ring::one()],
        }
    }

    pub fn var_pow(n: usize) -> Self {
        Polynomial {
            coeffs: (0..n + 1)
                .map(|i| if i < n { Ring::zero() } else { Ring::one() })
                .collect(),
        }
    }

    pub fn mul_var_pow(&self, n: usize) -> Self {
        let mut coeffs = vec![];
        for _i in 0..n {
            coeffs.push(Ring::zero());
        }
        for c in &self.coeffs {
            coeffs.push(c.clone());
        }
        Polynomial { coeffs }
    }

    pub fn mul_scalar(&self, x: &Ring) -> Self {
        let mut ans = Polynomial {
            coeffs: self.coeffs.iter().map(|c| Ring::mul_refs(c, x)).collect(),
        };
        ans.reduce();
        ans
    }

    //getting stuff
    pub fn coeff(&self, i: usize) -> Ring {
        if i < self.coeffs.len() {
            self.coeffs[i].clone()
        } else {
            Ring::zero()
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

    pub fn as_constant(&self) -> Option<Ring> {
        if self.coeffs.len() == 0 {
            Some(Ring::zero())
        } else if self.coeffs.len() == 1 {
            Some(self.coeffs[0].clone())
        } else {
            None
        }
    }

    pub fn is_monic(&self) -> bool {
        match self.degree() {
            Some(d) => self.coeff(d) == Ring::one(),
            None => false,
        }
    }

    pub fn evaluate(&self, x: &Ring) -> Ring {
        // f(x) = a + bx + cx^2 + dx^3
        // evaluate as f(x) = a + x(b + x(c + x(d)))
        let mut y = Ring::zero();
        for c in self.coeffs.iter().rev() {
            Ring::mul_mut(&mut y, x);
            Ring::add_mut(&mut y, c)
        }
        y
    }

    pub fn derivative(mut self) -> Self {
        if self.coeffs.len() > 0 {
            for i in 0..self.coeffs.len() - 1 {
                self.coeffs[i] = self.coeffs[i + 1].clone();
                Ring::mul_mut(&mut self.coeffs[i], &Ring::from_int(&Integer::from(i + 1)));
            }
            self.coeffs.pop();
        }
        self
    }
}

impl<Ring: IntegralDomain> Polynomial<Ring> {
    pub fn pseudorem(a: Self, b: Self) -> Option<Result<Self, &'static str>> {
        Self::pseudorem_rref(a, &b)
    }

    pub fn pseudorem_lref(a: &Self, b: Self) -> Option<Result<Self, &'static str>> {
        Self::pseudorem_refs(a, &b)
    }

    pub fn pseudorem_refs(a: &Self, b: &Self) -> Option<Result<Self, &'static str>> {
        Self::pseudorem_rref(a.clone(), b)
    }

    //None if b = 0
    //error if deg(a) < deg(b)
    pub fn pseudorem_rref(mut a: Self, b: &Self) -> Option<Result<Self, &'static str>> {
        let m = a.coeffs.len();
        let n = b.coeffs.len();

        if n == 0 {
            None
        } else if m < n {
            Some(Err("Should have deg(a) >= deg(b) for pseudo remainder"))
        } else {
            Self::mul_mut(
                &mut a,
                &Self::constant(Ring::nat_pow(
                    &Self::coeff(b, n - 1),
                    &Natural::from(m - n + 1),
                )),
            );

            match Self::try_quorem_rref(a, b) {
                Ok((_q, r)) => Some(Ok(r)),
                Err(_) => panic!(),
            }
        }
    }

    //efficiently compute the gcd of a and b up to scalar multipication using pseudo remainder subresultant sequence
    //the returned polynomial should the smallest non-zero subresultant polynomial
    //https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Trivial_pseudo-remainder_sequence
    pub fn subresultant_gcd(mut a: Self, mut b: Self) -> Self {
        match a.degree() {
            None => b,
            Some(mut a_deg) => match b.degree() {
                None => a,
                Some(mut b_deg) => {
                    if a_deg < b_deg {
                        (a, b) = (b, a);
                        (a_deg, b_deg) = (b_deg, a_deg);
                    }
                    let mut beta = {
                        if (a_deg - b_deg) % 2 == 0 {
                            Ring::from_int(&Integer::from(-1))
                        } else {
                            Ring::from_int(&Integer::from(1))
                        }
                    };
                    let mut psi = Ring::from_int(&Integer::from(-1));
                    loop {
                        let d = a.degree().unwrap() - b.degree().unwrap();
                        let gamma = Self::coeff(&b, b.degree().unwrap());
                        let r = Self::div(
                            Self::pseudorem_rref(a, &b).unwrap().unwrap(),
                            Self::constant(beta),
                        )
                        .unwrap();
                        (a, b) = (b, r);

                        if b == Self::zero() {
                            break;
                        }

                        if d == 0 {
                            //can only happen in the first loop
                            debug_assert_eq!(psi, Ring::one().neg());
                            psi = Ring::one();
                        } else {
                            psi = Ring::div(
                                Ring::nat_pow(&Ring::neg_ref(&gamma), &Natural::from(d)),
                                Ring::nat_pow(&psi, &Natural::from(d - 1)),
                            )
                            .unwrap();
                        }
                        beta = Ring::mul(
                            Ring::neg(gamma),
                            Ring::nat_pow(
                                &psi,
                                &Natural::from(a.degree().unwrap() - b.degree().unwrap()),
                            ),
                        );
                    }
                    a
                }
            },
        }
    }

    pub fn resultant(a: Self, b: Self) -> Ring {
        let subresultant_gcd = Self::subresultant_gcd(a, b);
        match Self::as_constant(&subresultant_gcd) {
            Some(res) => res,
            None => Ring::zero(),
        }
    }
}

impl<Ring: IntegralDomain> IntegralDomain for Polynomial<Ring> {}

// impl<R: UniqueFactorizationDomain> UniqueFactorizationDomain for Polynomial<R> {}

impl<Ring: GreatestCommonDivisorDomain> Polynomial<Ring> {
    pub fn factor_primitive(mut self) -> Option<(Ring, Self)> {
        if self == Self::zero() {
            None
        } else {
            let g = Ring::gcd_list(self.coeffs.iter().collect());
            for i in 0..self.coeffs.len() {
                self.coeffs[i] = Ring::div_refs(&self.coeffs[i], &g).unwrap()
            }
            Some((g, self))
        }
    }

    pub fn primitive_part(self) -> Option<Self> {
        match self.factor_primitive() {
            Some((_unit, prim)) => Some(prim),
            None => None,
        }
    }
}

impl<Ring: GreatestCommonDivisorDomain + CharacteristicZero> Polynomial<Ring> {
    pub fn primitive_squarefree_part(self) -> Self {
        let f = self;
        if f == Self::zero() {
            f
        } else {
            let g = Self::subresultant_gcd(f.clone(), f.clone().derivative());
            let (_c, g_prim) = g.factor_primitive().unwrap();
            let (_c, f_prim) = f.factor_primitive().unwrap();
            let f_prim_sqfree = Self::div(f_prim, g_prim).unwrap();
            f_prim_sqfree
        }
    }
}

impl<Ring: FavoriteAssociate + IntegralDomain> FavoriteAssociate for Polynomial<Ring> {
    fn factor_fav_assoc(mut self) -> (Self, Self) {
        if self == Self::zero() {
            (Self::one(), Self::zero())
        } else {
            let (u, _c) = Ring::factor_fav_assoc(self.coeffs[self.coeffs.len() - 1].clone());
            for i in 0..self.coeffs.len() {
                self.coeffs[i] = Ring::div_refs(&self.coeffs[i], &u).unwrap()
            }
            (Self::constant(u), self)
        }
    }
}

impl<Ring: CharacteristicZero> CharacteristicZero for Polynomial<Ring> {}

impl<Ring: IntegralDomain + FiniteUnits> FiniteUnits for Polynomial<Ring> {
    fn all_units() -> Vec<Self> {
        Ring::all_units()
            .into_iter()
            .map(|u| Self::constant(u))
            .collect()
    }
}

// pub trait InterpolatablePolynomials: ComRing {
//     fn interpolate(points: &Vec<(Self::ElemT, Self::ElemT)>) -> Option<Polynomial<Self>>;
// }

impl<F: Field> EuclideanDomain for Polynomial<F> {
    fn norm(elem: &Self) -> Option<Natural> {
        if elem == &Self::zero() {
            None
        } else {
            Some(Natural::from(elem.coeffs.len() - 1))
        }
    }

    fn quorem(a: Self, b: Self) -> Option<(Self, Self)> {
        Self::quorem_rref(a, &b)
    }

    fn quorem_lref(a: &Self, b: Self) -> Option<(Self, Self)> {
        Self::quorem_refs(a, &b)
    }

    fn quorem_refs(a: &Self, b: &Self) -> Option<(Self, Self)> {
        Self::quorem_rref(a.clone(), b)
    }

    fn quorem_rref(a: Self, b: &Self) -> Option<(Self, Self)> {
        match Self::try_quorem_rref(a, b) {
            Ok((q, r)) => Some((q, r)),
            Err(RingDivisionError::NotDivisible) => panic!(),
            Err(RingDivisionError::DivideByZero) => None,
        }
    }
}

impl<Ring: IntegralDomain> Polynomial<Ring> {
    pub fn interpolate_by_lagrange_basis(points: &Vec<(Ring, Ring)>) -> Option<Self> {
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
        let mut numerator = Self::zero();
        for i in 0..points.len() {
            let (_xi, yi) = &points[i];
            let mut term = Self::constant(yi.clone());

            for j in i + 1..points.len() {
                // (x - xj) for j<i
                let (xj, _yj) = &points[j];
                Self::mul_mut(
                    &mut term,
                    &Self::from_coeffs(vec![Ring::neg_ref(xj), Ring::one()]),
                );
            }
            for j in 0..i {
                // (xj - x) for i<j
                let (xj, _yj) = &points[j];
                Self::mul_mut(
                    &mut term,
                    &Self::from_coeffs(vec![xj.clone(), Ring::neg(Ring::one())]),
                );
            }

            for j in 0..points.len() {
                for k in j + 1..points.len() {
                    if i != j && i != k {
                        let (xj, _yj) = &points[j];
                        let (xk, _yk) = &points[k];
                        Self::mul_mut(
                            &mut term,
                            &Self::constant(Ring::add_ref(Ring::neg_ref(xk), xj)),
                        );
                    }
                }
            }
            Self::add_mut(&mut numerator, &term);
        }

        let mut denominator = Self::one();
        for i in 0..points.len() {
            let (xi, _yi) = &points[i];
            for j in i + 1..points.len() {
                let (xj, _yj) = &points[j];
                Self::mul_mut(
                    &mut denominator,
                    &Self::constant(Ring::add_ref(Ring::neg_ref(xj), xi)),
                );
            }
        }

        match Self::div(numerator, denominator) {
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

impl<Ring: PrincipalIdealDomain> Polynomial<Ring> {
    pub fn interpolate_by_linear_system(points: &Vec<(Ring, Ring)>) -> Option<Self> {
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
        let mut mat = Matrix::<Ring>::zero(n, n);
        for r in 0..n {
            let (x, _y) = &points[r];
            let mut x_pow = Ring::one();
            for c in 0..n {
                *mat.at_mut(r, c).unwrap() = x_pow.clone();
                Ring::mul_mut(&mut x_pow, x);
            }
        }

        let mut output_vec = Matrix::zero(n, 1);
        for r in 0..n {
            let (_x, y) = &points[r];
            *output_vec.at_mut(r, 0).unwrap() = y.clone();
        }

        match mat.col_solve(output_vec) {
            Some(coeff_vec) => Some(Self::from_coeffs(
                (0..n)
                    .map(|i| coeff_vec.at(i, 0).unwrap().clone())
                    .collect(),
            )),
            None => None,
        }
    }
}

impl<Ring: UniqueFactorizationDomain + GreatestCommonDivisorDomain + FiniteUnits> Polynomial<Ring> {
    fn factor_primitive_linear_part(
        &self,
        mut f: Polynomial<Ring>,
    ) -> (Factored<Polynomial<Ring>>, Polynomial<Ring>)
    where
        Polynomial<Ring>: UniqueFactorizationDomain,
    {
        //f should be a primitive polynomial over R
        //linear factor (a+bx) of c0 + c1*x + c2*x^2 + ... + cn*x^n
        //must be such that a divides c0 and b divides cn
        //so just factor c0 and cn and check all divisors

        let mut linear_factors = Factored::factored_one();
        'seek_linear_factor: while f.degree().unwrap() > 0 {
            let c0 = Self::coeff(&f, 0);
            if c0 == Ring::zero() {
                //linear factor of x
                f = Self::div(f, Self::var()).unwrap();
                linear_factors = Factored::mul(
                    linear_factors,
                    Factored::factored_irreducible_unchecked(Self::var()),
                );
                continue 'seek_linear_factor;
            } else {
                //look for linear factors of the form (a+bx)
                let c0fs = Ring::factor(&Self::coeff(&f, 0)).unwrap();
                let cnfs = Ring::factor(&f.coeff(f.degree().unwrap())).unwrap();
                for a_assoc in Ring::divisors(&c0fs) {
                    for u in Ring::all_units() {
                        let a = Ring::mul_ref(u, &a_assoc);
                        for b in Ring::divisors(&cnfs) {
                            //a ranges over all divisors of c0
                            //b ranges over all divisors factors of cn up to associates
                            //try the linear factor (a+bx)
                            let lin = Self::from_coeffs(vec![a.clone(), b]);
                            match Self::div_refs(&f, &lin) {
                                Ok(new_f) => {
                                    f = new_f;
                                    linear_factors = Factored::mul(
                                        linear_factors,
                                        Factored::factored_irreducible_unchecked(lin),
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
        Ring: UniqueFactorizationDomain + GreatestCommonDivisorDomain + CharacteristicZero + FiniteUnits,
    > Polynomial<Ring>
{
    fn factorize_primitive_polynomial_by_kroneckers_method(&self, f: Self) -> Factored<Self>
    where
        Self: UniqueFactorizationDomain,
    {
        debug_assert_ne!(f, Self::zero());
        fn partial_factor<
            'a,
            Ring: UniqueFactorizationDomain
                + GreatestCommonDivisorDomain
                + CharacteristicZero
                + FiniteUnits,
        >(
            f: Polynomial<Ring>,
        ) -> Option<(Polynomial<Ring>, Polynomial<Ring>)> {
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
            let f_deg = f.degree().unwrap();
            if f_deg == 1 {
                //linear factor is irreducible
                None
            } else {
                let max_factor_degree = f_deg / 2;
                let mut f_points = vec![];
                let mut elem_gen = Ring::generate_distinct_elements();
                //take more samples than necessary, then take the subset with the smallest number of divisors
                while f_points.len() < 3 * (max_factor_degree + 1) {
                    //loop terminates because polynomial over integral domain has finitely many roots
                    let x = elem_gen.next().unwrap();
                    let y = Polynomial::evaluate(&f, &x);
                    if y != Ring::zero() {
                        f_points.push((x, Ring::factor(&y).unwrap()));
                    }
                }

                //compute all factors of each y value. choose the y with the most divisors to only factor up to units
                f_points.sort_by_cached_key(|(_x, yf)| Ring::count_divisors(yf));
                let _ = f_points.split_off(max_factor_degree + 1);
                //possible_g_points is (x, possible_y_values)
                let all_possible_g_points: Vec<(Ring, Vec<Ring>)> = f_points
                    .into_iter()
                    .rev()
                    .enumerate()
                    .map(|(i, (x, yf))| {
                        let mut y_divs = vec![];
                        for d in Ring::divisors(&yf) {
                            if i == 0 {
                                //take divisors up to associates for one, because we only care about g up to associates
                                y_divs.push(d);
                            } else {
                                //take _all_ divisors for the rest
                                for u in Ring::all_units() {
                                    y_divs.push(Ring::mul_ref(u, &d));
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
                    match Polynomial::interpolate_by_lagrange_basis(&possible_g_points) {
                        Some(g) => {
                            if g.degree().unwrap() >= 1 {
                                //g is a possible proper divisor of f
                                match Polynomial::div_refs(&f, &g) {
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
        full_factor_using_partial_factor(f, &partial_factor::<Ring>)
    }

    pub fn factorize_by_kroneckers_method(&self) -> Option<Factored<Self>>
    where
        Self: UniqueFactorizationDomain,
    {
        let f = self;
        if f == &Self::zero() {
            return None;
        }

        // println!();
        // println!("factorize_by_kroneckers_method");
        // println!("f = {:?}", f);
        let (scalar_part, f) = f.clone().factor_primitive().unwrap();
        // println!("scalar_part = {:?}", scalar_part);
        // println!("primitive_part = {:?}", f);
        let factored_scalar_part = Ring::factor(&scalar_part).unwrap();
        // println!("factored_scalar_part = {:?}", factored_scalar_part);
        let factored_scalar_part_poly = Factored::new_unchecked(
            Self::constant(factored_scalar_part.unit().clone()),
            factored_scalar_part
                .factors()
                .iter()
                .map(|(p, k)| (Self::constant(p.clone()), k.clone()))
                .collect(),
        );
        let (linear_factors, f) = self.factor_primitive_linear_part(f);

        let f_sqfree = f.clone().primitive_squarefree_part();
        let f_sqfree_factors = self.factorize_primitive_polynomial_by_kroneckers_method(f_sqfree);

        let mut f = f;
        let mut f_factors = Factored::factored_one();
        for (sqfree_factor, k) in f_sqfree_factors.factors() {
            debug_assert_eq!(k, &Natural::from(1u8)); //the factorization should be squarefree
            loop {
                match Self::div_refs(&f, sqfree_factor) {
                    Ok(new_f) => {
                        f = new_f;
                        f_factors = Factored::mul(
                            f_factors,
                            Factored::new_unchecked(
                                Self::one(),
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
        debug_assert_eq!(f.degree().unwrap(), 0);
        debug_assert!(Self::is_unit(f.clone()));
        f_factors = Factored::mul(f_factors, Factored::factored_unit_unchecked(f));

        return Some(Factored::mul(
            Factored::mul(linear_factors, factored_scalar_part_poly),
            f_factors,
        ));
    }
}

impl<Ring: Field + FiniteUnits> Polynomial<Ring> {
    pub fn factorize_by_trying_all_factors(self) -> Option<Factored<Polynomial<Ring>>>
    where
        Self: UniqueFactorizationDomain,
    {
        let f = self;
        if f == Self::zero() {
            None
        } else {
            fn partial_factor<'a, Ring: Field + FiniteUnits>(
                f: Polynomial<Ring>,
            ) -> Option<(Polynomial<Ring>, Polynomial<Ring>)> {
                let f_deg = f.degree().unwrap();
                let max_factor_degree = f_deg / 2;
                for d in 0..max_factor_degree {
                    for mut coeffs in itertools::Itertools::multi_cartesian_product(
                        (0..d + 1).into_iter().map(|_d| {
                            let mut all_elems = vec![Ring::zero()];
                            all_elems.append(&mut Ring::all_units());
                            all_elems
                        }),
                    ) {
                        coeffs.push(Ring::one());
                        let g = Polynomial::from_coeffs(coeffs);
                        // println!("{}", self.to_string(&g));
                        match Polynomial::div_refs(&f, &g) {
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
            Some(full_factor_using_partial_factor(f, &partial_factor::<Ring>))
        }
    }
}

impl<F: FieldOfFractions> Polynomial<F>
where
    F::R: EuclideanDomain + FavoriteAssociate,
{
    pub fn factor_primitive_fof(&self) -> (F, Polynomial<F::R>) {
        let div = F::R::lcm_list(self.coeffs.iter().map(|c| F::denominator(&c)).collect());

        let (mul, prim) = self
            .apply_map(|c| F::as_base_ring(F::mul_ref(F::from_base_ring(div.clone()), c)).unwrap())
            .factor_primitive()
            .unwrap();

        (
            F::div(F::from_base_ring(mul), F::from_base_ring(div)).unwrap(),
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

/*
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
    //https://github.com/JeffreyDF/LLL_factorization/blob/master/factorization.py

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

    //check that g has invertible leading coefficient
    debug_assert!(old_ring.is_unit(old_poly_ring.coeff(&h, old_poly_ring.degree(&h).unwrap())));

    //check that h is monic
    debug_assert!(old_ring.equal(
        &old_poly_ring.coeff(&h, old_poly_ring.degree(&h).unwrap()),
        &old_ring.one()
    ));

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

    // Lifting the factorization:

    // (f - g*h)
    let e = new_poly_ring.add(f, new_poly_ring.neg(new_poly_ring.mul_refs(&g, &h)));
    // Division of se by h in Z/m^2. Doesnt fail because h is monic
    let (q, r) = new_poly_ring
        .try_quorem_rref(new_poly_ring.mul_refs(&s, &e), &h)
        .unwrap();
    let (g_, h_) = (
        new_poly_ring.sum(vec![
            &g.clone(),
            &new_poly_ring.mul_ref(e, &t),
            &new_poly_ring.mul_ref(q, &g),
        ]),
        new_poly_ring.add(h, r),
    );

    //Lifting the Bezout coefficients:
    let b = new_poly_ring.sum(vec![
        &new_poly_ring.mul_refs(&g_, &s),
        &new_poly_ring.mul_refs(&h_, &t),
        &new_poly_ring.neg(new_poly_ring.one()),
    ]);
    let (c, d) = new_poly_ring
        .try_quorem_rref(new_poly_ring.mul_refs(&b, &s), &h_)
        .unwrap();
    let (s_, t_) = (
        new_poly_ring.add(s, new_poly_ring.neg(d)),
        new_poly_ring.sum(vec![
            &t.clone(),
            &new_poly_ring.neg(new_poly_ring.mul(t, b)),
            &new_poly_ring.neg(new_poly_ring.mul(c, g)),
        ]),
    );

    // Return the polynomials as embedded in Z[x].
    (g_, h_, s_, t_)
}
*/

// impl<'a> PolynomialRing<'a, IntegerRing> {
//     pub fn factorize_by_hensel_lifting(
//         &self,
//         f: Polynomial<Integer>,
//     ) -> Factored<Polynomial<Integer>> {
//         debug_assert!(!self.equal(&f, &self.zero()));

//         println!("{:?}", ZZ_POLY.factorize_by_kroneckers_method(&f));

//         let pgen = NaturalPrimeGenerator::new();
//         for p in pgen {
//             let modp = EuclideanQuotient::new_field_unchecked(ZZ, Integer::from(&p));
//             let poly_modp = PolynomialRing::new(&modp);
//             let gs = poly_modp.factorize_by_trying_all_factors(f.clone());

//             if gs.factors().iter().all(|(_g, k)| k == &Natural::from(1u8)) {
//                 let gs: Vec<_> = gs.factors().iter().map(|(g, _k)| g).collect();
//                 println!("{:?}", p);
//                 for g in gs.into_iter() {
//                     println!();
//                     let h = poly_modp.div_refs(&f, g).unwrap();
//                     println!(
//                         "gh = {}    {}",
//                         poly_modp.to_string(g),
//                         poly_modp.to_string(&h)
//                     );

//                     let (u, s, t) = poly_modp.xgcd(g.clone(), h);
//                     debug_assert!(poly_modp.equal(&u, &poly_modp.one()));

//                     println!(
//                         "st = {}    {}",
//                         poly_modp.to_string(&s),
//                         poly_modp.to_string(&t)
//                     );
//                 }

//                 todo!();
//             }
//         }
//         panic!()
//     }
// }

// impl PartialEq for Polynomial<Integer> {
//     fn eq(&self, other: &Self) -> bool {
//         self.coeffs == other.coeffs
//     }
// }

// impl Eq for Polynomial<Integer> {}

impl<Ring: ComRing + Hash> Hash for Polynomial<Ring> {
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
        unreduced.reduce();
        assert_eq!(unreduced, reduced);

        let mut unreduced = Polynomial {
            coeffs: vec![
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
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
        match Polynomial::div(a.elem(), b.elem()) {
            Ok(c) => {
                println!("{:?} {:?} {:?}", a, b, c);
                assert_eq!(a, b * Ergonomic::new(c))
            }
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) + 1;
        match Polynomial::div(a.elem(), b.elem()) {
            Ok(_c) => panic!(),
            Err(RingDivisionError::NotDivisible) => {}
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        match Polynomial::div(a.elem(), b.elem()) {
            Ok(_c) => panic!(),
            Err(RingDivisionError::NotDivisible) => {}
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        let b = 0 * x;
        match Polynomial::div(a.elem(), b.elem()) {
            Ok(_c) => panic!(),
            Err(RingDivisionError::DivideByZero) => {}
            Err(_) => panic!(),
        }

        let a = 0 * x;
        let b = (x - x) + 5;
        match Polynomial::div(a.elem(), b.elem()) {
            Ok(c) => {
                assert_eq!(c, Polynomial::zero())
            }
            Err(RingDivisionError::DivideByZero) => panic!(),
            Err(_) => panic!(),
        }

        let a = 3087 * x - 8805 * x.pow(2) + 607 * x.pow(3) + x.pow(4);
        let b = (x - x) + 1;
        match Polynomial::div(a.elem(), b.elem()) {
            Ok(c) => {
                assert_eq!(c, a.elem())
            }
            Err(RingDivisionError::DivideByZero) => panic!(),
            Err(_) => panic!(),
        }
    }

    #[test]
    fn euclidean() {
        let x = &Ergonomic::new(Polynomial::<Rational>::var());

        let a = 1 + x + 3 * x.pow(2) + x.pow(3) + 7 * x.pow(4) + x.pow(5);
        let b = 1 + x + 3 * x.pow(2) + 2 * x.pow(3);
        let (q, r) = Polynomial::quorem_refs(&a.elem(), &b.elem()).unwrap();
        let (q, r) = (Ergonomic::new(q), Ergonomic::new(r));
        println!("{:?} = {:?} * {:?} + {:?}", a, b, q, r);
        assert_eq!(a, &b * &q + &r);

        let a = 3 * x;
        let b = 2 * x;
        let (q, r) = Polynomial::quorem_refs(&a.elem(), &b.elem()).unwrap();
        let (q, r) = (Ergonomic::new(q), Ergonomic::new(r));
        println!("{:?} = {:?} * {:?} + {:?}", a, b, q, r);
        assert_eq!(a, &b * &q + &r);

        let a = 3 * x + 5;
        let b = 2 * x + 1;
        let c = 1 + x + x.pow(2);
        let x = &a * &b;
        let y = &b * &c;
        let g = Polynomial::gcd(x.elem(), y.elem());

        println!("gcd({:?} , {:?}) = {:?}", x, y, g);
        Polynomial::div_refs(&g, &b.elem()).unwrap();
        Polynomial::div_refs(&b.elem(), &g).unwrap();
    }

    #[test]
    fn test_pseudo_remainder() {
        let x = &Ergonomic::new(Polynomial::<Integer>::var());
        {
            let f = (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5)
                .elem();
            let g = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).elem();

            println!("f = {}", f.to_string());
            println!("g = {}", g.to_string());

            let r1 = Polynomial::pseudorem_refs(&f, &g).unwrap().unwrap();
            println!("r1 = {}", r1.to_string());
            assert_eq!(r1, (-15 * x.pow(4) + 3 * x.pow(2) - 9).elem());

            let r2 = Polynomial::pseudorem_refs(&g, &r1).unwrap().unwrap();
            println!("r2 = {}", r2.to_string());
            assert_eq!(r2, (15795 * x.pow(2) + 30375 * x - 59535).elem());
        }
        println!();
        {
            let f = (4 * x.pow(3) + 2 * x - 7).elem();
            let g = Polynomial::zero();

            println!("f = {}", f.to_string());
            println!("g = {}", g.to_string());

            if let None = Polynomial::pseudorem_refs(&f, &g) {
            } else {
                assert!(false);
            }
        }
        println!();
        {
            let f = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).elem();
            let g = (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5)
                .elem();

            println!("f = {}", f.to_string());
            println!("g = {}", g.to_string());

            if let Err(_msg) = Polynomial::pseudorem_refs(&f, &g).unwrap() {
            } else {
                assert!(false);
            }
        }
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

    #[test]
    fn test_evaluate() {
        let x = &Ergonomic::new(Polynomial::<Integer>::var());
        let f = (1 + x + 3 * x.pow(2) + x.pow(3) + 7 * x.pow(4) + x.pow(5)).elem();
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
    fn test_factor_by_kroneckers_method_over_integers() {
        let x = &Ergonomic::new(Polynomial::<Integer>::var());

        //primitive cases
        let f = ((1 + x).pow(2)).elem();
        assert!(Factored::equal(
            &Polynomial::factorize_by_kroneckers_method(&f).unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![((1 + x).elem(), Natural::from(2u8))]
            )
        ));

        let f = (-1 - 2 * x).elem();
        let fs1 = f.factorize_by_kroneckers_method().unwrap();
        let fs2 = &Factored::new_unchecked(
            Polynomial::neg(Polynomial::one()),
            vec![((1 + 2 * x).elem(), Natural::from(1u8))],
        );
        println!("fs1={:?} fs2={:?}", fs1, fs2);
        assert!(Factored::equal(&fs1, &fs2));

        let f = (x.pow(5) + x.pow(4) + x.pow(2) + x + 2).elem();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![
                    ((1 + x + x.pow(2)).elem(), Natural::from(1u8)),
                    ((2 - x + x.pow(3)).elem(), Natural::from(1u8))
                ]
            )
        ));

        let f = (1 + x + x.pow(2)).pow(2).elem();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![((1 + x + x.pow(2)).elem(), Natural::from(2u8))]
            )
        ));

        //non-primitive cases
        let f = (2 + 2 * x).elem();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![
                    (Polynomial::from_int(&Integer::from(2)), Natural::from(1u8)),
                    ((1 + x).elem(), Natural::from(1u8))
                ]
            )
        ));

        let f = (12 * (2 + 3 * x) * (x - 1).pow(2)).elem();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![
                    (Polynomial::from_int(&Integer::from(2)), Natural::from(2u8)),
                    (Polynomial::from_int(&Integer::from(3)), Natural::from(1u8)),
                    ((2 + 3 * x).elem(), Natural::from(1u8)),
                    ((x - 1).elem(), Natural::from(2u8))
                ]
            )
        ));

        let f = Polynomial::<Integer>::one();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(Polynomial::one(), vec![])
        ));

        let f = ((x.pow(4) + x + 1) * (x.pow(3) + x + 1)).elem();
        assert!(Factored::equal(
            &f.factorize_by_kroneckers_method().unwrap(),
            &Factored::new_unchecked(
                Polynomial::one(),
                vec![
                    ((x.pow(4) + x + 1).elem(), Natural::from(1u8)),
                    ((x.pow(3) + x + 1).elem(), Natural::from(1u8))
                ]
            )
        ));
    }

    #[test]
    fn test_derivative() {
        let x = &Ergonomic::new(Polynomial::<Integer>::var());
        let f = (2 + 3 * x - x.pow(2) + 7 * x.pow(3)).elem();
        let g = (3 - 2 * x + 21 * x.pow(2)).elem();
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
        let x = &Ergonomic::new(Polynomial::<Integer>::var());
        let f = (2 + 3 * x - x.pow(2) + 7 * x.pow(3) + x.pow(4)).elem();
        let g = (3 - 2 * x + 21 * x.pow(2)).elem();
        assert!(f.is_monic());
        assert!(!g.is_monic());
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
        let x = &Ergonomic::new(Polynomial::<Integer>::var());

        let f =
            (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5).elem();
        let g = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).elem();
        assert_eq!(
            Polynomial::subresultant_gcd(f, g),
            Polynomial::constant(Integer::from(260708))
        );

        let f = (3 * x.pow(6) + 5 * x.pow(4) - 4 * x.pow(2) - 9 * x + 21).elem();
        let g =
            (x.pow(8) + x.pow(6) - 3 * x.pow(4) - 3 * x.pow(3) + 8 * x.pow(2) + 2 * x - 5).elem();
        assert_eq!(
            Polynomial::subresultant_gcd(f, g),
            Polynomial::constant(Integer::from(260708))
        );

        let f = ((x + 2).pow(2) * (2 * x - 3).pow(2)).elem();
        let g = ((3 * x - 1) * (2 * x - 3).pow(2)).elem();
        assert_eq!(
            Polynomial::subresultant_gcd(f, g),
            (7056 - 9408 * x + 3136 * x.pow(2)).elem()
        );
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
                Polynomial::from_coeffs(vec![
                    Rational::from_signeds(1, 2),
                    Rational::from_signeds(1, 3),
                ]),
                Polynomial::from_coeffs(vec![Integer::from(3), Integer::from(2)]),
            ),
            (
                Polynomial::from_coeffs(vec![
                    Rational::from_signeds(4, 1),
                    Rational::from_signeds(6, 1),
                ]),
                Polynomial::from_coeffs(vec![Integer::from(2), Integer::from(3)]),
            ),
        ] {
            let (mul, ans) = f.factor_primitive_fof();
            assert!(Polynomial::are_associate(&ans, &exp));
            assert_eq!(
                Polynomial::mul_scalar(&ans.apply_map(|c| Rational::from(c)), &mul),
                f
            );
        }
    }
}
