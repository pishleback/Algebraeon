#![allow(dead_code)]

use malachite_nz::integer::Integer;
use malachite_nz::natural::Natural;

use std::hash::Hash;

use super::matrix::*;
use super::ring::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Polynomial<R: ComRing> {
    //vec![c0, c1, c2, c3, ..., cn] represents the polynomial c0 + c1*x + c2*x^2 + c3*x^3 + ... + cn * x^n
    //if non-empty, the last item must not be zero
    coeffs: Vec<R>,
}

impl<R: ComRing> ToString for Polynomial<R> {
    fn to_string(&self) -> String {
        if self == &Self::zero() {
            String::from("0")
        } else {
            let mut s = String::new();
            let mut first = true;
            for (k, c) in self.coeffs.iter().enumerate() {
                if c != &R::zero() {
                    if c == &R::one() {
                        if k == 0 {
                            s += "1";
                        } else {
                            s += "+"
                        }
                    } else if c == &R::one().neg() {
                        if k == 0 {
                            s += "-1";
                        } else {
                            s += "-"
                        }
                    } else {
                        if !first {
                            s += "+";
                        }
                        s += "(";
                        s += &c.to_string();
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

    pub fn evaluate(&self, x: &R) -> R {
        // f(x) = a + bx + cx^2 + dx^3
        // evaluate as f(x) = a + x(b + x(c + x(d)))
        let mut y = R::zero();
        for c in self.coeffs.iter().rev() {
            y.mul_mut(x);
            y.add_mut(c)
        }
        y
    }

    pub fn apply_map<S: ComRing>(&self, f: impl Fn(&R) -> S) -> Polynomial<S> {
        Polynomial::new(self.coeffs.iter().map(f).collect())
    }

    pub fn derivative(mut self) -> Self {
        if self.coeffs.len() > 0 {
            for i in 0..self.coeffs.len() - 1 {
                self.coeffs[i] = self.coeffs[i + 1].clone();
                self.coeffs[i].mul_mut(&R::from_int(&Integer::from(i + 1)));
            }
            self.coeffs.pop();
        }
        self
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

impl<R: ComRing> From<&R> for Polynomial<R> {
    fn from(x: &R) -> Self {
        if x == &R::zero() {
            Self { coeffs: vec![] }
        } else {
            Self {
                coeffs: vec![x.clone()],
            }
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

    fn div(a: Self, b: Self) -> Result<Self, RingDivisionError> {
        Self::div_rref(a, &b)
    }

    fn div_lref(a: &Self, b: Self) -> Result<Self, RingDivisionError> {
        Self::div_refs(a, &b)
    }

    fn div_refs(a: &Self, b: &Self) -> Result<Self, RingDivisionError> {
        let q_res = Self::div_rref(a.clone(), b);
        match &q_res {
            Ok(q) => debug_assert_eq!(&Self::mul_refs(q, b), a),
            Err(_) => {}
        };
        q_res
    }

    fn div_rref(mut a: Self, b: &Self) -> Result<Self, RingDivisionError> {
        //try to find q such that q*b == a
        // a0 + a1*x + a2*x^2 + ... + am*x^m = (q0 + q1*x + q2*x^2 + ... + qk*x^k) * (b0 + b1*x + b2*x^2 + ... + bn*x^n)
        // 1 + x + x^2 + x^3 + x^4 + x^5 = (?1 + ?x + ?x^2) * (1 + x + x^2 + x^3)      m=6 k=3 n=4
        let m = a.coeffs.len();
        let n = b.coeffs.len();
        if n == 0 {
            Err(RingDivisionError::DivideByZero)
        } else if m == 0 {
            Ok(Self::zero())
        } else if m < n {
            Err(RingDivisionError::NotDivisible)
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
                    Err(RingDivisionError::NotDivisible) => {
                        return Err(RingDivisionError::NotDivisible);
                    }
                    Err(_) => panic!(),
                }
            }
            if a != Self::zero() {
                return Err(RingDivisionError::NotDivisible);
            }
            Ok(q)
        }
    }
}

impl<R: IntegralDomain> IntegralDomain for Polynomial<R> {}

impl<R: UniqueFactorizationDomain> UniqueFactorizationDomain for Polynomial<R> {}

impl<R: GreatestCommonDivisorDomain> Polynomial<R> {
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

impl<R: FavoriteAssociate + IntegralDomain> FavoriteAssociate for Polynomial<R> {
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

impl<R: ComRing + Hash> Hash for Polynomial<R> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.coeffs.hash(state);
    }
}

impl<R: IntegralDomain + FiniteUnits> FiniteUnits for Polynomial<R> {
    fn all_units() -> Vec<Self> {
        R::all_units().into_iter().map(|u| Self::from(u)).collect()
    }
}

pub trait InterpolatablePolynomials: ComRing {
    fn interpolate(points: &Vec<(Self, Self)>) -> Option<Polynomial<Self>>;
}

impl<F: Field> EuclideanDomain for Polynomial<F> {
    fn norm(&self) -> Option<Natural> {
        if self == &Self::zero() {
            None
        } else {
            Some(Natural::from(self.coeffs.len() - 1))
        }
    }

    fn quorem(a: Self, b: Self) -> Option<(Self, Self)> {
        Self::quorem_rref(a, &b)
    }

    fn quorem_lref(a: &Self, b: Self) -> Option<(Self, Self)> {
        Self::quorem_refs(a, &b)
    }

    fn quorem_refs(a: &Self, b: &Self) -> Option<(Self, Self)> {
        let res = Self::quorem_rref(a.clone(), b);
        match &res {
            Some((q, r)) => debug_assert_eq!(&Self::add_ref(Self::mul_refs(q, b), r), a),
            None => {}
        };
        res
    }

    fn quorem_rref(mut a: Self, b: &Self) -> Option<(Self, Self)> {
        //try to find q such that q*b == a
        // a0 + a1*x + a2*x^2 + ... + am*x^m = (q0 + q1*x + q2*x^2 + ... + qk*x^k) * (b0 + b1*x + b2*x^2 + ... + bn*x^n)
        // 1 + x + x^2 + x^3 + x^4 + x^5 = (?1 + ?x + ?x^2) * (1 + x + x^2 + x^3)      m=6 k=3 n=4
        let m = a.coeffs.len();
        let n = b.coeffs.len();

        if n == 0 {
            None
        } else if m < n {
            Some((Self::zero(), a))
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
            Some((q, a))
        }
    }
}

pub fn interpolate_by_lagrange_basis<R: IntegralDomain>(
    points: &Vec<(R, R)>,
) -> Option<Polynomial<R>> {
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

    let mut numerator = Polynomial::zero();
    for i in 0..points.len() {
        let (_xi, yi) = &points[i];
        let mut term = Polynomial::from(yi);

        for j in i + 1..points.len() {
            // (x - xj) for j<i
            let (xj, _yj) = &points[j];
            term.mul_mut(&Polynomial::new(vec![xj.neg_ref(), R::one()]));
        }
        for j in 0..i {
            // (xj - x) for i<j
            let (xj, _yj) = &points[j];
            term.mul_mut(&Polynomial::new(vec![xj.clone(), R::one().neg()]));
        }

        for j in 0..points.len() {
            for k in j + 1..points.len() {
                if i != j && i != k {
                    let (xj, _yj) = &points[j];
                    let (xk, _yk) = &points[k];
                    term.mul_mut(&Polynomial::from(R::add_ref(xk.neg_ref(), xj)));
                }
            }
        }
        numerator.add_mut(&term);
    }

    let mut denominator = Polynomial::one();
    for i in 0..points.len() {
        let (xi, _yi) = &points[i];
        for j in i + 1..points.len() {
            let (xj, _yj) = &points[j];
            denominator.mul_mut(&Polynomial::from(R::add_ref(xj.neg_ref(), xi)));
        }
    }

    match Polynomial::div(numerator, denominator) {
        Ok(interp_poly) => Some(interp_poly),
        Err(RingDivisionError::NotDivisible) => None,
        Err(RingDivisionError::DivideByZero) => {
            panic!("are the input points distinct?")
        }
    }
}

pub fn interpolate_by_linear_system<R: PrincipalIdealDomain>(
    points: &Vec<(R, R)>,
) -> Option<Polynomial<R>> {
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
    let mut mat = Matrix::zero(n, n);
    for r in 0..n {
        let (x, _y) = &points[r];
        let mut x_pow = R::one();
        for c in 0..n {
            *mat.at_mut(r, c).unwrap() = x_pow.clone();
            x_pow.mul_mut(x);
        }
    }

    let mut output_vec = Matrix::zero(n, 1);
    for r in 0..n {
        let (_x, y) = &points[r];
        *output_vec.at_mut(r, 0).unwrap() = y.clone();
    }

    match mat.col_solve(output_vec) {
        Some(coeff_vec) => Some(Polynomial::new(
            (0..n)
                .map(|i| coeff_vec.at(i, 0).unwrap().clone())
                .collect(),
        )),
        None => None,
    }
}

fn factor_primitive_linear_part<
    R: UniquelyFactorable + GreatestCommonDivisorDomain + FiniteUnits,
>(
    mut f: Polynomial<R>,
) -> (Factored<Polynomial<R>>, Polynomial<R>) {
    //f should be a primitive polynomial over R
    //linear factor (a+bx) of c0 + c1*x + c2*x^2 + ... + cn*x^n
    //must be such that a divides c0 and b divides cn
    //so just factor c0 and cn and check all divisors

    let mut linear_factors = Factored::one();
    'seek_linear_factor: while f.degree().unwrap() > 0 {
        let c0 = f.coeff(0);
        if c0 == R::zero() {
            //linear factor of x
            f = Polynomial::div(f, Polynomial::var()).unwrap();
            linear_factors = Factored::mul(
                linear_factors,
                Factored::new_irreducible_unchecked(Polynomial::var()),
            );
            continue 'seek_linear_factor;
        } else {
            //look for linear factors of the form (a+bx)
            let c0fs = f.coeff(0).factor().unwrap();
            let cnfs = f.coeff(f.degree().unwrap()).factor().unwrap();
            for a_assoc in c0fs.divisors() {
                for u in R::all_units() {
                    let a = R::mul_ref(u, &a_assoc);
                    for b in cnfs.divisors() {
                        //a ranges over all divisors of c0
                        //b ranges over all divisors factors of cn up to associates
                        //try the linear factor (a+bx)
                        let lin = Polynomial::new(vec![a.clone(), b]);
                        match Polynomial::div_refs(&f, &lin) {
                            Ok(new_f) => {
                                f = new_f;
                                linear_factors = Factored::mul(
                                    linear_factors,
                                    Factored::new_irreducible_unchecked(lin),
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

// pub fn squarefree_part_by_yuns<R: ComRing>(f: &Polynomial<R>) -> Polynomial<R>
// where
//     Polynomial<R>: GreatestCommonDivisorDomain,
// {
//     todo!();
//     //need fast gcd over the integers
//     //need pseudo-remainder sequences
// }

pub fn factor_by_kroneckers_method<
    R: UniquelyFactorable
        + GreatestCommonDivisorDomain
        + InfiniteRing
        + FiniteUnits
        + InterpolatablePolynomials,
>(
    f: &Polynomial<R>,
) -> Option<Factored<Polynomial<R>>>
where
    Polynomial<R>: UniqueFactorizationDomain,
{
    if f == &Polynomial::<R>::zero() {
        return None;
    }

    let (scalar_part, f) = f.clone().factor_primitive().unwrap();
    let factored_scalar_part = scalar_part.factor().unwrap();
    let factored_scalar_part_poly = Factored::new_unchecked(
        Polynomial::from(scalar_part),
        Polynomial::from(factored_scalar_part.unit()),
        factored_scalar_part
            .factors()
            .iter()
            .map(|(p, k)| (Polynomial::from(p), k.clone()))
            .collect(),
    );
    let (linear_factors, f) = factor_primitive_linear_part(f);
    return Some(Factored::mul(
        linear_factors,
        Factored::mul(factored_scalar_part_poly, factor_primitive(f)),
    ));

    fn factor_primitive<
        R: UniquelyFactorable
            + GreatestCommonDivisorDomain
            + InfiniteRing
            + FiniteUnits
            + InterpolatablePolynomials,
    >(
        f: Polynomial<R>,
    ) -> Factored<Polynomial<R>>
    where
        Polynomial<R>: UniqueFactorizationDomain,
    {
        debug_assert_ne!(f, Polynomial::<R>::zero());
        let f_deg = f.degree().unwrap();
        if f_deg == 0 {
            debug_assert!(f.coeff(0).is_unit());
            Factored::new_unit_unchecked(f)
        } else if f_deg == 1 {
            //linear primitive factors are always irreducible
            Factored::new_irreducible_unchecked(f.clone())
        } else {
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

            let max_factor_degree = f_deg / 2;
            let mut f_points = vec![];
            let mut elem_gen = R::generate_distinct_elements();
            //take more samples than necessary, then take the subset with the smallest number of divisors
            while f_points.len() < 3 * (max_factor_degree + 1) {
                //loop terminates because polynomial over integral domain has finitely many roots
                let x = elem_gen.next().unwrap();
                let y = f.evaluate(&x);
                if y != R::zero() {
                    f_points.push((x, y.factor().unwrap()));
                }
            }

            //compute all factors of each y value. choose the y with the most divisors to only factor up to units
            f_points.sort_by_cached_key(|(_x, yf)| yf.count_divisors());
            let _ = f_points.split_off(max_factor_degree + 1);
            //possible_g_points is (x, possible_y_values)
            let all_possible_g_points: Vec<(R, Vec<R>)> = f_points
                .into_iter()
                .rev()
                .enumerate()
                .map(|(i, (x, yf))| {
                    let mut y_divs = vec![];
                    for d in yf.divisors() {
                        if i == 0 {
                            //take divisors up to associates for one, because we only care about g up to associates
                            y_divs.push(d);
                        } else {
                            //take _all_ divisors for the rest
                            for u in R::all_units() {
                                y_divs.push(R::mul_ref(u, &d));
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
                match R::interpolate(&possible_g_points) {
                    Some(g) => {
                        if g.degree().unwrap() >= 1 {
                            //g is a possible proper divisor of f
                            match Polynomial::div_refs(&f, &g) {
                                Ok(h) => {
                                    //g really is a proper divisor of f
                                    return Factored::mul(factor_primitive(g), factor_primitive(h));
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
            Factored::new_irreducible_unchecked(f.clone())
        }
    }
}

// pub fn factor_over_field_of_fractions<F: FieldOfFractions>(
//     poly: &Polynomial<F>,
//     factor_over_base_ring: Box<dyn Fn(Polynomial<F::R>) -> Factored<Polynomial<F::R>>>,
// ) -> Factored<Polynomial<F>>
// where
//     Polynomial<F>: UniqueFactorizationDomain,
//     Polynomial<F::R>: UniqueFactorizationDomain,
// {
//     todo!();
// }

// impl<R: UniqueFactorizationDomain + Hash>
//     UniqueFactorizationDomain for Polynomial<R>
// {
//     default fn factor(&self) -> Option<UniqueFactorization<Polynomial<R>>> {
//         todo!();
//     }
// }

// impl<R: UniqueFactorizationDomain + FiniteUnits + CharacteristicZero + Hash>
//     UniqueFactorizationDomain for Polynomial<R>
// {
//     default fn factor(&self) -> Option<UniqueFactorization<Polynomial<R>>> {
//         /*
//         Kronecker's method for factoring polynomials over infinite ufds with finitely many units (eg the integers)
//         Suppose we want to factor f(x) = 2 + x + x^2 + x^4 + x^5
//         Assume it has a proper factor g(x). wlog g(x) has degree <= 2
//         g(x) is determined by its value at 3 points, say at x=0, x=1, x=-1
//         f(0)=2, f(1)=6, f(-1)=2     if one of these was zero, then we would have found a linear factor
//         g(0) divides 2, g(1) divides 6, g(-1) divides 2
//         there are finitely many possible values of g(0), g(1) and g(-1) which satisfy these
//         infact there are 4*8*4=128 possible tripples
//         however, only 64 need to be checked as the other half are their negatives
//         more abstractly, some possibilities can be avoided because we only care about g up to multiplication by a unit
//          */
//         match self.degree() {
//             None => None,
//             Some(self_degree) => {
//                 println!("factor poly {}", self.to_string());
//                 let factor_degree = self_degree / 2;
//                 println!("factor degree {}", factor_degree);
//                 let mut points = vec![];
//                 let mut elem_gen = R::generate_distinct_elements();
//                 for _i in 0..factor_degree {
//                     points.push(elem_gen.next().unwrap());
//                 }
//                 println!("{:?}", points);
//                 todo!()
//             }
//         }
//     }
// }

// impl<F: FieldOfFractions + Hash> UniqueFactorizationDomain for Polynomial<F> {
//     fn factor(&self) -> Option<UniqueFactorization<Polynomial<F>>> {
//         println!("factor poly field {}", self.to_string());
//         todo!()
//     }
// }

pub fn subsylvester_matrix<R: ComRing>(
    f_deg: usize,
    g_deg: usize,
    f: &Polynomial<R>,
    g: &Polynomial<R>,
    k: usize,
) -> Matrix<R> {
    match f.degree() {
        Some(d) => assert!(d <= f_deg),
        None => {}
    }
    match g.degree() {
        Some(d) => assert!(d <= g_deg),
        None => {}
    }
    assert!(k <= f_deg);
    assert!(k <= g_deg);

    let mut smat = Matrix::zero(f_deg + g_deg - k, f_deg + g_deg - 2 * k);
    for j in 0..g_deg - k {
        for i in 0..f_deg + 1 {
            *smat.at_mut(j + i, j).unwrap() = f.coeff(f_deg - i);
        }
    }
    for i in 0..f_deg - k {
        for j in 0..g_deg + 1 {
            *smat.at_mut(j + i, g_deg - k + i).unwrap() = g.coeff(g_deg - j);
        }
    }
    smat
}

pub fn sylvester_matrix<R: ComRing>(
    f_deg: usize,
    g_deg: usize,
    f: &Polynomial<R>,
    g: &Polynomial<R>,
) -> Matrix<R> {
    subsylvester_matrix(f_deg, g_deg, f, g, 0)
}

pub fn resultant_naive<R: ComRing>(
    f_deg: usize,
    g_deg: usize,
    f: &Polynomial<R>,
    g: &Polynomial<R>,
) -> R {
    sylvester_matrix(f_deg, g_deg, f, g).det_naive().unwrap()
}

//determinant of this is the subresultant polynomial
pub fn subresultant_matrix<R: ComRing>(
    f_deg: usize,
    g_deg: usize,
    f: &Polynomial<R>,
    g: &Polynomial<R>,
    k: usize,
) -> Matrix<Polynomial<R>> {
    assert!(k <= f_deg);
    assert!(k <= g_deg);
    let tmat = subsylvester_matrix(f_deg, g_deg, f, g, k).apply_map(|a| Polynomial::from(a));
    let mut vmat = Matrix::<Polynomial<R>>::zero(f_deg + g_deg - 2 * k, f_deg + g_deg - k);
    for i in 0..f_deg + g_deg - 2 * k - 1 {
        *vmat.at_mut(i, i).unwrap() = Polynomial::one();
    }
    for i in 0..k + 1 {
        *vmat
            .at_mut(f_deg + g_deg - 2 * k - 1, f_deg + g_deg - 2 * k - 1 + i)
            .unwrap() = Polynomial::var_pow(k - i);
    }
    let prod_mat = Matrix::mul_refs(&vmat, &tmat).unwrap();
    prod_mat
}

//bad way to compute subresultant polynomials
pub fn subresultant_naive<R: ComRing>(
    f_deg: usize,
    g_deg: usize,
    f: &Polynomial<R>,
    g: &Polynomial<R>,
    k: usize,
) -> Polynomial<R> {
    subresultant_matrix(f_deg, g_deg, f, g, k)
        .det_naive()
        .unwrap()
}

#[cfg(test)]
mod tests {
    use core::panic;
    use std::collections::HashMap;

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
            Err(RingDivisionError::NotDivisible) => {}
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        let b = (2 * x + 1) * (3 * x + 2) * (4 * x + 5) * (5 * x + 6) * (6 * x + 7);
        match Polynomial::<Integer>::div(a.elem(), b.elem()) {
            Ok(_c) => panic!(),
            Err(RingDivisionError::NotDivisible) => {}
            Err(_) => panic!(),
        }

        let a = (2 * x + 1) * (3 * x + 2) * (4 * x + 5);
        let b = 0 * x;
        match Polynomial::<Integer>::div(a.elem(), b.elem()) {
            Ok(_c) => panic!(),
            Err(RingDivisionError::DivideByZero) => {}
            Err(_) => panic!(),
        }

        let a = 0 * x;
        let b = (x - x) + 5;
        match Polynomial::<Integer>::div(a.elem(), b.elem()) {
            Ok(c) => {
                assert_eq!(c, Polynomial::zero())
            }
            Err(RingDivisionError::DivideByZero) => panic!(),
            Err(_) => panic!(),
        }

        let a = 3087 * x - 8805 * x.pow(2) + 607 * x.pow(3) + x.pow(4);
        let b = (x - x) + 1;
        match Polynomial::<Integer>::div(a.elem(), b.elem()) {
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

    #[test]
    fn test_evaluate() {
        let x = &Ergonomic::new(Polynomial::<Integer>::var());
        let f = (1 + x + 3 * x.pow(2) + x.pow(3) + 7 * x.pow(4) + x.pow(5)).elem();
        assert_eq!(f.evaluate(&Integer::from(3)), Integer::from(868));

        let f = Polynomial::<Integer>::zero();
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
            let f = interpolate_by_lagrange_basis(&points).unwrap();
            for (inp, out) in &points {
                assert_eq!(&f.evaluate(&inp), out);
            }
        }

        //f(x)=2x
        match interpolate_by_lagrange_basis(&vec![
            (Integer::from(0), Integer::from(0)),
            (Integer::from(1), Integer::from(2)),
        ]) {
            Some(f) => {
                assert_eq!(f, Polynomial::new(vec![Integer::from(0), Integer::from(2)]))
            }
            None => panic!(),
        }

        //f(x)=1/2x does not have integer coefficients
        match interpolate_by_lagrange_basis(&vec![
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
            let f = interpolate_by_linear_system(&points).unwrap();
            for (inp, out) in &points {
                assert_eq!(&f.evaluate(&inp), out);
            }
        }

        //f(x)=2x
        match interpolate_by_linear_system(&vec![
            (Integer::from(0), Integer::from(0)),
            (Integer::from(1), Integer::from(2)),
        ]) {
            Some(f) => {
                assert_eq!(f, Polynomial::new(vec![Integer::from(0), Integer::from(2)]))
            }
            None => panic!(),
        }

        //f(x)=1/2x does not have integer coefficients
        match interpolate_by_linear_system(&vec![
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
        assert_eq!(
            factor_by_kroneckers_method(&f).unwrap(),
            Factored::new_unchecked(
                f.clone(),
                Polynomial::one(),
                HashMap::from([((1 + x).elem(), Natural::from(2u8))])
            )
        );

        let f = (-1 - 2 * x).elem();
        assert_eq!(
            factor_by_kroneckers_method(&f).unwrap(),
            Factored::new_unchecked(
                f.clone(),
                Polynomial::one().neg(),
                HashMap::from([((1 + 2 * x).elem(), Natural::from(1u8))])
            )
        );

        let f = (x.pow(5) + x.pow(4) + x.pow(2) + x + 2).elem();
        assert_eq!(
            factor_by_kroneckers_method(&f).unwrap(),
            Factored::new_unchecked(
                f.clone(),
                Polynomial::one(),
                HashMap::from([
                    ((1 + x + x.pow(2)).elem(), Natural::from(1u8)),
                    ((2 - x + x.pow(3)).elem(), Natural::from(1u8))
                ])
            )
        );

        let f = (1 + x + x.pow(2)).pow(2).elem();
        assert_eq!(
            factor_by_kroneckers_method(&f).unwrap(),
            Factored::new_unchecked(
                f.clone(),
                Polynomial::one(),
                HashMap::from([((1 + x + x.pow(2)).elem(), Natural::from(2u8)),])
            )
        );

        //non-primitive cases
        let f = (2 + 2 * x).elem();
        assert_eq!(
            factor_by_kroneckers_method(&f).unwrap(),
            Factored::new_unchecked(
                f.clone(),
                Polynomial::one(),
                HashMap::from([
                    (Polynomial::from_int(&Integer::from(2)), Natural::from(1u8)),
                    ((1 + x).elem(), Natural::from(1u8))
                ])
            )
        );

        let f = (12 * (2 + 3 * x) * (x - 1).pow(2)).elem();
        assert_eq!(
            factor_by_kroneckers_method(&f).unwrap(),
            Factored::new_unchecked(
                f.clone(),
                Polynomial::one(),
                HashMap::from([
                    (Polynomial::from_int(&Integer::from(2)), Natural::from(2u8)),
                    (Polynomial::from_int(&Integer::from(3)), Natural::from(1u8)),
                    ((2 + 3 * x).elem(), Natural::from(1u8)),
                    ((x - 1).elem(), Natural::from(2u8))
                ])
            )
        );

        let f = Polynomial::<Integer>::one();
        assert_eq!(
            factor_by_kroneckers_method(&f).unwrap(),
            Factored::new_unchecked(f.clone(), Polynomial::one(), HashMap::new())
        );
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
    fn test_sylvester_matrix() {
        let f = Polynomial::new(vec![
            Integer::from(1),
            Integer::from(2),
            Integer::from(3),
            Integer::from(4),
        ]);
        let g = Polynomial::new(vec![Integer::from(5), Integer::from(6), Integer::from(7)]);
        assert_eq!(
            sylvester_matrix(3, 2, &f, &g),
            Matrix::from_rows(vec![
                vec![
                    Integer::from(4),
                    Integer::from(0),
                    Integer::from(7),
                    Integer::from(0),
                    Integer::from(0)
                ],
                vec![
                    Integer::from(3),
                    Integer::from(4),
                    Integer::from(6),
                    Integer::from(7),
                    Integer::from(0)
                ],
                vec![
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(5),
                    Integer::from(6),
                    Integer::from(7)
                ],
                vec![
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(0),
                    Integer::from(5),
                    Integer::from(6)
                ],
                vec![
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(5)
                ]
            ])
        );
        assert_eq!(
            subsylvester_matrix(3, 2, &f, &g, 1),
            Matrix::from_rows(vec![
                vec![Integer::from(4), Integer::from(7), Integer::from(0)],
                vec![Integer::from(3), Integer::from(6), Integer::from(7)],
                vec![Integer::from(2), Integer::from(5), Integer::from(6)],
                vec![Integer::from(1), Integer::from(0), Integer::from(5)],
            ])
        );

        let f = Polynomial::new(vec![Integer::from(1)]);
        let g = Polynomial::new(vec![Integer::from(2)]);
        let mat = sylvester_matrix(0, 0, &f, &g);
        assert_eq!(mat, Matrix::zero(0, 0));

        let f = Polynomial::new(vec![Integer::from(1)]);
        let g = Polynomial::new(vec![Integer::from(2), Integer::from(3)]);
        let mat = sylvester_matrix(0, 1, &f, &g);
        assert_eq!(mat, Matrix::from_rows(vec![vec![Integer::from(1)],]));
    }

    #[test]
    fn test_resultant() {
        //naive algorithm (naive det of sylvester matrix)
        let f = Polynomial::new(vec![
            Integer::from(1),
            Integer::from(2),
            Integer::from(3),
            Integer::from(4),
        ]);
        let g = Polynomial::new(vec![Integer::from(5), Integer::from(6), Integer::from(7)]);
        assert_eq!(resultant_naive(3, 2, &f, &g), Integer::from(832));

        let f = Polynomial::new(vec![
            Integer::from(1),
            Integer::from(2),
            Integer::from(3),
            Integer::from(4),
        ]);
        let g = Polynomial::new(vec![Integer::from(5), Integer::from(6), Integer::from(7)]);
        assert_eq!(resultant_naive(4, 3, &f, &g), Integer::from(0));

        //subresultant
        let f = Polynomial::new(vec![
            Integer::from(1),
            Integer::from(2),
            Integer::from(3),
            Integer::from(4),
        ]);
        let g = Polynomial::new(vec![Integer::from(5), Integer::from(6), Integer::from(7)]);
        let subres = subresultant_naive(3, 2, &f, &g, 1);
        assert_eq!(
            subres,
            Polynomial::new(vec![Integer::from(64), Integer::from(-24)])
        );
    }

    // #[test]
    // fn test_squarefree_part_by_yuns() {
    //     let x = &Ergonomic::new(Polynomial::<Integer>::var());
    //     let f = ((x + 1).pow(3) * (2 * x + 3).pow(2)).elem();
    //     let g = ((x + 1) * (2 * x + 3)).elem();
    //     assert_eq!(squarefree_part_by_yuns(&f), g);
    // }
}
