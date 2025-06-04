use algebraeon_nzq::{Integer, Natural, Rational};

use crate::parsing::polynomial::{parse_integer_polynomial, parse_rational_polynomial};

#[derive(Debug, Clone)]
pub struct Polynomial<Set> {
    //vec![c0, c1, c2, c3, ..., cn] represents the polynomial c0 + c1*x + c2*x^2 + c3*x^3 + ... + cn * x^n
    //if non-empty, the last item must not be zero
    pub coeffs: Vec<Set>,
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
        #[allow(clippy::redundant_closure_for_method_calls)]
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

    pub fn apply_map_into<ImgSet>(self, f: impl Fn(Set) -> ImgSet) -> Polynomial<ImgSet> {
        Polynomial::from_coeffs(self.coeffs.into_iter().map(f).collect())
    }

    pub fn apply_map_with_powers<ImgSet>(
        &self,
        f: impl Fn((usize, &Set)) -> ImgSet,
    ) -> Polynomial<ImgSet> {
        Polynomial::from_coeffs(self.coeffs.iter().enumerate().map(f).collect())
    }

    pub fn apply_map_into_with_powers<ImgSet>(
        self,
        f: impl Fn((usize, Set)) -> ImgSet,
    ) -> Polynomial<ImgSet> {
        Polynomial::from_coeffs(self.coeffs.into_iter().enumerate().map(f).collect())
    }
}

pub trait PolynomialFromStr: Sized {
    type Err;

    fn from_str(polynomial_str: &str, var: &str) -> Result<Self, Self::Err>;
}

impl PolynomialFromStr for Polynomial<Natural> {
    type Err = ();

    fn from_str(polynomial_str: &str, var: &str) -> Result<Self, ()> {
        Ok(Self::from_coeffs(
            Polynomial::<Integer>::from_str(polynomial_str, var)?
                .into_coeffs()
                .into_iter()
                .map(Natural::try_from)
                .collect::<Result<Vec<_>, _>>()?,
        ))
    }
}

impl PolynomialFromStr for Polynomial<Integer> {
    type Err = ();

    fn from_str(polynomial_str: &str, var: &str) -> Result<Self, ()> {
        parse_integer_polynomial(polynomial_str, var).map_err(|_| ())
    }
}

impl PolynomialFromStr for Polynomial<Rational> {
    type Err = ();

    fn from_str(polynomial_str: &str, var: &str) -> Result<Self, ()> {
        parse_rational_polynomial(polynomial_str, var).map_err(|_| ())
    }
}
