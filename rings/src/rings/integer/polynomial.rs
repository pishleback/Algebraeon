use crate::{rings::integer::*, polynomial::*};
use std::rc::Rc;

impl GreatestCommonDivisorSignature for PolynomialStructure<IntegerCanonicalStructure> {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        self.gcd_by_primitive_subresultant(x.clone(), y.clone())
    }
}

impl FactorableSignature for PolynomialStructure<IntegerCanonicalStructure> {
    fn factor(&self, p: &Self::Set) -> Option<FactoredElement<Self>> {
        use berlekamp_zassenhaus::factorize_by_berlekamp_zassenhaus_algorithm;
        // self.factorize_by_kroneckers_method(p)
        factorize_by_berlekamp_zassenhaus_algorithm(p.clone())
    }
}

impl Polynomial<Integer> {
    //https://en.wikipedia.org/wiki/Landau-Mignotte_bound
    /// Return the Mignotte bound for the coefficients of any factor of $f(x)$: If
    /// $$f(x) = \sum_{i=0}^n \lambda_i x^i$$
    /// then the absolute value of the coefficients of any factor of \(f\) are bounded by
    /// $${n \choose \lfloor n/2 \rfloor} \max \lbrace |\lambda_i| : i = 0, \dots, n \rbrace$$
    pub fn mignotte_factor_coefficient_bound(&self) -> Option<Natural> {
        let deg = self.degree()?;
        let mut bound = Natural::ZERO;
        for coeff in self.coeffs() {
            bound += coeff.abs();
        }
        bound *= choose(Natural::from(deg), Natural::from(deg / 2));
        Some(bound)
    }

    //https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Lagrange's_and_Cauchy's_bounds
    pub fn cauchys_root_bound(&self) -> Option<Rational> {
        let d = self.degree()?;
        if d == 0 {
            None
        } else {
            Some(
                Rational::ONE
                    + (0..d)
                        .map(|i| Rational::from_integers(self.coeff(i), self.coeff(d)).abs())
                        .max()
                        .unwrap(),
            )
        }
    }
}

impl FactorableSignature for MultiPolynomialStructure<IntegerCanonicalStructure> {
    fn factor(&self, p: &Self::Set) -> Option<FactoredElement<Self>> {
        self.factor_by_yuns_and_kroneckers_inductively(
            Rc::new(Integer::factor),
            Rc::new(Polynomial::factor),
            p,
        )
    }
}

#[cfg(test)]
mod tests {
    use berlekamp_zassenhaus::*;

    use crate::structure::IntoErgonomic;

    use super::*;

    #[test]
    fn test_zassenhaus_against_kroneckers() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();

        let f =
            ((2 * x.pow(3) + 6 * x.pow(2) - 4) * (6 * x.pow(5) + 7 * x.pow(4) - 4)).into_verbose();
        let fs = f
            .clone()
            .factorize_by_kroneckers_method(Integer::factor)
            .unwrap();
        println!("{}", f);
        // println!("{}", f.clone().factorize_by_kroneckers_method().unwrap());
        // println!("{}", f.clone().factorize_by_zassenhaus_algorithm().unwrap());
        assert!(FactoredElement::equal(
            &fs,
            &factorize_by_berlekamp_zassenhaus_algorithm_naive(f.clone()).unwrap()
        ));
        assert!(FactoredElement::equal(
            &fs,
            &factorize_by_berlekamp_zassenhaus_algorithm(f.clone()).unwrap()
        ));

        let f = (49 * x.pow(2) - 10000).into_verbose();
        let fs = f
            .clone()
            .factorize_by_kroneckers_method(Integer::factor)
            .unwrap();
        println!("{}", f);
        // println!("{}", f.clone().factorize_by_kroneckers_method().unwrap());
        // println!("{}", f.clone().factorize_by_zassenhaus_algorithm().unwrap());
        assert!(FactoredElement::equal(
            &fs,
            &factorize_by_berlekamp_zassenhaus_algorithm_naive(f.clone()).unwrap()
        ));
        assert!(FactoredElement::equal(
            &fs,
            &factorize_by_berlekamp_zassenhaus_algorithm(f.clone()).unwrap()
        ));
    }
}
