use crate::{number::integer::*, number::natural::functions::*, polynomial::polynomial::*};

impl GreatestCommonDivisorStructure for PolynomialStructure<CannonicalStructure<Integer>> {
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        self.gcd_by_primitive_subresultant(x.clone(), y.clone())
    }
}

impl UniqueFactorizationStructure for PolynomialStructure<CannonicalStructure<Integer>> {
    fn factor(&self, p: &Self::Set) -> Option<Factored<Self>> {
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
        bound *= choose_usize(deg, deg / 2);
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

#[cfg(test)]
mod tests {
    use berlekamp_zassenhaus::*;

    use crate::structure::elements::IntoErgonomic;

    use super::*;

    #[test]
    fn test_zassenhaus_against_kroneckers() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();

        let f =
            ((2 * x.pow(3) + 6 * x.pow(2) - 4) * (6 * x.pow(5) + 7 * x.pow(4) - 4)).into_verbose();
        let fs = f.clone().factorize_by_kroneckers_method().unwrap();
        println!("{}", f);
        // println!("{}", f.clone().factorize_by_kroneckers_method().unwrap());
        // println!("{}", f.clone().factorize_by_zassenhaus_algorithm().unwrap());
        assert!(Factored::equal(
            &fs,
            &factorize_by_berlekamp_zassenhaus_algorithm_naive(f.clone()).unwrap()
        ));
        assert!(Factored::equal(
            &fs,
            &factorize_by_berlekamp_zassenhaus_algorithm(f.clone()).unwrap()
        ));

        let f = (49 * x.pow(2) - 10000).into_verbose();
        let fs = f.clone().factorize_by_kroneckers_method().unwrap();
        println!("{}", f);
        // println!("{}", f.clone().factorize_by_kroneckers_method().unwrap());
        // println!("{}", f.clone().factorize_by_zassenhaus_algorithm().unwrap());
        assert!(Factored::equal(
            &fs,
            &factorize_by_berlekamp_zassenhaus_algorithm_naive(f.clone()).unwrap()
        ));
        assert!(Factored::equal(
            &fs,
            &factorize_by_berlekamp_zassenhaus_algorithm(f.clone()).unwrap()
        ));
    }
}
