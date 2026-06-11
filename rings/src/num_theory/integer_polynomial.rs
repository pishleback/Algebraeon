use crate::{
    num_theory::berlekamp_zassenhaus::factorize_by_berlekamp_zassenhaus_algorithm,
    polynomial::*,
    structure::{
        Factored, FactoringMonoidSignature, GreatestCommonDivisorSignature, MetaFactoringMonoid,
    },
};
use algebraeon_structures::*;
use std::rc::Rc;

impl<B: BorrowedStructure<IntegerCanonicalStructure>> GreatestCommonDivisorSignature
    for PolynomialStructure<IntegerCanonicalStructure, B>
{
    fn gcd(&self, x: &Self::Elem, y: &Self::Elem) -> Self::Elem {
        self.gcd_by_primitive_subresultant(x.clone(), y.clone())
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> FactoringMonoidSignature
    for PolynomialStructure<IntegerCanonicalStructure, B>
{
    fn factor_unchecked(&self, p: &Self::Elem) -> Factored<Polynomial<Integer>, Natural> {
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
            bound += Abs::abs(coeff);
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
                        .map(|i| {
                            Rational::from_integers(self.coeff(i).as_ref(), self.coeff(d).as_ref())
                                .abs()
                        })
                        .max()
                        .unwrap(),
            )
        }
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure> + 'static> FactoringMonoidSignature
    for MultiPolynomialStructure<IntegerCanonicalStructure, B>
{
    fn factor_unchecked(&self, p: &Self::Elem) -> Factored<Self::Elem, Natural> {
        self.factor_by_yuns_and_kroneckers_inductively(
            Rc::new(Integer::factor),
            Rc::new(Polynomial::factor),
            p,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        num_theory::berlekamp_zassenhaus::factorize_by_berlekamp_zassenhaus_algorithm_naive,
        structure::{IntoErgonomic, UniqueFactorizationMonoidSignature},
    };

    #[test]
    fn test_zassenhaus_against_kroneckers() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();

        let f =
            ((2 * x.pow(3) + 6 * x.pow(2) - 4) * (6 * x.pow(5) + 7 * x.pow(4) - 4)).into_verbose();
        let fs = Integer::structure()
            .polynomials()
            .factorize_by_kroneckers_method(f.clone(), Integer::factor);
        println!("{}", f);
        // println!("{}", f.clone().factorize_by_kroneckers_method().unwrap());
        // println!("{}", f.clone().factorize_by_zassenhaus_algorithm().unwrap());
        assert!(Polynomial::<Integer>::structure().factorizations().equal(
            &fs,
            &factorize_by_berlekamp_zassenhaus_algorithm_naive(f.clone())
        ));
        assert!(
            Polynomial::<Integer>::structure()
                .factorizations()
                .equal(&fs, &factorize_by_berlekamp_zassenhaus_algorithm(f.clone()))
        );

        let f = (49 * x.pow(2) - 10000).into_verbose();
        let fs = Integer::structure()
            .polynomials()
            .factorize_by_kroneckers_method(f.clone(), Integer::factor);
        println!("{}", f);
        // println!("{}", f.clone().factorize_by_kroneckers_method().unwrap());
        // println!("{}", f.clone().factorize_by_zassenhaus_algorithm().unwrap());
        assert!(Polynomial::<Integer>::structure().factorizations().equal(
            &fs,
            &factorize_by_berlekamp_zassenhaus_algorithm_naive(f.clone())
        ));

        println!("{}", fs);
        println!("{}", factorize_by_berlekamp_zassenhaus_algorithm(f.clone()));

        assert!(
            Polynomial::<Integer>::structure()
                .factorizations()
                .equal(&fs, &factorize_by_berlekamp_zassenhaus_algorithm(f.clone()))
        );
    }
}
