use crate::{integer::*, polynomial::*};
use std::rc::Rc;

impl<B: BorrowedStructure<IntegerCanonicalStructure>> GreatestCommonDivisorSignature
    for PolynomialStructure<IntegerCanonicalStructure, B>
{
    fn gcd(&self, x: &Self::Set, y: &Self::Set) -> Self::Set {
        self.gcd_by_primitive_subresultant(x.clone(), y.clone())
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> FactoringMonoidSignature
    for PolynomialStructure<IntegerCanonicalStructure, B>
{
    fn factor_unchecked(&self, p: &Self::Set) -> Factored<Polynomial<Integer>, Natural> {
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

// #[derive(Debug, Clone, PartialEq, Eq)]
// struct IrreducibleIntegerMultiPolynomialStructure {}

// impl Signature for IrreducibleIntegerMultiPolynomialStructure {}

// impl SetSignature for IrreducibleIntegerMultiPolynomialStructure {
//     type Set = MultiPolynomial<Integer>;

//     fn is_element(&self, x: &Self::Set) -> bool {
//         todo!()
//     }
// }

// impl EqSignature for IrreducibleIntegerMultiPolynomialStructure {
//     fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
//         Integer::structure().multivariable_polynomials().equal(a, b)
//     }
// }

// impl OrdSignature for IrreducibleIntegerMultiPolynomialStructure {
//     fn cmp(&self, a: &Self::Set, b: &Self::Set) -> Ordering {
//         todo!()
//     }
// }

// impl<B: BorrowedStructure<IntegerCanonicalStructure> + 'static> UniqueFactorizationSignature
//     for MultiPolynomialStructure<IntegerCanonicalStructure, B>
// {
//     type Irreducibles = IrreducibleIntegerMultiPolynomialStructure;

//     type Factorizations<SelfB: BorrowedStructure<Self>> = FactoredRingElementStructure<Self, SelfB>;

//     fn factorizations<'a>(&'a self) -> Self::Factorizations<&'a Self> {
//         FactoredRingElementStructure::new(self)
//     }

//     fn into_factorizations(self) -> Self::Factorizations<Self> {
//         FactoredRingElementStructure::new(self)
//     }

//     fn irreducibles(&self) -> impl std::borrow::Borrow<Self::Irreducibles> {
//         IrreducibleIntegerMultiPolynomialStructure {}
//     }
// }

// impl<B: BorrowedStructure<IntegerCanonicalStructure> + 'static> UniqueFactorizationSignature
//     for MultiPolynomialStructure<PolynomialStructure<IntegerCanonicalStructure, B>, B>
// {

// }

impl<B: BorrowedStructure<IntegerCanonicalStructure> + 'static> FactoringMonoidSignature
    for MultiPolynomialStructure<IntegerCanonicalStructure, B>
{
    fn factor_unchecked(&self, p: &Self::Set) -> Factored<Self::Set, Natural> {
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
        assert!(
            Polynomial::<Integer>::structure()
                .factorizations()
                .equal(&fs, &factorize_by_berlekamp_zassenhaus_algorithm(f.clone()))
        );
    }
}
