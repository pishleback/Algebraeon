use itertools::Itertools;
use malachite_base::num::basic::traits::{One, Zero};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

use crate::rings::{
    linear::matrix::Matrix,
    number::natural::nat_to_usize,
    polynomial::polynomial::*,
    ring_structure::{
        cannonical::{
            FieldOfFractionsDomain, GreatestCommonDivisorDomain, Ring, UniqueFactorizationDomain,
        },
        quotient::*,
        structure::RingStructure,
    },
    structure::*,
};

pub type ANFStructure = QuotientStructure<PolynomialStructure<CannonicalStructure<Rational>>, true>;

pub fn new_anf(f: Polynomial<Rational>) -> ANFStructure {
    ANFStructure::new(PolynomialStructure::new(Rational::structure()).into(), f)
}

impl ANFStructure {
    pub fn degree(&self) -> usize {
        self.modulus().degree().unwrap()
    }

    pub fn trace_form_matrix(&self, elems: &Vec<Polynomial<Rational>>) -> Matrix<Rational> {
        let n = self.degree();
        assert_eq!(n, elems.len());
        Matrix::construct(n, n, |r, c| {
            self.trace(&Polynomial::mul(&elems[r], &elems[c]))
        })
    }

    pub fn discriminant(&self, elems: &Vec<Polynomial<Rational>>) -> Rational {
        self.trace_form_matrix(elems).det().unwrap()
    }

    pub fn compute_integral_basis(&self) -> Vec<Polynomial<Rational>> {
        //https://www.ucl.ac.uk/~ucahmki/intbasis.pdf
        println!("compute_basis_ring_of_integers");
        let n = self.degree();
        let mut guess = (0..n)
            .map(|i| self.integral_multiple(&Polynomial::<Rational>::var_pow(i)))
            .collect_vec();

        'search: loop {
            for algint in &guess {
                debug_assert!(algint.num_coeffs() <= n); //lets keep our basis alg ints reduced
            }

            let disc = self.discriminant(&guess);
            debug_assert_eq!(disc.denominator_ref(), &Natural::ONE); //discriminant of algebraic integers is an integer
            let disc = Rational::numerator(&disc);
            debug_assert_ne!(disc, Integer::ZERO); //discriminant of a basis is non-zero
            // println!("{}", disc);
            let (_sign, disc_factors) = disc.factor().unwrap().unit_and_factors();
            // If p is a prime such that p^2 divides Disc
            // then can find an alg int of the form
            // 1/p (x_1a_1 + ... + x_na_n)
            // 0 <= x_i <= p - 1 and x_i in Z
            // where {a_i} is the current guess at an integral basis
            // If no algebraic integers of this form exist then we have an actual integral basis
            // If one does exist then we can add it to the integral basis & reduce to get a new guess at a basis

            // println!("guess = {:?}", guess);
            // println!("disc = {:?}", disc);
            // println!("disc_factors = {:?}", disc_factors);
            for (p, k) in disc_factors {
                debug_assert!(p >= Integer::ZERO);
                let p = nat_to_usize(p.unsigned_abs_ref()).unwrap(); //if p is too big for usize then this algorithm was doomed to take longer than my lifespan anyway

                if k >= 2 {
                    // println!("p = {}", p);

                    for coeffs in (0..n).map(|_i| 0..p).multi_cartesian_product() {
                        let alpha = Polynomial::from_coeffs(
                            Polynomial::sum(
                                (0..n)
                                    .map(|i| {
                                        Polynomial::mul(
                                            &Polynomial::constant(Rational::from(coeffs[i])),
                                            &guess[i],
                                        )
                                    })
                                    .collect(),
                            )
                            .into_coeffs()
                            .into_iter()
                            .map(|c| c / Rational::from(p))
                            .collect(),
                        );

                        // println!("coeffs = {:?}  alpha = {:?}  min_poly = {}", coeffs, alpha, self.min_poly(&alpha));

                        if !self.is_zero(&alpha) && self.is_algebraic_integer(&alpha) {
                            // println!("alpha = {:?} {}", alpha, self.min_poly(&alpha));

                            guess.push(alpha);
                            let guess_mat = Matrix::construct(n + 1, n, |r, c| guess[r].coeff(c));
                            let (mul, guess_mat_prim) = guess_mat.factor_primitive_fof();
                            let guess_mat_prim_hnf = guess_mat_prim
                                .flip_cols()
                                .row_reduced_hermite_normal_form()
                                .flip_cols();

                            // guess_mat.pprint();
                            // guess_mat_prim_hnf.pprint();

                            // println!("{:?}", mul);

                            guess = (0..n)
                                .rev()
                                .map(|i| {
                                    self.from_row_vector(
                                        guess_mat_prim_hnf
                                            .get_row(i)
                                            .apply_map(|v| Rational::from(v) * &mul),
                                    )
                                })
                                .collect();

                            // println!("new_guess = {:?}", guess);
                            continue 'search;
                        }
                    }
                }
            }
            return guess;
        }
    }

    //matrix representing column vector multiplication by a on the left
    pub fn col_multiplication_matrix(&self, a: &Polynomial<Rational>) -> Matrix<Rational> {
        let a_reduced = self.reduce(a);
        let deg = self.degree();
        Matrix::from_cols(
            (0..self.degree())
                .map(|i| {
                    let mut coeffs = self
                        .reduce(&Polynomial::mul(a, &Polynomial::var_pow(i)))
                        .into_coeffs();
                    debug_assert!(coeffs.len() <= deg);
                    while coeffs.len() < deg {
                        coeffs.push(Rational::zero());
                    }
                    coeffs
                })
                .collect(),
        )
    }
    //matrix representing row vector multiplication by a on the right
    pub fn row_multiplication_matrix(&self, a: &Polynomial<Rational>) -> Matrix<Rational> {
        self.col_multiplication_matrix(a).transpose()
    }

    pub fn to_col_vector(&self, a: &Polynomial<Rational>) -> Matrix<Rational> {
        let a_reduced = self.reduce(a);
        Matrix::construct(self.degree(), 1, |r, c| a_reduced.coeff(r))
    }
    pub fn to_row_vector(&self, a: &Polynomial<Rational>) -> Matrix<Rational> {
        self.to_col_vector(a).transpose()
    }

    pub fn from_col_vector(&self, v: Matrix<Rational>) -> Polynomial<Rational> {
        assert_eq!(v.cols(), 1);
        assert_eq!(v.rows(), self.degree());
        Polynomial::from_coeffs(
            (0..self.degree())
                .map(|i| v.at(i, 0).unwrap().clone())
                .collect(),
        )
    }
    pub fn from_row_vector(&self, v: Matrix<Rational>) -> Polynomial<Rational> {
        self.from_col_vector(v.transpose())
    }

    pub fn min_poly(&self, a: &Polynomial<Rational>) -> Polynomial<Rational> {
        self.col_multiplication_matrix(a)
            .minimal_polynomial()
            .unwrap()
    }

    pub fn norm(&self, a: &Polynomial<Rational>) -> Rational {
        self.col_multiplication_matrix(a).det().unwrap()
    }

    pub fn trace(&self, a: &Polynomial<Rational>) -> Rational {
        self.col_multiplication_matrix(a).trace().unwrap()
    }

    pub fn is_algebraic_integer(&self, a: &Polynomial<Rational>) -> bool {
        if self.trace(a).denominator_ref() != &Natural::ONE {
            return false;
        }
        if self.norm(a).denominator_ref() != &Natural::ONE {
            return false;
        }
        self.min_poly(a)
            .coeffs()
            .into_iter()
            .all(|c| c.denominator_ref() == &Natural::ONE)
    }

    //return a scalar multiple of a which is an algebraic integer
    fn integral_multiple(&self, a: &Polynomial<Rational>) -> Polynomial<Rational> {
        let m = Integer::lcm_list(
            self.min_poly(a)
                .coeffs()
                .into_iter()
                .map(|c| Integer::from(c.denominator_ref()))
                .collect(),
        );
        let b = Polynomial::mul(&Polynomial::constant(Rational::from(m)), a);
        debug_assert!(self.is_algebraic_integer(&b));
        b
    }
}

struct RingOfIntegers {
    anf: ANFStructure,
    basis: Vec<Polynomial<Rational>>,
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_anf_to_and_from_vector() {
        let x = &Polynomial::<Rational>::var().into_ring();
        let anf = new_anf((x.pow(5) - x + 1).into_set());
        let alpha = (x.pow(9) + 5).into_set();

        println!("{}", alpha);
        println!("{}", anf.reduce(&alpha));
        println!("{}", anf.min_poly(&alpha));
        anf.to_col_vector(&alpha).pprint();

        assert_eq!(
            anf.to_col_vector(&alpha),
            Matrix::from_cols(vec![vec![
                Rational::from(4),
                Rational::from(1),
                Rational::from(0),
                Rational::from(0),
                Rational::from(-1)
            ]])
        );

        assert!(anf.equal(
            &anf.from_col_vector(Matrix::from_cols(vec![vec![
                Rational::from(4),
                Rational::from(1),
                Rational::from(0),
                Rational::from(0),
                Rational::from(-1)
            ]])),
            &alpha
        ));
    }
}
