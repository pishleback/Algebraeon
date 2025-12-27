use super::integer_lattice_ring_of_integers::RingOfIntegersWithIntegralBasisStructure;
use crate::{
    algebraic_number_field::{
        embedded::anf_multi_primitive_element_theorem, structure::AlgebraicNumberFieldSignature,
    },
    matrix::*,
    polynomial::*,
    structure::*,
};
use algebraeon_nzq::{
    Integer, Natural, Rational, RationalCanonicalStructure,
    traits::{Abs, Fraction},
};
use algebraeon_sets::structure::*;
use itertools::Itertools;
use std::borrow::{Borrow, Cow};

pub type AlgebraicNumberFieldPolynomialQuotientStructure = PolynomialQuotientRingStructure<
    RationalCanonicalStructure,
    RationalCanonicalStructure,
    PolynomialStructure<RationalCanonicalStructure, RationalCanonicalStructure>,
    true,
>;

impl Polynomial<Rational> {
    pub fn algebraic_number_field(
        self,
    ) -> Result<AlgebraicNumberFieldPolynomialQuotientStructure, ()> {
        Rational::structure()
            .into_polynomial_ring()
            .into_quotient_field(self)
    }

    pub fn algebraic_number_field_unchecked(
        self,
    ) -> AlgebraicNumberFieldPolynomialQuotientStructure {
        Rational::structure()
            .into_polynomial_ring()
            .into_quotient_field_unchecked(self)
    }

    //return the splitting field and the roots of f in the splitting field
    pub fn splitting_field(
        &self,
    ) -> (
        AlgebraicNumberFieldPolynomialQuotientStructure,
        Vec<Polynomial<Rational>>,
    ) {
        let roots = self.primitive_part_fof().all_complex_roots();
        let (g, roots_rel_g) = anf_multi_primitive_element_theorem(roots.iter().collect());
        (g.generated_algebraic_number_field(), roots_rel_g)
    }
}

impl CharZeroFieldSignature for AlgebraicNumberFieldPolynomialQuotientStructure {
    fn try_to_rat(&self, x: &Self::Set) -> Option<Rational> {
        let x = self.reduce(x);
        match x.degree() {
            None => Some(Rational::ZERO),
            Some(0) => Some(x.coeff(0).into_owned()),
            Some(_) => None,
        }
    }
}

impl<'h, B: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>>
    FreeModuleSignature<RationalCanonicalStructure>
    for RingHomomorphismRangeModuleStructure<
        'h,
        RationalCanonicalStructure,
        AlgebraicNumberFieldPolynomialQuotientStructure,
        PrincipalRationalSubfieldInclusion<AlgebraicNumberFieldPolynomialQuotientStructure, B>,
    >
{
    type Basis = EnumeratedFiniteSetStructure;

    fn basis_set(&self) -> impl std::borrow::Borrow<Self::Basis> {
        self.module()
            .coefficient_ring_inclusion()
            .range_module_structure()
            .basis_set()
            .borrow()
            .clone()
    }

    fn to_component<'a>(&self, b: &usize, v: &'a Polynomial<Rational>) -> Cow<'a, Rational> {
        self.module()
            .coefficient_ring_inclusion()
            .range_module_structure()
            .to_component(b, v)
    }

    fn from_component(&self, b: &usize, r: &Rational) -> Polynomial<Rational> {
        self.module()
            .coefficient_ring_inclusion()
            .range_module_structure()
            .from_component(b, r)
    }
}

impl AlgebraicNumberFieldSignature for AlgebraicNumberFieldPolynomialQuotientStructure {
    type Basis = EnumeratedFiniteSetStructure;
    type RingOfIntegers = RingOfIntegersWithIntegralBasisStructure;
    type RationalInclusion<B: BorrowedStructure<Self>> =
        PrincipalRationalSubfieldInclusion<Self, B>;

    fn roi(&self) -> Self::RingOfIntegers {
        self.compute_ring_of_integers()
    }

    fn finite_dimensional_rational_extension<'a>(&'a self) -> Self::RationalInclusion<&'a Self> {
        self.rational_extension()
    }
    fn into_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self> {
        self.into_rational_extension()
    }

    fn is_algebraic_integer(&self, a: &Polynomial<Rational>) -> bool {
        if self.trace(a).denominator() != Natural::ONE {
            return false;
        }
        if self.norm(a).denominator() != Natural::ONE {
            return false;
        }
        self.min_poly(a)
            .coeffs()
            .into_iter()
            .all(|c| c.denominator() == Natural::ONE)
    }
}

impl AlgebraicNumberFieldPolynomialQuotientStructure {
    pub fn compute_integral_basis_and_discriminant(&self) -> (Vec<Polynomial<Rational>>, Integer) {
        //https://www.ucl.ac.uk/~ucahmki/intbasis.pdf
        // println!("compute_basis_ring_of_integers");
        let n = self.degree();
        let mut guess = (0..n)
            .map(|i| self.integral_multiple(&Polynomial::<Rational>::var_pow(i)))
            .collect::<Vec<_>>();

        'search: loop {
            for algint in &guess {
                debug_assert!(algint.num_coeffs() <= n); //lets keep our basis alg ints reduced
            }

            let disc = self.rational_extension().discriminant(&guess);
            debug_assert_eq!(disc.clone().denominator(), Natural::ONE); //discriminant of algebraic integers is an integer
            let disc = disc.numerator();
            debug_assert_ne!(disc, Integer::ZERO); //discriminant of a basis is non-zero
            //    println!("{}", disc);
            let (_sign, mut disc_factors) = disc.factor().unwrap().into_unit_and_powers();
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
            disc_factors.sort_by_key(|(p, _k)| p.clone()); //try small primes first

            for (p, k) in disc_factors {
                debug_assert!(p >= Integer::ZERO);
                let p = p.abs().try_into().unwrap(); //if p is too big for usize then this algorithm was doomed to take longer than my lifespan anyway

                if k >= Natural::TWO {
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
                            let guess_mat =
                                Matrix::construct(n + 1, n, |r, c| guess[r].coeff(c).into_owned());
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
                                    self.from_row(
                                        guess_mat_prim_hnf
                                            .get_row_submatrix(i)
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
            return (guess, disc);
        }
    }

    pub fn compute_ring_of_integers(&self) -> RingOfIntegersWithIntegralBasisStructure {
        let (integral_basis, discriminant) = self.compute_integral_basis_and_discriminant();
        RingOfIntegersWithIntegralBasisStructure::new(self.clone(), integral_basis, discriminant)
    }
}

impl
    IntegralDomainExtensionAllPolynomialRoots<
        RationalCanonicalStructure,
        AlgebraicNumberFieldPolynomialQuotientStructure,
    >
    for PrincipalRationalSubfieldInclusion<
        AlgebraicNumberFieldPolynomialQuotientStructure,
        AlgebraicNumberFieldPolynomialQuotientStructure,
    >
{
    fn all_roots(
        &self,
        polynomial: &Polynomial<Rational>,
    ) -> Vec<<AlgebraicNumberFieldPolynomialQuotientStructure as SetSignature>::Set> {
        let anf = self.range();
        anf.polynomial_ring()
            .factor(&polynomial.apply_map(|x| self.image(x)))
            .unwrap()
            .into_powers()
            .into_iter()
            .filter_map(|(factor, power)| {
                match anf.polynomial_ring().degree(&factor) {
                    None | Some(0) => unreachable!(),
                    Some(1) => {
                        // factor = a + bx
                        // so root = -a/b
                        let a = anf.polynomial_ring().coeff(&factor, 0);
                        let b = anf.polynomial_ring().coeff(&factor, 1);
                        Some(vec![
                            anf.neg(&anf.div(a.as_ref(), b.as_ref()).unwrap());
                            power.try_into().unwrap()
                        ])
                    }
                    Some(n) => {
                        debug_assert!(n >= 2);
                        None
                    }
                }
            })
            .flatten()
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::IntoErgonomic;

    #[test]
    fn test_anf_integral_multiple() {
        let anf = Polynomial::<Rational>::from_str("200 * x^2 - 1", "x")
            .unwrap()
            .algebraic_number_field()
            .unwrap();

        println!("anf = {:?}", anf);

        let a = Polynomial::<Rational>::from_str("x", "x").unwrap();

        println!("{:?}", anf.integral_multiple(&a));
    }

    #[test]
    fn test_anf_to_and_from_vector() {
        let x = &Polynomial::<Rational>::var().into_ergonomic();
        let anf = (x.pow(5) - x + 1)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let alpha = (x.pow(9) + 5).into_verbose();

        println!("{}", alpha);
        println!("{}", anf.reduce(&alpha));
        println!("{}", anf.min_poly(&alpha));
        anf.to_col(&alpha).pprint();

        assert_eq!(
            anf.to_col(&alpha),
            Matrix::from_cols(vec![vec![
                Rational::from(4),
                Rational::from(1),
                Rational::from(0),
                Rational::from(0),
                Rational::from(-1)
            ]])
        );

        assert!(anf.equal(
            &anf.from_col(Matrix::from_cols(vec![vec![
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
