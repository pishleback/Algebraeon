// https://en.wikipedia.org/wiki/Berlekamp%27s_algorithm
/*!
 * # The Berlekamp Zassenhaus algorithm
 * The Berlekamp Zassenhaus algorithm (BZA) is a method for factoring integer polynomials.
 * BZA works on primitive squarefree polynomials, so other algorithms must be used first to reduce the factorization problem to this case.
 *
 * Let $f(x) \in \mathbb{Z}\[x\]$ be non-zero primitive squarefree polynomial of degree at least $1$
 * 1. A bound $B \in \mathbb{N}$ is found such that the absolute value of the coefficients of any factor of $f$ are bounded above by $B$.
 * 2. A prime number $p$ is found such that the reduction $\tilde{f}(x) \in \mathbb{F}_p(x)$ of $f$ modulo $p$ is has the same degree as $f$ and that $\tilde{f}$ is squarefree.
 * 3. $\tilde{f}(x)$ is factored using methods for polynomial factorization over finite fields.
 *    Let $k$ be the number of factors, and let $\tilde{g}_1(x), \dots, \tilde{g}_k(x)$ denote the factors, so
 *    $$\tilde{f}(x) = \tilde{g}_1(x) \cdots \tilde{g}_k(x)$$
 *    The polynomials $\tilde{g}_1(x), \dots, \tilde{g}_k(x)$ are called the modular factors.
 * 4. The factorization $\tilde{f}(x) = \tilde{g}_1(x) \cdots \tilde{g}_k(x)$ of $f$ modulo $p$ is lifted to factorizations modulo $p^t$ for sufficiently large $t$ such that $p^t \ge 2B$.
 * 5. Any factorization of $f$ over $\mathbb{Z}$ would yield a factorization modulo $p^t$, so to factor $f(x)$ it suffices to check all products of subsets of the modular factors and check if they yield a true factor of $f$ when the coefficients are lifted to $\mathbb{Z}$.
 *    The inequality $p^t \ge 2B$ ensures that the lifting to $\mathbb{Z}$ is unique, so all factors of $f$ will be found in this way.
 *
 * # Implementation and Optimizations
 * Step 5 of the algorithm is the bottleneck, since the time taken to check all subsets of the modular factors grows exponentially with the number of modular factors.
 * Despite this exponential growth, even a naive implementation of the algorithm performs well.
 *
 * There are also ways to speed up the search process in step 5.
 *  - Complimentary subsets: Only check half the subsets of the modular factors by ignoring one from each complimentary pair.
 *  - The d-1 test: (TODO)
 *  - The d-2 test: (TODO)
 *  - Degree sets: (TODO)
 *  - Memory stacks: (TODO)
 *  - LLL basis reduction methods: (TODO)
 *
 * The naive implementation of BZA loops over half the subsets of the modular factors (excluding complimentary pairs) and sees if they produce a factor of $f$ by performing a polynomial division.
 *
 */

/*
Justifying 2B <= p^t:
   Let B be the root bound of f(x).
   For each modular coefficient c we can make
       0 <= c < p^t
       or
       -floor(p^t/2) <= c <= floor(p^t/2)
   For each coefficient c of a factor of f(x)
       -B <= c <= B
   So for unique lifting we need
       B <= floor(p^t/2)
   If t is such that 2B <= p^t then this is the case since
       2B <= p^t
    -> B <= p^t/2
    -> B <= floor(p^t/2) (since B is an integer)
*/

/*
some improvements
 - knapsack improvement https://www.math.fsu.edu/~hoeij/knapsack/paper/knapsack.pdf
   This might be of use when I do the LLL bit, especially 10 and 11: http://people.csail.mit.edu/madhu/FT98/
 - brute search optimizations https://www.shoup.net/papers/asz.pdf
 - berlekamp_zassenhaus_algorithm
 - don't use factor by find factor, just do it all in one go and partition the modular factors into true factors

*/

use crate::number::integer::*;
use crate::number::natural::*;
use crate::number::rational::*;
use crate::{polynomial::polynomial::*, structure::quotient::*};
use algebraeon_sets::combinations::LexicographicCombinationsWithRemovals;
use primes::PrimeGenerator;

fn compute_polynomial_factor_bound(poly: &Polynomial<Integer>) -> Natural {
    poly.mignotte_factor_coefficient_bound().unwrap()
}

#[derive(Debug)]
struct BerlekampAassenhausAlgorithmState {
    poly: Polynomial<Integer>,
    degree: usize,
    factor_coeff_bound: Natural,
    minimum_modolus: Natural,
    prime_gen: PrimeGenerator,
}

impl BerlekampAassenhausAlgorithmState {
    fn new(poly: Polynomial<Integer>) -> Self {
        let degree = poly.degree().unwrap();
        debug_assert_ne!(degree, 0);
        let factor_coeff_bound = compute_polynomial_factor_bound(&poly);
        let minimum_modolus = Natural::TWO * &factor_coeff_bound;
        let prime_gen = PrimeGenerator::new();
        Self {
            poly,
            degree,
            factor_coeff_bound,
            minimum_modolus,
            prime_gen,
        }
    }

    fn next_prime(&mut self) -> BerlekampZassenhausAlgorithmStateAtPrime {
        loop {
            let p = self.prime_gen.next().unwrap();
            if let Some(state_at_p) =
                BerlekampZassenhausAlgorithmStateAtPrime::new_at_prime(self, p)
            {
                return state_at_p;
            }
        }
    }
}

#[derive(Debug)]
struct BerlekampZassenhausAlgorithmStateAtPrime {
    poly: Polynomial<Integer>,
    leading_coeff: Integer,
    degree: usize,
    modulus: Integer,
    modular_factors: Vec<Polynomial<Integer>>,
}

impl BerlekampZassenhausAlgorithmStateAtPrime {
    fn new_at_prime(state: &BerlekampAassenhausAlgorithmState, p: Natural) -> Option<Self> {
        let mod_p = QuotientStructure::new_field(Integer::structure(), Integer::from(&p));
        let poly_mod_p = PolynomialStructure::new(mod_p.into());
        if poly_mod_p.degree(&state.poly) == Some(state.degree) {
            let facotred_f_mod_p = poly_mod_p.factor(&state.poly).unwrap();
            match facotred_f_mod_p.into_hensel_factorization(state.poly.clone()) {
                Some(hensel_factorization_f_over_p) => {
                    let mut hensel_factorization_f_over_p =
                        hensel_factorization_f_over_p.dont_lift_bezout_coeffs();
                    while hensel_factorization_f_over_p.modolus() < state.minimum_modolus {
                        hensel_factorization_f_over_p.linear_lift();
                    }
                    let modulus = hensel_factorization_f_over_p.modolus();
                    let modular_factors = hensel_factorization_f_over_p
                        .factors()
                        .into_iter()
                        .map(|g| g.clone())
                        .collect();
                    Some(BerlekampZassenhausAlgorithmStateAtPrime {
                        poly: state.poly.clone(),
                        leading_coeff: state.poly.leading_coeff().unwrap(),
                        degree: state.degree.clone(),
                        modulus,
                        modular_factors,
                    })
                }
                None => None,
            }
        } else {
            None
        }
    }
}

trait SemigroupStructure: Structure {
    fn compose(&self, a: &Self::Set, b: &Self::Set) -> Self::Set;
}

struct MemoryStack<SG: SemigroupStructure> {
    semigroup: SG,
    modular_factor_values: Vec<SG::Set>,
    // Store the partial products of a previous calculation
    // Since subsets are visited in lexographic order, if a test is performed frequently, values towards the right will need to be updated
    /*
    For example if this vector contains
        [(0, a), (2, b), (3, c)]
    That means the subset [0, 2, 3] has partial products
        v(0) = a
        v(0)v(2) = b
        v(0)v(2)v(3) = c
     */
    prev_calc: Vec<(usize, SG::Set)>,
}

impl<SG: SemigroupStructure> MemoryStack<SG> {
    fn new(semigroup: SG, modular_factor_values: Vec<SG::Set>) -> Self {
        Self {
            semigroup,
            modular_factor_values,
            prev_calc: vec![],
        }
    }

    fn get_val(&self, i: usize) -> &<SG as Structure>::Set {
        &self.modular_factor_values[i]
    }

    fn get_product(&mut self, subset: &Vec<usize>) -> &<SG as Structure>::Set {
        debug_assert!(!subset.is_empty());
        let mut i = 0;
        loop {
            if i == self.prev_calc.len() {
                break;
            }
            if i == subset.len() {
                break;
            }
            if self.prev_calc[i].0 != subset[i] {
                break;
            }
            i += 1;
        }
        /*
        By way of example, suppose
            subset = abcxyz
            memory = a ab acd abcd abcde
        Then we use abc from memory and compose with the remaining items xyz
         */

        // subset and memory agree in the first i places
        self.prev_calc.truncate(i);
        if i == 0 {
            let ss = subset[0];
            self.prev_calc.push((ss, self.get_val(ss).clone()));
            i += 1;
        }
        // subset and memory still agree in the first i >=1 places
        for j in i..subset.len() {
            let ss = subset[j];
            self.prev_calc.push((
                ss,
                self.semigroup
                    .compose(&self.prev_calc[j - 1].1, self.get_val(ss)),
            ));
        }

        &self.prev_calc[subset.len() - 1].1
    }
}

/// The (d-1) test
/// A quick test allowing some subsets of modular factors to be ruled out from yielding true factors.
///
/// A bound is found for the (d-1)st coefficient of any factors of $f(x)$.
/// The (d-1)st coefficient of the product of modular factors is easily computed as the sum of the (d-1)st coefficients of each modular factor.
/// The (d-1)st coefficients of the modular factors can be quickly summed and checked whether they are in the possible range for true factors.
/// Even better, this summing can be translated to machine arithmetic taking advantage of the wrapping behavour of binary addition. This is at the cost of a slight lossening of the bound on the (d-1)st coefficient.
mod dminusone_test {
    use super::*;

    #[derive(Debug, Clone)]
    struct DMinusOneTestSemigroupElem {
        approx_coeff_lower_bound: usize, // An upper bound is given by this +1
        degree: usize,
    }
    #[derive(Debug, Clone, PartialEq, Eq)]
    struct DMinusOneTestSemigroup {}
    impl Structure for DMinusOneTestSemigroup {
        type Set = DMinusOneTestSemigroupElem;
    }
    impl SemigroupStructure for DMinusOneTestSemigroup {
        fn compose(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
            DMinusOneTestSemigroupElem {
                approx_coeff_lower_bound: a
                    .approx_coeff_lower_bound
                    .wrapping_add(b.approx_coeff_lower_bound),
                degree: a.degree.wrapping_add(b.degree),
            }
        }
    }

    pub struct DMinusOneTest {
        memory_stack: MemoryStack<DMinusOneTestSemigroup>,
        factor_dminusone_coeff_bound_divby_gdeg: usize,
    }
    impl DMinusOneTest {
        pub fn new(
            modulus: &Integer,
            f: &Polynomial<Integer>,
            modular_factors: &Vec<Polynomial<Integer>>,
        ) -> Self {
            use malachite_base::num::arithmetic::traits::{Abs, Mod};

            // Probably 2^64
            let machine_range = Natural::from(usize::MAX) + Natural::ONE;

            // Conversion multiplier from 0,...,p^t to 0,...,2^64
            let conversion_mult =
                Rational::from_integers(Integer::from(&machine_range), modulus.clone());

            /*
            Compute a bound for the (d-1)st coefficient of any factor of f using a root bound.

            If g(x) = M(x-a)(x-b)(x-c)(x-d) is a factor of f(x) then it has the same roots and g(x) = Mx^4 - M(a + b + c + d)x^3 + ...
            The (d-1)st = 3rd coefficient of g(x) is M(a + b + c + d)
            So if B is a root bound for f(x) then the (d-1)st coefficient of g(x) is bounded by M*B*deg(g)

            Here we compute M*B and leave the multiplication by deg(g) for later
            */
            let factor_dminusone_coeff_bound_divby_gdeg =
                Rational::from(f.leading_coeff().unwrap().abs()) * f.cauchys_root_bound().unwrap();

            Self {
                factor_dminusone_coeff_bound_divby_gdeg: nat_to_usize(
                    &(&factor_dminusone_coeff_bound_divby_gdeg * &conversion_mult)
                        .ceil()
                        .unsigned_abs(),
                )
                .unwrap_or(usize::MAX),
                memory_stack: MemoryStack::new(
                    DMinusOneTestSemigroup {},
                    modular_factors
                        .iter()
                        .map(|g| {
                            let d = g.degree().unwrap();
                            let coeff = (f.leading_coeff().unwrap() * g.coeff(d - 1)).rem(modulus);
                            DMinusOneTestSemigroupElem {
                                approx_coeff_lower_bound: nat_to_usize(
                                    &(Rational::from(coeff) * &conversion_mult)
                                        .floor()
                                        .unsigned_abs()
                                        .rem(&machine_range),
                                )
                                .unwrap(),
                                degree: d,
                            }
                        })
                        .collect(),
                ),
            }
        }

        /// Return true if the subset definitely wont yield a true factor.
        /// Return false if the subset might yield a true factor.
        pub fn test(&mut self, subset: &Vec<usize>) -> bool {
            let DMinusOneTestSemigroupElem {
                approx_coeff_lower_bound: range_bot,
                degree: d,
            } = *self.memory_stack.get_product(subset);
            let range_top = range_bot.wrapping_add(subset.len());
            /*
            In the machine coordinates 0-usize::MAX we have:
                The (d-1) coefficient of the lifted product lies in the range [range_bot, range_top]
                The (d-1) coefficient of any true factor is bounded above by b := r*d where r is a root bound for f
            So we need to check whether [0, b] union [usize::MAX - b + 1, usize::MAX] is disjoint from [range_bot, range_top]
            */
            let b1 = self
                .factor_dminusone_coeff_bound_divby_gdeg
                .checked_mul(d)
                .unwrap_or(usize::MAX);
            if b1 == 0 {
                0 < range_bot && range_bot <= range_top
            } else {
                let b2 = usize::MAX - b1 + 1;
                b1 < range_bot && range_bot <= range_top && range_top < b2
            }
        }
    }
}

// Polynomial division test. This test is never wrong.
type ModularFactorMultSemigrp =
    PolynomialStructure<QuotientStructure<CannonicalStructure<Integer>, false>>;
impl SemigroupStructure for ModularFactorMultSemigrp {
    fn compose(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.mul(a, b)
    }
}

impl BerlekampZassenhausAlgorithmStateAtPrime {
    fn factor_by_try_all_subsets<'a>(
        &'a self,
    ) -> Factored<PolynomialStructure<CannonicalStructure<Integer>>> {
        let n = self.modular_factors.len();

        let mut dminusone_test =
            dminusone_test::DMinusOneTest::new(&self.modulus, &self.poly, &self.modular_factors);
        let mut modular_factor_product_memory_stack = MemoryStack::new(
            PolynomialStructure::new(
                QuotientStructure::new_ring(Integer::structure(), self.modulus.clone()).into(),
            ),
            self.modular_factors.iter().map(|g| g.clone()).collect(),
        );

        let mut factored = Factored::factored_one(Polynomial::<Integer>::structure());
        let mut excluded_modular_factors = vec![];
        let mut f = self.poly.clone();
        let mut k = 1; // The cardinality of the subset to search for each loop
        let mut m = n; // Keep track of the number of remaining modular factors i.e. m = n - excluded_modular_factors.len()

        // Loop over subsets of modular factors in order of increasing cardinality, then in lexographic order for each cardinality
        // Only half of the cardinalities need to be checked since complimentary subsets need not be checked
        // Keep searching for cardinalities up to and including half the number of remaining modular factors
        while k <= m / 2 {
            let mut k_combinations = LexicographicCombinationsWithRemovals::new(n, k);
            if 2 * k == m {
                // When m is even and k = m/2 we only need to iterate over half the subsets of size k
                k_combinations.exclude(n - 1); // Since modular factor n-1 is checked last in the lexographic ordering, this trick works
            }
            for i in &excluded_modular_factors {
                // Exclude any previously found modular factors from the search
                k_combinations.exclude(*i);
            }
            loop {
                match k_combinations.next() {
                    Some(subset) => {
                        if dminusone_test.test(&subset) {
                            continue;
                        }

                        let g = Polynomial::mul(
                            &Polynomial::constant(self.leading_coeff.clone()),
                            modular_factor_product_memory_stack.get_product(&subset),
                        )
                        .apply_map(|c| {
                            let c = c.rem(&self.modulus);
                            if c > Integer::quo(&self.modulus, &Integer::TWO).unwrap() {
                                c - &self.modulus
                            } else {
                                c.clone()
                            }
                        })
                        .primitive_part() //factoring f(x) = 49x^2-10000 had possible_factor = 49x-700, which is only a factor over the rationals and not over the integers unless we take the primitive part which is 7x-100, soo this seems to make sense though I cant properly justify it right now.
                        .unwrap();
                        debug_assert_ne!(g.degree().unwrap(), 0);

                        match Polynomial::div(&f, &g) {
                            Ok(h) => {
                                //f = gh
                                // Update variables for next loop:
                                // Divide f by the found factor g
                                f = h;
                                // Add g to this list of found factors
                                factored.mul_mut(Factored::factored_irreducible_unchecked(
                                    Polynomial::<Integer>::structure(),
                                    g,
                                ));
                                // Remove the modular factors for g from future consideration
                                m -= k;
                                for i in subset {
                                    excluded_modular_factors.push(i);
                                    k_combinations.exclude(i);
                                }
                            }
                            Err(RingDivisionError::NotDivisible) => {}
                            Err(RingDivisionError::DivideByZero) => unreachable!(),
                        }
                    }
                    None => {
                        break;
                    }
                }
            }
            k += 1;
        }

        // The remaining modular factors must give the last irreducible factor in the factorization
        if m > 0 {
            factored.mul_mut(Factored::factored_irreducible_unchecked(
                Polynomial::<Integer>::structure(),
                f,
            ));
        }

        factored
    }
}

/// Factor an integer polynomial using a naive implementation of the Berlekamp-Zassenhaus algorithm.
/// No optimizations are used when searching for combinations of modular factors yielding true factors.
pub fn factorize_by_berlekamp_zassenhaus_algorithm(
    poly: Polynomial<Integer>,
) -> Option<Factored<PolynomialStructure<CannonicalStructure<Integer>>>> {
    if poly.is_zero() {
        None
    } else {
        Some(
            Polynomial::<Integer>::structure()
                .factorize_using_primitive_sqfree_factorize_by_yuns_algorithm(poly, &|f| {
                    if f.degree().unwrap() == 0 {
                        Factored::factored_unit_unchecked(Polynomial::<Integer>::structure(), f)
                    } else {
                        let state = BerlekampAassenhausAlgorithmState::new(f).next_prime();
                        state.factor_by_try_all_subsets()
                    }
                }),
        )
    }
}

/// Find a factor of a primitive squarefree integer polynomial by a naive implementation of the Berlekamp-Zassenhaus algorithm.
fn find_factor_primitive_sqfree_by_berlekamp_zassenhaus_algorithm_naive(
    f: Polynomial<Integer>,
) -> FindFactorResult<PolynomialStructure<CannonicalStructure<Integer>>> {
    let f_deg = f.degree().unwrap();
    debug_assert_ne!(f_deg, 0);
    let factor_coeff_bound = f.mignotte_factor_coefficient_bound().unwrap();
    let minimum_modolus = Natural::TWO * factor_coeff_bound;

    if f_deg == 1 {
        FindFactorResult::Irreducible
    } else {
        let prime_gen = PrimeGenerator::new();
        for p in prime_gen {
            let mod_p = QuotientStructure::new_field(Integer::structure(), Integer::from(&p));
            let poly_mod_p = PolynomialStructure::new(mod_p.into());
            if poly_mod_p.degree(&f).unwrap() == f_deg {
                let facotred_f_mod_p = poly_mod_p.factor(&f).unwrap();
                match facotred_f_mod_p.into_hensel_factorization(f.clone()) {
                    Some(hensel_factorization_f_over_p) => {
                        let mut hensel_factorization_f_over_p =
                            hensel_factorization_f_over_p.dont_lift_bezout_coeffs();
                        while hensel_factorization_f_over_p.modolus() < minimum_modolus {
                            hensel_factorization_f_over_p.linear_lift();
                        }

                        let modulus = hensel_factorization_f_over_p.modolus();

                        let modular_factors = hensel_factorization_f_over_p.factors();

                        for subset in (0..modular_factors.len())
                            .map(|_i| vec![false, true])
                            .multi_cartesian_product()
                        {
                            let possible_factor = Polynomial::mul(
                                &Polynomial::constant(f.leading_coeff().unwrap()),
                                &Polynomial::product(
                                    (0..modular_factors.len())
                                        .filter(|i| subset[*i])
                                        .map(|i| modular_factors[i])
                                        .collect(),
                                ),
                            )
                            .apply_map(|c| {
                                let c = c.rem(&modulus);
                                if c > Integer::quo(&modulus, &Integer::TWO).unwrap() {
                                    c - &modulus
                                } else {
                                    c.clone()
                                }
                            })
                            .primitive_part() //factoring f(x) = 49x^2-10000 had possible_factor = 49x-700, which is only a factor over the rationals and not over the integers unless we take the primitive part which is 7x-100, soo this seems to make sense though I cant properly justify it right now.
                            .unwrap();

                            if possible_factor.degree().unwrap() != 0
                                && possible_factor.degree().unwrap() != f_deg
                            {
                                match Polynomial::div(&f, &possible_factor) {
                                    Ok(other_factor) => {
                                        return FindFactorResult::Composite(
                                            possible_factor,
                                            other_factor,
                                        )
                                    }
                                    Err(RingDivisionError::NotDivisible) => {}
                                    Err(RingDivisionError::DivideByZero) => unreachable!(),
                                }
                            }
                        }
                        return FindFactorResult::Irreducible;
                    }
                    None => {}
                }
            }
        }
        unreachable!("Because there are infinitely many primes")
    }
}

/// Factor an integer polynomial using a naive implementation of the Berlekamp-Zassenhaus algorithm to find one factor at a time.
/// No optimizations are used when searching for combinations of modular factors yielding true factors.
pub fn factorize_by_berlekamp_zassenhaus_algorithm_naive(
    f: Polynomial<Integer>,
) -> Option<Factored<PolynomialStructure<CannonicalStructure<Integer>>>> {
    if f.is_zero() {
        None
    } else {
        Some(
            Polynomial::<Integer>::structure()
                .factorize_using_primitive_sqfree_factorize_by_yuns_algorithm(f, &|f| {
                    factorize_by_find_factor(&Polynomial::<Integer>::structure(), f, &|f| {
                        find_factor_primitive_sqfree_by_berlekamp_zassenhaus_algorithm_naive(f)
                    })
                }),
        )
    }
}
