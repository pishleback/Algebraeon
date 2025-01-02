// https://en.wikipedia.org/wiki/Berlekamp%27s_algorithm
/*!
 * # The Berlekamp Zassenhaus algorithm
 * The Berlekamp Zassenhaus algorithm (BZA) is a method for factoring integer polynomials.
 * BZA works on primitive squarefree polynomials, so other algorithms must be used first to reduce the factorization problem to this case.
 *
 * Let $f(x) \in \mathbb{Z}[x]$ be non-zero primitive squarefree polynomial of degree at least $1$
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

use crate::{number::integer::*, polynomial::polynomial::*, ring_structure::quotient::*};
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
                Some(mut hensel_factorization_f_over_p) => {
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

trait SemiGrpTestInfo {
    fn compose(a: &Self, b: &Self) -> Self;
}

struct BZATestManager<G: SemiGrpTestInfo> {
    modular_factors_semigrp: Vec<G>,
    prev_subset: Vec<usize>,
    prev_semigrp: G,
}

impl<G: SemiGrpTestInfo> BZATestManager<G> {
    fn test(&self, subset: &Vec<usize>) {}
}

impl BerlekampZassenhausAlgorithmStateAtPrime {
    fn factor_by_try_all_subsets<'a>(
        &'a self,
    ) -> Factored<PolynomialStructure<CannonicalStructure<Integer>>> {
        // println!("{:?}", self);

        /*
         * Loop over subsets of modular factors in order of increasing cardinality, then in lexographic order for each cardinality
         * Exclude half the subsets by ignoring complimentary subsets
         */

        // let mod_m = QuotientStructure::new_ring(Integer::structure(), self.modulus.clone());
        // let poly_mod_m = PolynomialStructure::new(mod_m.into());

        let modular_factors: Vec<_> = self.modular_factors.iter().collect();
        let n = modular_factors.len();

        // self.modular_factors
        //     .sort_by_cached_key(|mf| poly_mod_m.degree(mf).unwrap());
        // println!("{:?}", self.modular_factors);

        let subset_to_possible_factor = |subset: &Vec<usize>| {
            Polynomial::mul(
                &Polynomial::constant(self.leading_coeff.clone()),
                &Polynomial::product(subset.iter().map(|i| modular_factors[*i]).collect()),
            )
            .apply_map(|c| {
                let c = Integer::rem(c, &self.modulus);
                if c > Integer::quo(&self.modulus, &Integer::TWO).unwrap() {
                    c - &self.modulus
                } else {
                    c.clone()
                }
            })
            .primitive_part() //factoring f(x) = 49x^2-10000 had possible_factor = 49x-700, which is only a factor over the rationals and not over the integers unless we take the primitive part which is 7x-100, soo this seems to make sense though I cant properly justify it right now.
            .unwrap()
        };

        // println!("n = {:?}", n);
        let mut factored = Factored::factored_one(Polynomial::<Integer>::structure());
        let mut excluded_modular_factors = vec![];
        let mut f = self.poly.clone();
        let mut k = 1; // The cardinality of the subset to search for each loop
        let mut m = n; // Keep track of the number of remaining modular factors i.e. m = n - excluded_modular_factors.len()

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
                        // println!("{:?}", subset);
                        let g = subset_to_possible_factor(&subset);
                        debug_assert_ne!(g.degree().unwrap(), 0);

                        // println!("possible_factor = {}", g);
                        // println!("{:?}", Polynomial::div(&f, &g));

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
                .factorize_by_primitive_sqfree_factorize_by_yuns_algorithm(poly, &|f| {
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
                // println!("f mod {} = {:?}", p, facotred_f_mod_p);
                match facotred_f_mod_p.into_hensel_factorization(f.clone()) {
                    Some(mut hensel_factorization_f_over_p) => {
                        // println!("good prime = {}", p);
                        // println!("{:?}", hensel_factorization_f_over_p.factors());

                        while hensel_factorization_f_over_p.modolus() < minimum_modolus {
                            hensel_factorization_f_over_p.linear_lift();
                        }

                        let modulus = hensel_factorization_f_over_p.modolus();

                        let modular_factors = hensel_factorization_f_over_p.factors();

                        // println!("{:?}", hensel_factorization_f_over_p.modolus());
                        // println!("{:?}", lifted_factors);
                        // println!("lifted_factors.len() = {:?}", lifted_factors.len());

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
                                let c = Integer::rem(c, &modulus);
                                if c > Integer::quo(&modulus, &Integer::TWO).unwrap() {
                                    c - &modulus
                                } else {
                                    c.clone()
                                }
                            })
                            .primitive_part() //factoring f(x) = 49x^2-10000 had possible_factor = 49x-700, which is only a factor over the rationals and not over the integers unless we take the primitive part which is 7x-100, soo this seems to make sense though I cant properly justify it right now.
                            .unwrap();

                            // println!("possible_factor = {:?}", possible_factor);

                            // println!("{:?}", Polynomial::div(&f, &possible_factor));

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
                .factorize_by_primitive_sqfree_factorize_by_yuns_algorithm(f, &|f| {
                    factorize_by_find_factor(&Polynomial::<Integer>::structure(), f, &|f| {
                        find_factor_primitive_sqfree_by_berlekamp_zassenhaus_algorithm_naive(f)
                    })
                }),
        )
    }
}
