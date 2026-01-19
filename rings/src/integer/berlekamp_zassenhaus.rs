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

use crate::polynomial::*;
use crate::structure::*;
use algebraeon_groups::structure::AssociativeCompositionSignature;
use algebraeon_groups::structure::CompositionSignature;
use algebraeon_nzq::primes;
use algebraeon_nzq::*;
use algebraeon_sets::combinatorics::LexicographicSubsetsWithRemovals;
use algebraeon_sets::structure::*;
use itertools::Itertools;
use std::collections::BTreeSet;
use std::ops::Rem;
use std::usize;

fn compute_polynomial_factor_bound(poly: &Polynomial<Integer>) -> Natural {
    poly.mignotte_factor_coefficient_bound().unwrap()
}

struct BerlekampZassenhausPrimeSearch<'a> {
    sqfree_prim_poly: &'a Polynomial<Integer>,
    prime_gen: Box<dyn Iterator<Item = usize>>,
    at_primes: Vec<BerlekampZassenhausAlgorithmStateAtPrime<'a>>,
    possible_proper_factor_degrees: BTreeSet<usize>,
    best_idx: Option<usize>,
    smallest_n: usize, // the smallest number of modular factors, achived by `at_primes[best_idx]`.
}

impl<'a> BerlekampZassenhausPrimeSearch<'a> {
    fn new(sqfree_prim_poly: &'a Polynomial<Integer>) -> Self {
        debug_assert_ne!(sqfree_prim_poly.degree().unwrap(), 0);
        let prime_gen = Box::new(primes());
        Self {
            sqfree_prim_poly,
            prime_gen,
            at_primes: vec![],
            possible_proper_factor_degrees: (1..sqfree_prim_poly.degree().unwrap())
                .collect::<BTreeSet<_>>(),
            best_idx: None,
            smallest_n: usize::MAX,
        }
    }

    // Yield primes `p` for which `self.poly` mod `p` is of full degree and the factorization of `self.poly` mod `p` remains squarefree
    fn next_good_prime(&mut self) {
        loop {
            let p = Natural::from(self.prime_gen.next().unwrap());
            if let Some(state_at_p) = BerlekampZassenhausAlgorithmStateAtPrime::try_new_at_prime(
                self.sqfree_prim_poly,
                p.clone(),
            ) {
                let modular_factors = state_at_p.hensel_factorization.factors();
                let n = modular_factors.len();

                // is `n` the smallest so far?
                if self.best_idx.is_none_or(|_| n < self.smallest_n) {
                    self.smallest_n = n;
                    self.best_idx = Some(self.at_primes.len());
                }

                // update the list of possible degrees of true factors
                let mut possible_proper_factor_degrees_at_p = BTreeSet::new();
                for mf in &modular_factors {
                    let d = mf.degree().unwrap();
                    for e in possible_proper_factor_degrees_at_p.clone() {
                        possible_proper_factor_degrees_at_p.insert(d + e);
                    }
                    possible_proper_factor_degrees_at_p.insert(d);
                }
                self.possible_proper_factor_degrees = self
                    .possible_proper_factor_degrees
                    .intersection(&possible_proper_factor_degrees_at_p)
                    .copied()
                    .collect();

                self.at_primes.push(state_at_p);
                return;
            }
        }
    }

    fn at_most_recent_prime(
        &mut self,
    ) -> Option<&mut BerlekampZassenhausAlgorithmStateAtPrime<'a>> {
        self.at_primes.last_mut()
    }

    fn at_best_prime(&self) -> Option<&BerlekampZassenhausAlgorithmStateAtPrime<'a>> {
        Some(&self.at_primes[self.best_idx?])
    }
}

#[derive(Default, Debug, Clone)]
enum BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization {
    #[default]
    Empty,
    Quadratic(HenselFactorization<true, IntegerCanonicalStructure>),
    Linear(HenselFactorization<false, IntegerCanonicalStructure>),
}

impl BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization {
    fn modulus(&self) -> Integer {
        match self {
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Empty => unreachable!(),
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Quadratic(hf) => {
                hf.modulus()
            }
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Linear(hf) => hf.modulus(),
        }
    }

    fn factors(&self) -> Vec<&Polynomial<Integer>> {
        match self {
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Empty => unreachable!(),
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Quadratic(hf) => {
                hf.factors()
            }
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Linear(hf) => hf.factors(),
        }
    }

    fn lift_to_modulus(&mut self, target_modulus: &Natural) {
        // For hensel lifting it is more efficient (due to how the coefficients grow during Hensel lifting) to begin
        // with quadratic lifts until the modulus is just below the size limit for primitive integers, and then proceed with linear lifts

        // quadratic lifting
        match self {
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Empty => unreachable!(),
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Quadratic(
                hensel_factorization,
            ) => {
                while hensel_factorization
                    .factorization_modulus_base()
                    .nat_pow(&(Natural::TWO * hensel_factorization.factorization_modulus_power()))
                    < Natural::from(u64::MAX)
                {
                    if hensel_factorization.modulus() >= target_modulus {
                        return;
                    }
                    hensel_factorization.quadratic_lift();
                }
            }
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Linear(_) => {}
        }

        *self = match std::mem::take(self) {
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Empty => unreachable!(),
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Quadratic(
                hensel_factorization,
            ) => BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Linear(
                hensel_factorization.dont_lift_bezout_coeffs(),
            ),
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Linear(hf) => {
                BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Linear(hf)
            }
        };

        // linear lifting
        match self {
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Empty => unreachable!(),
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Quadratic(_) => {
                unreachable!()
            }
            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Linear(
                hensel_factorization,
            ) => loop {
                if hensel_factorization.modulus() >= target_modulus {
                    return;
                }
                hensel_factorization.linear_lift();
            },
        }
    }
}

#[derive(Debug, Clone)]
struct BerlekampZassenhausAlgorithmStateAtPrime<'a> {
    p: Natural,
    sqfree_prim_poly: &'a Polynomial<Integer>,
    hensel_factorization: BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization,
}

impl<'a> BerlekampZassenhausAlgorithmStateAtPrime<'a> {
    fn try_new_at_prime(sqfree_prim_poly: &'a Polynomial<Integer>, p: Natural) -> Option<Self> {
        debug_assert!(p.is_irreducible());
        let mod_p = Integer::structure().into_quotient_field_unchecked(Integer::from(&p));
        let poly_mod_p = mod_p.polynomials();
        if poly_mod_p.degree(sqfree_prim_poly) == sqfree_prim_poly.degree() {
            poly_mod_p
                .factorizations()
                .into_hensel_factorization(
                    poly_mod_p.factor(sqfree_prim_poly).unwrap_nonzero(),
                    sqfree_prim_poly.clone(),
                )
                .map(
                    |hensel_factorization| BerlekampZassenhausAlgorithmStateAtPrime {
                        p,
                        sqfree_prim_poly,
                        hensel_factorization:
                            BerlekampZassenhausAlgorithmStateAtPrimeHenselFactorization::Quadratic(
                                hensel_factorization,
                            ),
                    },
                )
        } else {
            None
        }
    }

    fn lift_to_modulus(&mut self, target_modulus: &Natural) {
        self.hensel_factorization.lift_to_modulus(target_modulus);
    }
}

struct MemoryStack<SG: AssociativeCompositionSignature> {
    semigroup: SG,
    modular_factor_values: Vec<SG::Set>,
    // Store the partial products of a previous calculation
    // Since subsets are visited in lexcographic order, if a test is performed frequently, values towards the right will need to be updated
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

impl<SG: AssociativeCompositionSignature> MemoryStack<SG> {
    fn new(semigroup: SG, modular_factor_values: Vec<SG::Set>) -> Self {
        Self {
            semigroup,
            modular_factor_values,
            prev_calc: vec![],
        }
    }

    fn get_val(&self, i: usize) -> &<SG as SetSignature>::Set {
        &self.modular_factor_values[i]
    }

    fn get_product(&mut self, subset: &[usize]) -> &<SG as SetSignature>::Set {
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
        #[allow(clippy::needless_range_loop)]
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
/// Even better, this summing can be translated to machine arithmetic taking advantage of the wrapping behaviour of binary addition. This is at the cost of a slight lossening of the bound on the (d-1)st coefficient.
mod dminusone_test {
    use algebraeon_nzq::traits::Abs;

    use super::*;

    #[derive(Debug, Clone)]
    struct DMinusOneTestSemigroupElem {
        approx_coeff_lower_bound: usize, // An upper bound is given by this +1
        degree: usize,
    }
    #[derive(Debug, Clone, PartialEq, Eq)]
    struct DMinusOneTestSemigroup {}
    impl Signature for DMinusOneTestSemigroup {}
    impl SetSignature for DMinusOneTestSemigroup {
        type Set = DMinusOneTestSemigroupElem;

        fn is_element(&self, _x: &Self::Set) -> Result<(), String> {
            Ok(())
        }
    }
    impl CompositionSignature for DMinusOneTestSemigroup {
        fn compose(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
            DMinusOneTestSemigroupElem {
                approx_coeff_lower_bound: a
                    .approx_coeff_lower_bound
                    .wrapping_add(b.approx_coeff_lower_bound),
                degree: a.degree.wrapping_add(b.degree),
            }
        }
    }
    impl AssociativeCompositionSignature for DMinusOneTestSemigroup {}

    pub struct DMinusOneTest {
        memory_stack: MemoryStack<DMinusOneTestSemigroup>,
        factor_dminusone_coeff_bound_divideby_gdeg: usize,
    }
    impl DMinusOneTest {
        pub fn new(
            modulus: &Integer,
            f: &Polynomial<Integer>,
            modular_factors: &[Polynomial<Integer>],
        ) -> Self {
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
            let factor_dminusone_coeff_bound_divideby_gdeg =
                Rational::from(f.leading_coeff().unwrap().abs()) * f.cauchys_root_bound().unwrap();

            Self {
                factor_dminusone_coeff_bound_divideby_gdeg:
                    (factor_dminusone_coeff_bound_divideby_gdeg * &conversion_mult)
                        .ceil()
                        .abs()
                        .try_into()
                        .unwrap_or(usize::MAX),
                memory_stack: MemoryStack::new(
                    DMinusOneTestSemigroup {},
                    modular_factors
                        .iter()
                        .map(|g| {
                            let d = g.degree().unwrap();
                            let coeff =
                                (f.leading_coeff().unwrap() * g.coeff(d - 1).as_ref()).rem(modulus);
                            DMinusOneTestSemigroupElem {
                                approx_coeff_lower_bound: (Rational::from(coeff)
                                    * &conversion_mult)
                                    .floor()
                                    .abs()
                                    .rem(&machine_range)
                                    .try_into()
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
        pub fn test(&mut self, subset: &[usize]) -> bool {
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
                .factor_dminusone_coeff_bound_divideby_gdeg
                .saturating_mul(d);
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
type ModularFactorMultSemigrp = PolynomialStructure<
    EuclideanRemainderQuotientStructure<
        IntegerCanonicalStructure,
        IntegerCanonicalStructure,
        false,
    >,
    EuclideanRemainderQuotientStructure<
        IntegerCanonicalStructure,
        IntegerCanonicalStructure,
        false,
    >,
>;
impl CompositionSignature for ModularFactorMultSemigrp {
    fn compose(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.mul(a, b)
    }
}
impl AssociativeCompositionSignature for ModularFactorMultSemigrp {}

// Polynomial degree
#[derive(Debug, Clone, PartialEq, Eq)]
struct ModularFactorDegreeSumSemigrp {}
impl Signature for ModularFactorDegreeSumSemigrp {}
impl SetSignature for ModularFactorDegreeSumSemigrp {
    type Set = usize;
    fn is_element(&self, _x: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}
impl CompositionSignature for ModularFactorDegreeSumSemigrp {
    fn compose(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        a + b
    }
}
impl AssociativeCompositionSignature for ModularFactorDegreeSumSemigrp {}

impl<'a> BerlekampZassenhausAlgorithmStateAtPrime<'a> {
    // try to find some factors of the polynomial by trying combinations of modular subsets
    fn partial_factor_by_test_modular_subsets(
        &self,
        max_subset_size: usize,
        possible_proper_factor_degrees: &BTreeSet<usize>,
    ) -> Vec<Polynomial<Integer>> {
        let modulus = self.hensel_factorization.modulus();
        let modular_factors = self
            .hensel_factorization
            .factors()
            .into_iter()
            .cloned()
            .collect::<Vec<_>>();
        let leading_coeff = self.sqfree_prim_poly.leading_coeff().unwrap();

        let n = modular_factors.len();

        let mut dminusone_test =
            dminusone_test::DMinusOneTest::new(&modulus, self.sqfree_prim_poly, &modular_factors);
        let mut modular_factor_product_memory_stack = MemoryStack::new(
            Integer::structure()
                .into_quotient_ring(modulus.clone())
                .into_polynomials(),
            modular_factors.clone(),
        );
        let mut modular_factor_degree_memory_stack = MemoryStack::new(
            ModularFactorDegreeSumSemigrp {},
            modular_factors
                .iter()
                .map(|mf| mf.degree().unwrap())
                .collect(),
        );

        let mut factors = vec![];
        let mut excluded_modular_factors = vec![];
        let mut f = self.sqfree_prim_poly.clone();
        let mut k = 1; // The cardinality of the subset to search for each loop
        let mut m = n; // Keep track of the number of remaining modular factors i.e. m = n - excluded_modular_factors.len()

        // Loop over subsets of modular factors in order of increasing cardinality, then in lexicographic order for each cardinality
        // Only half of the cardinalities need to be checked since complimentary subsets need not be checked
        // Keep searching for cardinalities up to and including half the number of remaining modular factors
        while k <= std::cmp::min(m / 2, max_subset_size) {
            let mut k_combinations = LexicographicSubsetsWithRemovals::new(n, k);
            if 2 * k == m {
                // When m is even and k = m/2 we only need to iterate over half the subsets of size k
                k_combinations.exclude(n - 1); // Since modular factor n-1 is checked last in the lexicographic ordering, this trick works
            }
            for i in &excluded_modular_factors {
                // Exclude any previously found modular factors from the search
                k_combinations.exclude(*i);
            }
            #[allow(clippy::while_let_loop)]
            loop {
                match k_combinations.next() {
                    Some(subset) => {
                        let d = modular_factor_degree_memory_stack.get_product(&subset);
                        debug_assert_eq!(
                            *d,
                            subset
                                .iter()
                                .map(|i| modular_factors[*i].degree().unwrap())
                                .sum::<usize>()
                        );

                        if dminusone_test.test(&subset) {
                            continue;
                        }

                        if !possible_proper_factor_degrees.contains(d) {
                            continue;
                        }

                        let g = Polynomial::mul(
                            &Polynomial::constant(leading_coeff.clone()),
                            modular_factor_product_memory_stack.get_product(&subset),
                        )
                        .apply_map(|c| {
                            let c = Rem::rem(c, &modulus);
                            if c > Integer::quo(&modulus, &Integer::TWO).unwrap() {
                                c - &modulus
                            } else {
                                c.clone()
                            }
                        })
                        .primitive_part() //factoring f(x) = 49x^2-10000 had possible_factor = 49x-700, which is only a factor over the rationals and not over the integers unless we take the primitive part which is 7x-100, soo this seems to make sense though I cant properly justify it right now.
                        .unwrap();
                        debug_assert_ne!(g.degree().unwrap(), 0);
                        if let Some(h) = Polynomial::try_divide(&f, &g) {
                            //f = gh
                            // Update variables for next loop:
                            // Divide f by the found factor g
                            f = h;
                            // Add g to this list of found factors
                            factors.push(g);
                            // Remove the modular factors for g from future consideration
                            m -= k;
                            for i in subset {
                                excluded_modular_factors.push(i);
                                k_combinations.exclude(i);
                            }
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
            factors.push(f);
        } else {
            // Or maybe there is just a unit left over.
            debug_assert!(f.is_unit());
            factors[0].mul_mut(&f);
        }
        factors
    }
}

fn factorize_primitive_squarefree_by_berlekamp_zassenhaus_algorithm(
    f: Polynomial<Integer>,
) -> Factored<Polynomial<Integer>, Natural> {
    debug_assert_ne!(f.degree().unwrap(), 0);
    debug_assert!(Integer::structure().polynomials().is_primitive(f.clone()));

    // We factor f into its modular factors at some primes p
    // There exists some partition of the modular factors yielding the true irreducible factors of f
    // For each prime we can take the set of possible sums of degrees of modular factors mod p - those are the only possible degrees of irreducible factors of f

    let mut prime_search = BerlekampZassenhausPrimeSearch::new(&f);
    for _ in 1..10 {
        prime_search.next_good_prime();
        if prime_search.possible_proper_factor_degrees.is_empty() {
            return Polynomial::<Integer>::structure()
                .factorizations()
                .new_irreducible_unchecked(f.clone());
        }

        if prime_search.smallest_n < 18 {
            break;
        }

        let possible_proper_factor_degrees = prime_search.possible_proper_factor_degrees.clone();
        let at_prime = prime_search.at_most_recent_prime().unwrap();
        println!("p={}", at_prime.p);
        at_prime.lift_to_modulus(&Natural::from(1000000000u64));
        let partially_factored =
            at_prime.partial_factor_by_test_modular_subsets(3, &possible_proper_factor_degrees);
        debug_assert!(!partially_factored.is_empty());
        if partially_factored.len() >= 2 {
            println!("awunga 1");
            let mut factored = Integer::structure().polynomials().factorizations().one();
            for f in partially_factored {
                Integer::structure().polynomials().factorizations().mul_mut(
                    &mut factored,
                    &factorize_primitive_squarefree_by_berlekamp_zassenhaus_algorithm(f),
                );
            }
            return factored;
        }
    }
    let mut at_prime = prime_search.at_best_prime().unwrap().clone();

    println!(
        "p = {}, n = {}, possible_proper_factor_degrees = {:?}",
        at_prime.p,
        at_prime.hensel_factorization.factors().len(),
        prime_search.possible_proper_factor_degrees
    );

    let factor_coeff_bound = compute_polynomial_factor_bound(at_prime.sqfree_prim_poly);
    let minimum_modolus = Natural::TWO * &factor_coeff_bound;
    at_prime.lift_to_modulus(&minimum_modolus);
    Integer::structure().polynomials().factorizations().product(
        at_prime
            .partial_factor_by_test_modular_subsets(
                usize::MAX,
                &prime_search.possible_proper_factor_degrees,
            )
            .into_iter()
            .map(|g| {
                Integer::structure()
                    .polynomials()
                    .factorizations()
                    .new_irreducible_unchecked(g)
            })
            .collect(),
    )
}

/// Factor an integer polynomial using a naive implementation of the Berlekamp-Zassenhaus algorithm.
/// No optimizations are used when searching for combinations of modular factors yielding true factors.
pub fn factorize_by_berlekamp_zassenhaus_algorithm(
    poly: Polynomial<Integer>,
) -> Factored<Polynomial<Integer>, Natural> {
    if poly.is_zero() {
        Factored::Zero
    } else {
        Polynomial::<Integer>::structure().factorize_by_primitive_factorize(
            poly,
            Integer::factor,
            &|poly| {
                // optimization for square-free inputs
                // if we can prove it's already square-free then we can skip yuns algorithm
                // test if it is square-free modulo some big random primes
                // if it is square-free then it's very likely square-free modulo one of the primes, proving the original is square-free too
                // this is a quicker test because computing the discriminant mod p does not suffer from coefficient blow-up
                const SQUARE_FREE_TEST_PRIMES: [u64; 3] = [179424691, 179424697, 179424719];
                #[cfg(debug_assertions)]
                let disc = poly.clone().discriminant().unwrap();
                for p in SQUARE_FREE_TEST_PRIMES {
                    let mod_p =
                        Integer::structure().into_quotient_field_unchecked(Integer::from(p));
                    let poly_mod_p = mod_p.polynomials();

                    if Integer::structure()
                        .polynomials()
                        .leading_coeff(&poly)
                        .unwrap()
                        .divisible(&Integer::from(p))
                    {
                        continue;
                    }

                    let disc_mod_p = poly_mod_p.discriminant(poly.clone()).unwrap();
                    let disc_mod_p_is_zero = mod_p.is_zero(&disc_mod_p);

                    #[cfg(debug_assertions)]
                    assert_eq!(disc.divisible(&Integer::from(p)), disc_mod_p_is_zero);

                    if !disc_mod_p_is_zero {
                        // the poly is square-free
                        return factorize_primitive_squarefree_by_berlekamp_zassenhaus_algorithm(
                            poly,
                        );
                    }
                }

                // it may already be square-free but we couldn't quickly prove it
                // use yuns algorithm to reduce to the case of square-free polynomials
                Polynomial::<Integer>::structure()
                    .factorize_using_primitive_sqfree_factorize_by_yuns_algorithm(poly, &|poly| {
                        factorize_primitive_squarefree_by_berlekamp_zassenhaus_algorithm(poly)
                    })
            },
        )
    }
}

/// Find a factor of a primitive squarefree integer polynomial by a naive implementation of the Berlekamp-Zassenhaus algorithm.
fn find_factor_primitive_sqfree_by_berlekamp_zassenhaus_algorithm_naive(
    f: Polynomial<Integer>,
) -> FindFactorResult<Polynomial<Integer>> {
    let f_deg = f.degree().unwrap();
    debug_assert_ne!(f_deg, 0);
    let factor_coeff_bound = f.mignotte_factor_coefficient_bound().unwrap();
    let minimum_modolus = Natural::TWO * factor_coeff_bound;

    if f_deg == 1 {
        FindFactorResult::Irreducible
    } else {
        let prime_gen = primes();
        for p in prime_gen {
            let mod_p = Integer::structure().into_quotient_field_unchecked(Integer::from(p));
            let poly_mod_p = mod_p.polynomials();
            if poly_mod_p.degree(&f).unwrap() == f_deg {
                let facotred_f_mod_p = poly_mod_p.factor(&f).unwrap_nonzero();
                if let Some(hensel_factorization_f_over_p) = poly_mod_p
                    .factorizations()
                    .into_hensel_factorization(facotred_f_mod_p, f.clone())
                {
                    let mut hensel_factorization_f_over_p =
                        hensel_factorization_f_over_p.dont_lift_bezout_coeffs();
                    while hensel_factorization_f_over_p.modulus() < minimum_modolus {
                        hensel_factorization_f_over_p.linear_lift();
                    }
                    let modulus = hensel_factorization_f_over_p.modulus();
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
                            let c = Rem::rem(c, &modulus);
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
                            debug_assert!(!possible_factor.is_zero());
                            if let Some(other_factor) = Polynomial::try_divide(&f, &possible_factor)
                            {
                                return FindFactorResult::Composite(possible_factor, other_factor);
                            }
                        }
                    }
                    return FindFactorResult::Irreducible;
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
) -> Factored<Polynomial<Integer>, Natural> {
    if f.is_zero() {
        Factored::Zero
    } else {
        Polynomial::<Integer>::structure().factorize_by_primitive_factorize(
            f,
            Integer::factor,
            &|f| {
                Polynomial::<Integer>::structure()
                    .factorize_using_primitive_sqfree_factorize_by_yuns_algorithm(f, &|f| {
                        factorize_by_find_factor(&Polynomial::<Integer>::structure(), f, &|f| {
                            find_factor_primitive_sqfree_by_berlekamp_zassenhaus_algorithm_naive(f)
                        })
                    })
            },
        )
    }
}
