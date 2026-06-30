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
*/

use crate::num_theory::modulo::montgomery::MontgomeryModuloOddPrimeStructure;
use crate::num_theory::natural_factorization::primes::is_prime_nat;
use crate::polynomial::*;
use crate::polynomial::{
    hensel_lifting_btree::HenselFactorization as ProductTreeHenselFactorization,
    hensel_lifting_linalg::HenselFactorization as LinearAlgebraHenselFactorization,
};
use crate::structure::*;
use crate::{
    matrix::{Matrix, StandardInnerProduct},
    structure::EuclideanRemainderQuotientStructure,
};
use algebraeon_sets::combinatorics::LexicographicSubsetsWithRemovals;
use algebraeon_structures::*;
use itertools::Itertools;
use std::collections::BTreeSet;
use std::ops::Rem;

fn compute_polynomial_factor_bound(poly: &Polynomial<Integer>) -> Natural {
    poly.mignotte_factor_coefficient_bound().unwrap()
}

/// Which Hensel-lifting implementation backs a factorization: the linear-algebra
/// lifter for low-degree inputs, or the product-tree lifter for high-degree inputs
/// with many modular factors (where the product tree gives subquadratic lifting).
// Only ever one backend is held at a time (one per factorization), so the size of
// the larger variant is not worth an extra heap indirection on every access.
#[allow(clippy::large_enum_variant)]
enum HenselBackend {
    LinearAlgebra(
        LinearAlgebraHenselFactorization<
            IntegerCanonicalStructure,
            EuclideanRemainderQuotientStructure<
                IntegerCanonicalStructure,
                IntegerCanonicalStructure,
                false,
            >,
        >,
    ),
    ProductTree(ProductTreeHenselFactorization<true, IntegerCanonicalStructure>),
}

/// A Hensel-lifted modular factorization of the squarefree primitive input `f` at
/// a good prime `p`, together with everything needed to recombine the modular
/// factors into the integer factors of `f`.
struct StateAtGoodPrime<'a> {
    p: usize,
    sqfree_prim_poly: &'a Polynomial<Integer>,
    hensel_factorization: HenselBackend,
    num_modular_factors: usize,
    can_quadratic_lift: bool,
}

/// The factorization of the squarefree primitive input `f` modulo a good prime
/// `p`, before any Hensel lifting. A prime is *good* when it does not divide the
/// leading coefficient (so the degree is preserved mod `p`) and `f mod p` is still
/// squarefree; the modular factors are then the distinct monic irreducibles of
/// `f mod p`.
struct FactorizationAtGoodPrime {
    p: usize,
    mod_p: MontgomeryModuloOddPrimeStructure,
    factors: NonZeroFactored<Polynomial<u64>, Natural>,
}

impl FactorizationAtGoodPrime {
    /// Factor `sqfree_prim_poly` modulo the prime `p`.
    ///
    /// Returns `None` when `p` is a bad prime: `p = 2` (unsupported by the
    /// Montgomery backend), `p` divides the leading coefficient (the degree drops
    /// mod `p`), or the reduction `f mod p` is not squarefree. Otherwise returns
    /// the distinct monic irreducible factors of `f mod p`.
    fn try_new(p: usize, sqfree_prim_poly: &Polynomial<Integer>) -> Option<Self> {
        if p == 2 {
            // Montgomery form can only handle odd primes
            // mod p=2 is efficient to implement so shouldn't just skip it here though... future optimization
            return None;
        }
        let p_nat = Natural::from(p);
        debug_assert!(is_prime_nat(&p_nat));
        let mut mod_p = MontgomeryModuloOddPrimeStructure::new_unchecked(p as u64);
        mod_p.populate_inv_cache();
        let poly_mod_p = mod_p.polynomials();
        let sqfree_prim_poly_mod_p = sqfree_prim_poly.apply_map(|c| mod_p.project_ref(c));
        if poly_mod_p.degree(&sqfree_prim_poly_mod_p) == sqfree_prim_poly.degree() {
            let fs = poly_mod_p
                .factorize_monic(&sqfree_prim_poly_mod_p)
                .unwrap()
                .factorize_squarefree();
            if !fs.is_squarefree() {
                return None;
            }

            // not sure at what p it makes most sense to switch algorithm?
            let fs = if p < 50 {
                fs.factorize_berlekamps()
            } else {
                fs.factorize_distinct_degree().factorize_cantor_zassenhaus()
            }
            .unwrap_nonzero();

            Some(Self {
                p,
                mod_p,
                factors: fs,
            })
        } else {
            None
        }
    }

    /// Degrees of the irreducible factors of `f mod p`. The degree multiset of any
    /// true integer factor of `f` is a sub-sum of these, which the caller uses to
    /// rule out impossible factor degrees.
    fn factor_degrees(&self) -> Vec<usize> {
        let polynomial_ring = self.mod_p.polynomials();
        self.factors
            .powers()
            .iter()
            .map(|(factor, _)| polynomial_ring.degree(factor).unwrap())
            .collect()
    }

    /// Pick a Hensel-lift backend and build the lifting state from this modular
    /// factorization. Degree `>= 600` uses the product-tree backend (subquadratic
    /// lifting when there are many factors); smaller degrees use the linear-algebra
    /// backend. The modular factors become the base (mod `p`) of the lift.
    fn into_hensel_state<'a>(
        self,
        sqfree_prim_poly: &'a Polynomial<Integer>,
    ) -> StateAtGoodPrime<'a> {
        let FactorizationAtGoodPrime { p, mod_p, factors } = self;
        let num_modular_factors = factors.powers().len();
        let hensel_factorization = if sqfree_prim_poly.degree().unwrap() >= 600 {
            let integer_factors = factors
                .into_powers()
                .into_iter()
                .map(|(factor, exponent)| {
                    debug_assert_eq!(exponent, Natural::ONE);
                    factor.apply_map(|c| mod_p.unproject_ref(c))
                })
                .collect::<Vec<_>>();
            HenselBackend::ProductTree(ProductTreeHenselFactorization::new_from_mod_field_factors(
                Integer::structure(),
                mod_p,
                sqfree_prim_poly.clone(),
                integer_factors,
            ))
        } else {
            let polynomial_ring_mod_p = mod_p.polynomials();
            HenselBackend::LinearAlgebra(
                LinearAlgebraHenselFactorization::from_mod_field_factorization(
                    |q| {
                        Integer::structure()
                            .into_euclidean_quotient_ring(q.clone())
                            .unwrap()
                    },
                    &polynomial_ring_mod_p.factorizations(),
                    factors,
                    sqfree_prim_poly.clone(),
                )
                .unwrap(),
            )
        };
        StateAtGoodPrime {
            p,
            sqfree_prim_poly,
            num_modular_factors,
            hensel_factorization,
            can_quadratic_lift: true,
        }
    }
}

impl<'a> StateAtGoodPrime<'a> {
    /// Hensel-lift the modular factorization until the modulus reaches
    /// `target_modulus`, i.e. until `f ≡ lc(f) · ∏ factors (mod p^k)` for some
    /// `p^k ≥ target_modulus`. Quadratic lifting doubles the exponent each step
    /// (preferred); the linear-algebra backend falls back to linear lifting
    /// (exponent `+1`) once its Bézout cofactors can no longer be lifted, after
    /// which only linear steps remain available.
    fn lift_to_modulus(&mut self, target_modulus: &Natural) {
        match &mut self.hensel_factorization {
            HenselBackend::LinearAlgebra(factorization) => {
                while self.can_quadratic_lift
                    && factorization.modulus() < target_modulus
                    && factorization.modulus() * factorization.modulus() < target_modulus
                {
                    factorization.quadratic_lift();
                }
                while factorization.modulus() < target_modulus {
                    factorization.linear_lift();
                    self.can_quadratic_lift = false;
                }
            }
            HenselBackend::ProductTree(factorization) => {
                while factorization.modulus() < target_modulus {
                    factorization.quadratic_lift();
                }
            }
        }
    }

    /// The current lifting modulus `p^k`.
    fn modulus(&self) -> Integer {
        match &self.hensel_factorization {
            HenselBackend::LinearAlgebra(factorization) => factorization.modulus().clone(),
            HenselBackend::ProductTree(factorization) => factorization.modulus(),
        }
    }

    /// The current lifted factors as integer polynomials with coefficients reduced
    /// modulo `p^k`.
    fn modular_factors(&self) -> Vec<Polynomial<Integer>> {
        match &self.hensel_factorization {
            HenselBackend::LinearAlgebra(factorization) => factorization.factors().clone(),
            HenselBackend::ProductTree(factorization) => {
                factorization.factors().into_iter().cloned().collect()
            }
        }
    }
}

/// Balanced (symmetric) residue of `x` modulo `modulus`: the representative lying
/// in `(-modulus/2, modulus/2]`. `Rem` returns the residue in `[0, modulus)`,
/// which this re-centres by subtracting `modulus` from the upper half.
fn symmetric_remainder(x: &Integer, modulus: &Integer) -> Integer {
    let r = Rem::rem(x, modulus);
    if r > Integer::quo(modulus, &Integer::TWO).unwrap() {
        r - modulus
    } else {
        r
    }
}

/// `p^exponent` as a `Natural`.
fn prime_power(p: usize, exponent: usize) -> Natural {
    Natural::from(p).pow(&Natural::from(exponent))
}

/// Power sums (Newton traces) of the roots of a monic factor, scaled by powers of
/// the input's leading coefficient and reduced modulo `modulus`.
///
/// For a monic factor `g` of degree `d` with roots `α_1, …, α_d`, this returns
/// `traces[i] = leading_coeff^i · (Σ_j α_j^i) (mod modulus)` for `1 ≤ i ≤
/// max_trace` (`traces[0]` is unused and left `0`). The power sums `Σ_j α_j^i` are
/// obtained from the coefficients of `g` by Newton's identities, so no roots are
/// ever computed. Scaling by `lc(f)^i` turns the power sums into p-adic integers
/// congruent to the corresponding traces of a true integer factor of `f`; these
/// scaled traces are the data van Hoeij's lattice reduction operates on.
fn scaled_newton_traces_mod(
    factor: &Polynomial<Integer>,
    max_trace: usize,
    leading_coeff: &Integer,
    modulus: &Integer,
) -> Vec<Integer> {
    let degree = factor.degree().unwrap();
    let mut raw_traces = vec![Integer::ZERO; max_trace + 1];
    let mut scaled_traces = vec![Integer::ZERO; max_trace + 1];
    let mut leading_coeff_power = Integer::ONE;

    for i in 1..=max_trace {
        let mut sum = Integer::ZERO;
        if i <= degree {
            sum += Integer::from(i) * factor.coeff(degree - i).as_ref();
        }
        for k in 1..=std::cmp::min(i - 1, degree) {
            sum += factor.coeff(degree - k).as_ref() * &raw_traces[i - k];
        }
        raw_traces[i] = Rem::rem(&(-sum), modulus);
        leading_coeff_power = Rem::rem(&(leading_coeff_power * leading_coeff), modulus);
        scaled_traces[i] = Rem::rem(&(&raw_traces[i] * &leading_coeff_power), modulus);
    }

    scaled_traces
}

/// Upper bound on the magnitude of the `trace_index`-th scaled trace of any factor
/// of `poly`.
///
/// If `g` divides `poly` and has roots `α_j`, then `|lc(poly)^i · Σ_j α_j^i| ≤
/// deg(poly) · R^i`, where `R ≥ |lc(poly) · α|` bounds every scaled root. `R` is
/// Fujiwara's root bound scaled by the leading coefficient and rounded up to a
/// power of two, so the bound is computed from bit lengths alone and never takes a
/// large integer root.
fn scaled_trace_bound(poly: &Polynomial<Integer>, trace_index: usize) -> Natural {
    let degree = poly.degree().unwrap();
    let leading_coeff_abs = Abs::abs(poly.leading_coeff().unwrap());
    // Fujiwara's root bound, scaled by the leading coefficient:
    // |lc(f) * alpha| <= 2 max_i (|a_i| |lc(f)|^(n-i-1))^(1/(n-i)).
    // A power-of-two overestimate avoids computing hundreds of large integer roots.
    let leading_coeff_bits = leading_coeff_abs.bits().count();
    let scaled_root_bound_exponent = (0..degree)
        .map(|i| {
            let codimension = degree - i;
            let coeff_bits = Abs::abs(poly.coeff(i).as_ref()).bits().count();
            let radicand_bits = coeff_bits + (codimension - 1) * leading_coeff_bits;
            radicand_bits.div_ceil(codimension)
        })
        .max()
        .unwrap_or(0)
        + 1;
    let scaled_root_bound = Natural::TWO.pow(&Natural::from(scaled_root_bound_exponent));
    Natural::from(degree) * scaled_root_bound.pow(&Natural::from(trace_index))
}

/// Smallest exponent `k` with `p^k > 2 · scaled_trace_bound(poly, trace_index)`.
///
/// A scaled trace lies in `(-bound, bound]`, so knowing it modulo such a `p^k`
/// determines its exact integer value. This is therefore the p-adic precision the
/// Hensel lift must reach before trace number `trace_index` can be recovered.
fn trace_bound_exponent(poly: &Polynomial<Integer>, p: usize, trace_index: usize) -> usize {
    let twice_bound = Natural::TWO * scaled_trace_bound(poly, trace_index);
    let mut exponent = 0;
    let mut power = Natural::ONE;
    let p = Natural::from(p);
    while power <= twice_bound {
        power *= &p;
        exponent += 1;
    }
    exponent
}

/// Extract the high-order p-adic digits of a scaled trace known modulo
/// `p^accuracy_exponent`.
///
/// Writes `trace_mod_pk = high · p^bound_exponent + low` with `low` the balanced
/// residue modulo `p^bound_exponent`, and returns `high` reduced to the balanced
/// residue modulo `p^(accuracy_exponent − bound_exponent)`. The `low` part is the
/// recoverable true value of the trace (it fits below the bound); the returned
/// `high` part is the error that must vanish for a 0/1 combination of modular
/// factors to be a genuine integer factor. These high parts are exactly the
/// entries van Hoeij's knapsack lattice drives to zero.
fn trace_cut(
    trace_mod_pk: &Integer,
    p: usize,
    bound_exponent: usize,
    accuracy_exponent: usize,
) -> Integer {
    debug_assert!(accuracy_exponent > bound_exponent);
    let low_modulus = Integer::from(prime_power(p, bound_exponent));
    let high_modulus = Integer::from(prime_power(p, accuracy_exponent - bound_exponent));
    let low = symmetric_remainder(trace_mod_pk, &low_modulus);
    let high_digits = (trace_mod_pk - low) / low_modulus;
    symmetric_remainder(&high_digits, &high_modulus)
}

/// Reduced row echelon form of `matrix` computed over the rationals, returned as
/// the list of its non-zero rows.
fn rational_rref(matrix: &Matrix<Integer>) -> Vec<Vec<Rational>> {
    let mut rows = (0..matrix.rows())
        .map(|r| {
            (0..matrix.cols())
                .map(|c| Rational::from(matrix.at(r, c).unwrap()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let mut pivot_row = 0;
    for pivot_col in 0..matrix.cols() {
        let Some(nonzero_row) =
            (pivot_row..rows.len()).find(|r| rows[*r][pivot_col] != Rational::ZERO)
        else {
            continue;
        };
        rows.swap(pivot_row, nonzero_row);

        let pivot = rows[pivot_row][pivot_col].clone();
        for entry in &mut rows[pivot_row] {
            *entry = &*entry / &pivot;
        }

        let pivot_row_data = rows[pivot_row].clone();
        for (r, row) in rows.iter_mut().enumerate() {
            if r == pivot_row || row[pivot_col] == Rational::ZERO {
                continue;
            }
            let multiplier = row[pivot_col].clone();
            for (entry, pivot_entry) in row.iter_mut().zip(pivot_row_data.iter()).skip(pivot_col) {
                *entry = &*entry - &multiplier * pivot_entry;
            }
        }

        pivot_row += 1;
        if pivot_row == rows.len() {
            break;
        }
    }

    rows.truncate(pivot_row);
    rows
}

/// Interpret a lattice basis whose columns are indexed by the modular factors as a
/// partition of those factors.
///
/// Returns `Some` only when the reduced row echelon form is a 0/1 block matrix —
/// every column has exactly one `1` and is otherwise `0` — in which case row `r`
/// collects the indices of the modular factors making up block `r`. Returns `None`
/// when the basis does not describe such a clean partition, i.e. the recombination
/// has not yet converged.
fn partition_from_lattice_basis(lattice_basis: &Matrix<Integer>) -> Option<Vec<Vec<usize>>> {
    let rref = rational_rref(lattice_basis);
    let mut partition = vec![vec![]; rref.len()];

    for col in 0..lattice_basis.cols() {
        let mut one_row = None;
        for (row_index, row) in rref.iter().enumerate() {
            if row[col] == Rational::ONE {
                if one_row.is_some() {
                    return None;
                }
                one_row = Some(row_index);
            } else if row[col] != Rational::ZERO {
                return None;
            }
        }
        partition[one_row?].push(col);
    }

    if partition.iter().any(Vec::is_empty) {
        None
    } else {
        Some(partition)
    }
}

/// One van Hoeij knapsack / LLL refinement round.
///
/// The current `lattice_basis` has one column per modular factor; its rows are the
/// candidate 0/1 combinations still believed to assemble into true factors. This
/// augments the basis with extra columns holding the high-order trace digits
/// (`trace_cuts`, the output of [`trace_cut`]) and extra rows for the `trace_moduli`,
/// LLL-reduces the result, keeps the vectors short enough to be genuine factor
/// combinations (those with `4 · ‖·‖² ≤ 4M²`), and projects them back to a basis
/// over the modular factors. Every round shrinks the lattice toward the true
/// factor partition. The modular-factor block is scaled by `M = max(1, #factors/2)`
/// so that a short vector keeps small 0/1 entries while still pricing the trace
/// coordinates against the moduli.
fn refine_factor_lattice(
    lattice_basis: Matrix<Integer>,
    trace_cuts: &[Vec<Integer>],
    trace_moduli: &[Integer],
) -> Matrix<Integer> {
    let modular_factor_count = lattice_basis.cols();
    let trace_count = trace_moduli.len();
    debug_assert_eq!(trace_cuts.len(), modular_factor_count);
    debug_assert!(trace_cuts.iter().all(|cuts| cuts.len() == trace_count));

    let scale = Integer::from(std::cmp::max(1, modular_factor_count / 2));
    let knapsack_basis = Matrix::construct(
        lattice_basis.rows() + trace_count,
        modular_factor_count + trace_count,
        |r, c| {
            if r < lattice_basis.rows() {
                if c < modular_factor_count {
                    &scale * lattice_basis.at(r, c).unwrap()
                } else {
                    let trace_index = c - modular_factor_count;
                    let mut value = Integer::ZERO;
                    for (factor_index, factor_cuts) in trace_cuts.iter().enumerate() {
                        value +=
                            lattice_basis.at(r, factor_index).unwrap() * &factor_cuts[trace_index];
                    }
                    value
                }
            } else if c == modular_factor_count + r - lattice_basis.rows() {
                trace_moduli[r - lattice_basis.rows()].clone()
            } else {
                Integer::ZERO
            }
        },
    );

    let (_, reduced_basis) = knapsack_basis.lll_integral_row_reduction_algorithm(
        &StandardInnerProduct::new(Integer::structure()),
        &Rational::from_integers(3, 4),
    );

    let gram_schmidt = reduced_basis
        .clone()
        .apply_map(|x| Rational::from(x))
        .gram_schmidt_row_orthogonalization(&StandardInnerProduct::new(Rational::structure()));
    let four_m_squared = Integer::from(4) * &scale * &scale * Integer::from(modular_factor_count)
        + Integer::from(trace_count * modular_factor_count * modular_factor_count);
    let four_m_squared = Rational::from(four_m_squared);

    let keep_count = (0..gram_schmidt.rows())
        .filter(|r| {
            let norm_squared = gram_schmidt
                .get_row(*r)
                .into_iter()
                .map(|x| &x * &x)
                .sum::<Rational>();
            Rational::from(4) * norm_squared <= four_m_squared
        })
        .max()
        .map(|r| r + 1)
        .expect("the full polynomial always gives a short knapsack vector");

    let projected_rows = (0..keep_count)
        .filter_map(|r| {
            let row = (0..modular_factor_count)
                .map(|c| {
                    let x = reduced_basis.at(r, c).unwrap();
                    debug_assert!(x.divisible(&scale));
                    x / &scale
                })
                .collect::<Vec<_>>();
            if row.iter().all(|x| x == &Integer::ZERO) {
                None
            } else {
                Some(row)
            }
        })
        .collect::<Vec<_>>();

    Matrix::from_rows(projected_rows)
        .row_span()
        .into_row_basis_matrix()
}

impl<'a> StateAtGoodPrime<'a> {
    /// Try to turn a (hopefully converged) lattice basis into integer factors of
    /// `f`.
    ///
    /// Reads the basis as a partition of the modular factors (`None` if it is not a
    /// clean partition). For each block it forms `lc(f) · ∏ modular factors`,
    /// reduces the coefficients to balanced residues modulo `p^k`, takes the
    /// primitive part as a candidate factor, and trial-divides it into the
    /// remaining part of `f`. It succeeds only if every block divides out exactly
    /// and the final quotient is a unit — a check over `Z` that certifies the
    /// factors regardless of the modulus, so it is always sound to attempt at the
    /// current precision.
    fn try_recombine_lattice_basis(
        &self,
        lattice_basis: &Matrix<Integer>,
    ) -> Option<Vec<Polynomial<Integer>>> {
        let partition = partition_from_lattice_basis(lattice_basis)?;
        let modulus = self.modulus();
        let modular_factors = self.modular_factors();
        let leading_coeff = self.sqfree_prim_poly.leading_coeff().unwrap();
        let mut remaining = self.sqfree_prim_poly.clone();
        let mut factors = vec![];

        for subset in partition {
            let modular_product = Polynomial::product(
                &subset
                    .iter()
                    .map(|i| &modular_factors[*i])
                    .collect::<Vec<_>>(),
            );
            let candidate = Polynomial::mul(
                &Polynomial::constant(leading_coeff.clone()),
                &modular_product,
            )
            .apply_map(|c| symmetric_remainder(c, &modulus))
            .primitive_part()
            .unwrap();
            let quotient = Polynomial::try_divide(&remaining, &candidate)?;
            factors.push(candidate);
            remaining = quotient;
        }

        if remaining.is_unit() {
            factors[0].mul_mut(&remaining);
            Some(factors)
        } else {
            None
        }
    }

    /// Factor `f` by van Hoeij's knapsack recombination of its modular factors.
    ///
    /// Each round lifts a batch of scaled traces to enough p-adic precision
    /// ([`trace_bound_exponent`] plus a safety margin), feeds their high-order
    /// digits ([`trace_cut`]) into the LLL lattice ([`refine_factor_lattice`]), and
    /// tries to read off and verify integer factors. It returns as soon as the
    /// lattice collapses to a single block — van Hoeij's certificate that `f` is
    /// irreducible, returned without further lifting — or a trial division
    /// certifies a non-trivial factorization. The exhaustive lift-to-`2·Mignotte`
    /// fallback is only taken for the small linear-algebra backend, where it is
    /// cheap; for the high-degree product-tree backend the bound is astronomically
    /// large, so convergence relies on the precision accumulated across rounds.
    fn factor_by_van_hoeij_knapsack(
        &mut self,
        minimum_modulus: &Natural,
    ) -> Vec<Polynomial<Integer>> {
        let modular_factor_count = self.num_modular_factors;
        let max_trace = self.sqfree_prim_poly.degree().unwrap() / 2;
        let leading_coeff = self.sqfree_prim_poly.leading_coeff().unwrap().clone();
        let mut lattice_basis = Matrix::<Integer>::ident(modular_factor_count);
        let mut trace_index = 1;
        let diagnostics = std::env::var_os("ALGEBRAEON_FACTOR_TRACE").is_some();
        // Two traces per round matches the precision schedule that reliably collapses
        // the recombination lattice (e.g. for the degree-384 case p7). The lift is now
        // fast enough that the extra p-adic precision this requires is affordable, even
        // for the high-degree product-tree backend.
        let trace_batch_size = 2;

        while trace_index <= max_trace {
            let trace_indices = (trace_index
                ..=std::cmp::min(trace_index + trace_batch_size - 1, max_trace))
                .collect::<Vec<_>>();
            let extra_accuracy = std::cmp::max(8, lattice_basis.rows() / trace_batch_size);
            let bound_exponents = trace_indices
                .iter()
                .map(|i| trace_bound_exponent(self.sqfree_prim_poly, self.p, *i))
                .collect::<Vec<_>>();
            let accuracy_exponents = bound_exponents
                .iter()
                .map(|b| b + extra_accuracy + 1)
                .collect::<Vec<_>>();
            let target_exponent = *accuracy_exponents.iter().max().unwrap();
            if diagnostics {
                eprintln!(
                    "van Hoeij: lifting traces {:?} to p-adic exponent {}",
                    trace_indices, target_exponent
                );
            }
            self.lift_to_modulus(&prime_power(self.p, target_exponent));
            if diagnostics {
                eprintln!("van Hoeij: lift complete");
            }

            let modulus = self.modulus();
            let traces = self
                .modular_factors()
                .iter()
                .map(|factor| {
                    scaled_newton_traces_mod(
                        factor,
                        *trace_indices.iter().max().unwrap(),
                        &leading_coeff,
                        &modulus,
                    )
                })
                .collect::<Vec<_>>();
            let trace_cuts = traces
                .iter()
                .map(|factor_traces| {
                    trace_indices
                        .iter()
                        .enumerate()
                        .map(|(j, trace_index)| {
                            trace_cut(
                                &factor_traces[*trace_index],
                                self.p,
                                bound_exponents[j],
                                accuracy_exponents[j],
                            )
                        })
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let trace_moduli = bound_exponents
                .iter()
                .zip(&accuracy_exponents)
                .map(|(b, a)| Integer::from(prime_power(self.p, a - b)))
                .collect::<Vec<_>>();

            let old_rank = lattice_basis.rows();
            if diagnostics {
                eprintln!("van Hoeij: starting LLL at rank {}", old_rank);
            }
            lattice_basis = refine_factor_lattice(lattice_basis, &trace_cuts, &trace_moduli);
            if diagnostics {
                eprintln!(
                    "van Hoeij: traces {:?}, extra accuracy {}, rank {} -> {}",
                    trace_indices,
                    extra_accuracy,
                    old_rank,
                    lattice_basis.rows()
                );
            }

            if let Some(partition) = partition_from_lattice_basis(&lattice_basis) {
                if partition.len() == 1 {
                    // The lattice collapsed to a single block: every modular factor
                    // belongs to the same true factor, so f is irreducible. This is
                    // van Hoeij's irreducibility certificate, and crucially it needs
                    // no further lifting - we already have f, so there is no need to
                    // reconstruct it from the modular factors (which would require
                    // lifting to ~2 * Mignotte, astronomically large for high-degree
                    // inputs with large coefficients such as the degree-972 case).
                    return vec![self.sqfree_prim_poly.clone()];
                }
                // A recombined candidate is accepted only when it divides f exactly
                // over the integers, so attempting recombination at the precision
                // already reached is always sound: a successful trial division is a
                // genuine factor regardless of the modulus. Try this first.
                if let Some(factors) = self.try_recombine_lattice_basis(&lattice_basis) {
                    return factors;
                }
                // Fallback: lift all the way to 2 * Mignotte bound and retry. This is
                // the textbook guarantee but the bound is ~p^(deg / log2 p), which is
                // astronomically large for the high-degree inputs handled by the
                // product-tree backend (lifting there is infeasible). Restrict the
                // exhaustive lift to the linear-algebra backend, where the degree is
                // small enough for it to be cheap; the product-tree backend instead
                // relies on the precision accumulated from successive trace lifts.
                if matches!(&self.hensel_factorization, HenselBackend::LinearAlgebra(_)) {
                    self.lift_to_modulus(minimum_modulus);
                    if let Some(factors) = self.try_recombine_lattice_basis(&lattice_basis) {
                        return factors;
                    }
                }
            }
            trace_index += trace_indices.len();
        }

        panic!("van Hoeij knapsack recombination did not converge")
    }
}

struct MemoryStack<SG: AssociativeCompositionSignature> {
    semigroup: SG,
    modular_factor_values: Vec<SG::Elem>,
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
    prev_calc: Vec<(usize, SG::Elem)>,
}

impl<SG: AssociativeCompositionSignature> MemoryStack<SG> {
    fn new(semigroup: SG, modular_factor_values: Vec<SG::Elem>) -> Self {
        Self {
            semigroup,
            modular_factor_values,
            prev_calc: vec![],
        }
    }

    fn get_val(&self, i: usize) -> &<SG as SetSignature>::Elem {
        &self.modular_factor_values[i]
    }

    fn get_product(&mut self, subset: &[usize]) -> &<SG as SetSignature>::Elem {
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
        type Elem = DMinusOneTestSemigroupElem;

        fn validate_element(&self, _x: &Self::Elem) -> Result<(), String> {
            Ok(())
        }
    }
    impl CompositionSignature for DMinusOneTestSemigroup {
        fn compose(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
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
                Rational::from(algebraeon_structures::Abs::abs(f.leading_coeff().unwrap()))
                    * f.cauchys_root_bound().unwrap();

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
    fn compose(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.mul(a, b)
    }
}
impl AssociativeCompositionSignature for ModularFactorMultSemigrp {}

// Polynomial degree
#[derive(Debug, Clone, PartialEq, Eq)]
struct ModularFactorDegreeSumSemigrp {}
impl Signature for ModularFactorDegreeSumSemigrp {}
impl SetSignature for ModularFactorDegreeSumSemigrp {
    type Elem = usize;
    fn validate_element(&self, _x: &Self::Elem) -> Result<(), String> {
        Ok(())
    }
}
impl CompositionSignature for ModularFactorDegreeSumSemigrp {
    fn compose(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        a + b
    }
}
impl AssociativeCompositionSignature for ModularFactorDegreeSumSemigrp {}

impl<'a> StateAtGoodPrime<'a> {
    // try to find some factors of the polynomial by trying combinations of modular subsets
    fn partial_factor_by_test_modular_subsets(
        &self,
        max_subset_size: usize,
        possible_proper_factor_degrees: &BTreeSet<usize>,
    ) -> Vec<Polynomial<Integer>> {
        let modulus = self.modulus();
        let modular_factors = self.modular_factors();
        let leading_coeff = self.sqfree_prim_poly.leading_coeff().unwrap();

        let n = modular_factors.len();

        let mut dminusone_test =
            dminusone_test::DMinusOneTest::new(&modulus, self.sqfree_prim_poly, &modular_factors);
        let mut modular_factor_product_memory_stack = MemoryStack::new(
            Integer::structure()
                .into_euclidean_quotient_ring(modulus.clone())
                .unwrap()
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
    f: &Polynomial<Integer>,
) -> Factored<Polynomial<Integer>, Natural> {
    debug_assert_ne!(f.degree().unwrap(), 0);
    debug_assert!(Integer::structure().polynomials().is_primitive(f.clone()));

    // We factor f into its modular factors at some primes p
    // There exists some partition of the modular factors yielding the true irreducible factors of f
    // For each prime we can take the set of possible sums of degrees of modular factors mod p - those are the only possible degrees of irreducible factors of f

    let factor_coeff_bound = compute_polynomial_factor_bound(f);
    let minimum_modulus = Natural::TWO * &factor_coeff_bound;

    let mut best_prime_factorization = None;
    let mut good_primes_checked = 0;
    let max_good_primes = 5;
    let diagnostics = std::env::var_os("ALGEBRAEON_FACTOR_TRACE").is_some();
    let mut possible_proper_factor_degrees = (1..f.degree().unwrap()).collect::<BTreeSet<_>>();
    for p in primes() {
        if let Some(prime_factorization) = FactorizationAtGoodPrime::try_new(p, f) {
            good_primes_checked += 1;
            let factor_degrees = prime_factorization.factor_degrees();
            let n = factor_degrees.len();
            if diagnostics {
                eprintln!("prime {p}: {n} modular factors");
            }
            let mut possible_proper_factor_degrees_at_p = BTreeSet::new();
            for d in factor_degrees {
                for e in possible_proper_factor_degrees_at_p.clone() {
                    possible_proper_factor_degrees_at_p.insert(d + e);
                }
                possible_proper_factor_degrees_at_p.insert(d);
            }
            possible_proper_factor_degrees = possible_proper_factor_degrees
                .intersection(&possible_proper_factor_degrees_at_p)
                .copied()
                .collect();

            if possible_proper_factor_degrees.is_empty() {
                return Polynomial::<Integer>::structure()
                    .factorizations()
                    .new_irreducible_unchecked(f.clone());
            }

            if best_prime_factorization
                .as_ref()
                .is_none_or(|best: &FactorizationAtGoodPrime| n < best.factors.powers().len())
            {
                best_prime_factorization = Some(prime_factorization);
            }

            if n <= 18 || good_primes_checked >= max_good_primes {
                break;
            }
        }
    }
    let mut best_p_state = best_prime_factorization.unwrap().into_hensel_state(f);
    let n = best_p_state.num_modular_factors;
    if diagnostics {
        eprintln!("selected prime {} with {n} modular factors", best_p_state.p);
    }
    if matches!(
        &best_p_state.hensel_factorization,
        HenselBackend::LinearAlgebra(_)
    ) {
        best_p_state.lift_to_modulus(&Natural::from(u64::MAX));
    }
    let factors = if best_p_state.num_modular_factors <= 18 {
        let partially_factored = best_p_state.partial_factor_by_test_modular_subsets(
            {
                if n < 10 {
                    4
                } else if n < 15 {
                    3
                } else {
                    2
                }
            },
            &possible_proper_factor_degrees,
        );
        debug_assert!(!partially_factored.is_empty());
        if partially_factored.len() >= 2 {
            let mut factored = Integer::structure().polynomials().factorizations().one();
            for f in partially_factored {
                Integer::structure().polynomials().factorizations().mul_mut(
                    &mut factored,
                    &factorize_primitive_squarefree_by_berlekamp_zassenhaus_algorithm(&f),
                );
            }
            return factored;
        }

        best_p_state.lift_to_modulus(&minimum_modulus);
        best_p_state
            .partial_factor_by_test_modular_subsets(usize::MAX, &possible_proper_factor_degrees)
    } else {
        best_p_state.factor_by_van_hoeij_knapsack(&minimum_modulus)
    };
    Integer::structure().polynomials().factorizations().product(
        &factors
            .into_iter()
            .map(|g| {
                Integer::structure()
                    .polynomials()
                    .factorizations()
                    .new_irreducible_unchecked(g)
            })
            .collect::<Vec<_>>(),
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
                            &poly,
                        );
                    }
                }

                // it may already be square-free but we couldn't quickly prove it
                // use yuns algorithm to reduce to the case of square-free polynomials
                Polynomial::<Integer>::structure()
                    .factorize_using_primitive_sqfree_factorize_by_yuns_algorithm(poly, &|poly| {
                        factorize_primitive_squarefree_by_berlekamp_zassenhaus_algorithm(&poly)
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
    let minimum_modulus = Natural::TWO * factor_coeff_bound;

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
                    while hensel_factorization_f_over_p.modulus() < minimum_modulus {
                        hensel_factorization_f_over_p.linear_lift();
                    }
                    let modulus = hensel_factorization_f_over_p.modulus();
                    let modular_factors = hensel_factorization_f_over_p.factors();
                    for subset in (0..modular_factors.len())
                        .map(|_i| vec![false, true])
                        .multi_cartesian_product()
                    {
                        let possible_factor = Polynomial::mul(
                            &Polynomial::constant(f.leading_coeff().unwrap().clone()),
                            &Polynomial::product(
                                &(0..modular_factors.len())
                                    .filter(|i| subset[*i])
                                    .map(|i| modular_factors[i])
                                    .collect::<Vec<_>>(),
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

#[cfg(test)]
mod van_hoeij_tests {
    use super::*;
    use crate::num_theory::zimmermann_polys::{p1, p2, p3, p4, p5, p6, p7, p8};

    fn poly(coeffs: &[i32]) -> Polynomial<Integer> {
        Polynomial::from_coeffs(coeffs.iter().map(|c| Integer::from(*c)).collect())
    }

    fn int_matrix(rows: &[&[i32]]) -> Matrix<Integer> {
        Matrix::from_rows(
            rows.iter()
                .map(|row| row.iter().map(|c| Integer::from(*c)).collect())
                .collect(),
        )
    }

    #[test]
    fn symmetric_remainder_recentres_into_balanced_range() {
        // [0, m) representatives in the upper half are shifted down by m.
        assert_eq!(
            symmetric_remainder(&Integer::from(3), &Integer::from(10)),
            Integer::from(3)
        );
        assert_eq!(
            symmetric_remainder(&Integer::from(7), &Integer::from(10)),
            Integer::from(-3)
        );
        // m/2 itself is the largest value kept positive.
        assert_eq!(
            symmetric_remainder(&Integer::from(5), &Integer::from(10)),
            Integer::from(5)
        );
        assert_eq!(
            symmetric_remainder(&Integer::from(10), &Integer::from(10)),
            Integer::from(0)
        );
        // Negative inputs reduce into [0, m) first, then re-centre.
        assert_eq!(
            symmetric_remainder(&Integer::from(-7), &Integer::from(10)),
            Integer::from(3)
        );
        assert_eq!(
            symmetric_remainder(&Integer::from(13), &Integer::from(10)),
            Integer::from(3)
        );
        // Odd modulus: balanced range is [-(m-1)/2, (m-1)/2].
        assert_eq!(
            symmetric_remainder(&Integer::from(4), &Integer::from(7)),
            Integer::from(-3)
        );
        assert_eq!(
            symmetric_remainder(&Integer::from(-1), &Integer::from(7)),
            Integer::from(-1)
        );
    }

    #[test]
    fn prime_power_is_p_to_the_exponent() {
        assert_eq!(prime_power(5, 0), Natural::ONE);
        assert_eq!(prime_power(5, 3), Natural::from(125u32));
        assert_eq!(prime_power(2, 10), Natural::from(1024u32));
    }

    #[test]
    fn scaled_newton_traces_match_power_sums() {
        // g = (x - 2)(x - 3) = x^2 - 5x + 6, roots 2 and 3.
        // Power sums: p1 = 5, p2 = 13, p3 = 35, p4 = 97.
        let g = poly(&[6, -5, 1]);
        let huge = Integer::from(1_000_000_000);
        let traces = scaled_newton_traces_mod(&g, 4, &Integer::ONE, &huge);
        assert_eq!(
            traces,
            vec![
                Integer::from(0),
                Integer::from(5),
                Integer::from(13),
                Integer::from(35),
                Integer::from(97),
            ]
        );

        // Scaling by lc = 2 multiplies trace i by 2^i.
        let scaled = scaled_newton_traces_mod(&g, 4, &Integer::from(2), &huge);
        assert_eq!(
            scaled,
            vec![
                Integer::from(0),
                Integer::from(10),
                Integer::from(52),
                Integer::from(280),
                Integer::from(1552),
            ]
        );

        // The traces are genuinely reduced modulo the given modulus.
        let reduced = scaled_newton_traces_mod(&g, 4, &Integer::ONE, &Integer::from(10));
        assert_eq!(
            reduced,
            vec![
                Integer::from(0),
                Integer::from(5),
                Integer::from(3),
                Integer::from(5),
                Integer::from(7),
            ]
        );
    }

    #[test]
    fn scaled_trace_bound_is_an_upper_bound() {
        // f = (x - 2)(x - 3); with lc = 1 the i-th scaled trace is 2^i + 3^i.
        let f = poly(&[6, -5, 1]);
        for (i, actual) in [(1usize, 5u32), (2, 13), (3, 35)] {
            assert!(scaled_trace_bound(&f, i) >= Natural::from(actual));
        }
    }

    #[test]
    fn trace_bound_exponent_brackets_twice_the_bound() {
        let f = poly(&[6, -5, 1]);
        for i in 1..=3 {
            let p = 5;
            let twice_bound = Natural::TWO * scaled_trace_bound(&f, i);
            let e = trace_bound_exponent(&f, p, i);
            // p^e is the smallest power of p strictly exceeding 2 * bound.
            assert!(prime_power(p, e) > twice_bound);
            assert!(e == 0 || prime_power(p, e - 1) <= twice_bound);
        }
    }

    #[test]
    fn trace_cut_extracts_high_padic_digits() {
        // p = 2, low part has 3 digits (mod 8), result kept mod 2^4 = 16.
        // 42 = 5 * 8 + 2  ->  low = 2, high = 5.
        assert_eq!(trace_cut(&Integer::from(42), 2, 3, 7), Integer::from(5));
        // 6 = 1 * 8 - 2   ->  low = -2 (balanced), high = 1.
        assert_eq!(trace_cut(&Integer::from(6), 2, 3, 7), Integer::from(1));
        // Values below the low modulus have no high part.
        assert_eq!(trace_cut(&Integer::from(2), 2, 3, 7), Integer::from(0));
        // -10 = (-1) * 8 - 2  ->  low = -2, high = -1.
        assert_eq!(trace_cut(&Integer::from(-10), 2, 3, 7), Integer::from(-1));
    }

    #[test]
    fn rational_rref_reduces_to_echelon_form() {
        // Full rank: reduces to the identity.
        let full = int_matrix(&[&[2, 4], &[1, 3]]);
        assert_eq!(
            rational_rref(&full),
            vec![
                vec![Rational::ONE, Rational::ZERO],
                vec![Rational::ZERO, Rational::ONE],
            ]
        );
        // Rank deficient: the dependent row is dropped.
        let deficient = int_matrix(&[&[1, 2], &[2, 4]]);
        assert_eq!(
            rational_rref(&deficient),
            vec![vec![Rational::ONE, Rational::from(2)]]
        );
    }

    #[test]
    fn partition_from_lattice_basis_reads_off_blocks() {
        // Identity basis: every factor is its own singleton block.
        assert_eq!(
            partition_from_lattice_basis(&Matrix::<Integer>::ident(3)),
            Some(vec![vec![0], vec![1], vec![2]])
        );
        // Two clean blocks {0,1} and {2,3}.
        let blocks = int_matrix(&[&[1, 1, 0, 0], &[0, 0, 1, 1]]);
        assert_eq!(
            partition_from_lattice_basis(&blocks),
            Some(vec![vec![0, 1], vec![2, 3]])
        );
        // Overlapping support is not a partition.
        let overlapping = int_matrix(&[&[1, 1, 0], &[0, 1, 1]]);
        assert_eq!(partition_from_lattice_basis(&overlapping), None);
    }

    #[test]
    fn refine_factor_lattice_merges_factors_with_cancelling_traces() {
        // Two modular factors whose scaled-trace high parts cancel only when taken
        // together (a + (-a) = 0) must be merged into a single block, while neither
        // alone gives an integer trace. One LLL round should collapse the identity
        // lattice to the single vector (1, 1).
        let lattice = Matrix::<Integer>::ident(2);
        let a = Integer::from(1000);
        let trace_cuts = vec![vec![a.clone()], vec![-a]];
        let trace_moduli = vec![Integer::from(1 << 20)];
        let refined = refine_factor_lattice(lattice, &trace_cuts, &trace_moduli);
        assert_eq!(refined.rows(), 1);
        assert_eq!(
            partition_from_lattice_basis(&refined),
            Some(vec![vec![0, 1]])
        );
    }

    #[test]
    fn factorization_at_good_prime_splits_mod_p() {
        // x^4 + 3x^2 + 2 = (x^2 + 1)(x^2 + 2) splits into four linear factors mod 17.
        let f = poly(&[2, 0, 3, 0, 1]);
        let factorization = FactorizationAtGoodPrime::try_new(17, &f).unwrap();
        let mut degrees = factorization.factor_degrees();
        degrees.sort_unstable();
        assert_eq!(degrees, vec![1, 1, 1, 1]);
        // p = 2 is rejected as a bad prime.
        assert!(FactorizationAtGoodPrime::try_new(2, &f).is_none());
    }

    #[test]
    fn lift_to_modulus_reaches_target_and_preserves_product() {
        let f = poly(&[2, 0, 3, 0, 1]);
        let mut state = FactorizationAtGoodPrime::try_new(17, &f)
            .unwrap()
            .into_hensel_state(&f);
        assert_eq!(state.num_modular_factors, 4);

        let target = prime_power(17, 8);
        state.lift_to_modulus(&target);
        assert!(state.modulus() >= target);

        // f is monic, so the lifted modular factors multiply back to f modulo p^k.
        let modulus = state.modulus();
        let modular_factors = state.modular_factors();
        let product = Polynomial::product(&modular_factors.iter().collect::<Vec<_>>())
            .apply_map(|c| symmetric_remainder(c, &modulus));
        assert!(Polynomial::<Integer>::structure().equal(&product, &f));
    }

    #[test]
    fn van_hoeij_recombines_a_reducible_polynomial() {
        // (x^2 + 1)(x^2 + 2) should be recovered as two degree-2 factors.
        let f = poly(&[2, 0, 3, 0, 1]);
        let minimum_modulus = Natural::TWO * compute_polynomial_factor_bound(&f);
        let mut state = FactorizationAtGoodPrime::try_new(17, &f)
            .unwrap()
            .into_hensel_state(&f);
        let factors = state.factor_by_van_hoeij_knapsack(&minimum_modulus);

        let mut degrees = factors
            .iter()
            .map(|g| g.degree().unwrap())
            .collect::<Vec<_>>();
        degrees.sort_unstable();
        assert_eq!(degrees, vec![2, 2]);

        let product = Polynomial::product(&factors.iter().collect::<Vec<_>>());
        assert!(Polynomial::<Integer>::structure().equal(&product, &f));
    }

    #[test]
    fn van_hoeij_certifies_irreducibility() {
        // x^4 + 1 is irreducible over Z but splits into four linear factors mod 17,
        // so the lattice must collapse to a single block and return f unchanged.
        let f = poly(&[1, 0, 0, 0, 1]);
        let minimum_modulus = Natural::TWO * compute_polynomial_factor_bound(&f);
        let mut state = FactorizationAtGoodPrime::try_new(17, &f)
            .unwrap()
            .into_hensel_state(&f);
        let factors = state.factor_by_van_hoeij_knapsack(&minimum_modulus);
        assert_eq!(factors.len(), 1);
        assert!(Polynomial::<Integer>::structure().equal(&factors[0], &f));
    }

    fn assert_factor_degrees(poly: Polynomial<Integer>, expected_degrees: &[usize]) {
        let polynomial_ring = Polynomial::<Integer>::structure();
        let factorization = polynomial_ring.factor(&poly);
        assert!(polynomial_ring.equal(
            &poly,
            &polynomial_ring.factorizations().expand(&factorization)
        ));

        let mut actual_degrees = factorization
            .powers()
            .unwrap()
            .iter()
            .flat_map(|(factor, exponent)| {
                let degree = factor.degree().unwrap();
                let multiplicity: usize = exponent.try_into().unwrap();
                std::iter::repeat_n(degree, multiplicity)
            })
            .collect::<Vec<_>>();
        actual_degrees.sort_unstable();
        assert_eq!(actual_degrees, expected_degrees);
    }

    macro_rules! zimmermann_factor_test {
        ($name:ident, $poly:ident, [$($degree:expr),* $(,)?]) => {
            #[test]
            #[ignore = "long-running integer polynomial factorization benchmark"]
            fn $name() {
                assert_factor_degrees($poly(), &[$($degree),*]);
            }
        };
    }

    zimmermann_factor_test!(
        factor_zimmermann_p1,
        p1,
        [
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8,
            8, 8, 8, 8, 8, 8, 8,
        ]
    );
    zimmermann_factor_test!(
        factor_zimmermann_p2,
        p2,
        [2, 2, 12, 12, 12, 12, 24, 24, 24, 24, 24, 24]
    );
    zimmermann_factor_test!(
        factor_zimmermann_p3,
        p3,
        [
            12, 12, 12, 12, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24
        ]
    );
    zimmermann_factor_test!(factor_zimmermann_p4, p4, [66, 396]);
    zimmermann_factor_test!(factor_zimmermann_p5, p5, [64]);
    zimmermann_factor_test!(factor_zimmermann_p6, p6, [12, 12, 12, 12, 48, 48]);
    zimmermann_factor_test!(factor_zimmermann_p7, p7, [384]);
    zimmermann_factor_test!(factor_zimmermann_p8, p8, [972]);
}
