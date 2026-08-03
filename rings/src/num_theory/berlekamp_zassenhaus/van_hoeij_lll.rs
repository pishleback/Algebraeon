use super::*;

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
    pub fn factor_by_van_hoeij_knapsack(
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

#[cfg(test)]
mod tests {
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
