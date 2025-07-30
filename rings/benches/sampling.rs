use algebraeon_nzq::{Integer, Rng};
use algebraeon_rings::polynomial::Polynomial;

pub fn sample_integer_polynomials(
    seed: u128,
    num: usize,
    degree: usize,
    max_coeff: usize,
) -> impl Iterator<Item = Polynomial<Integer>> {
    let mut rng = Rng::new(seed);
    (0..num).map(move |_| {
        let mut coeffs = (0..(degree + 1))
            .map(|_| {
                rng.uniform_random_integer_from_inclusive_range(
                    -Integer::from(max_coeff),
                    Integer::from(max_coeff),
                )
            })
            .collect::<Vec<_>>();
        while coeffs.last().unwrap() == &Integer::ZERO {
            *coeffs.last_mut().unwrap() = rng.uniform_random_integer_from_inclusive_range(
                -Integer::from(max_coeff),
                Integer::from(max_coeff),
            );
        }
        Polynomial::from_coeffs(coeffs)
    })
}
