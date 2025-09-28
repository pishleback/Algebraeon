use algebraeon_nzq::{Natural, Rational};
use algebraeon_rings::natural::factorization::factor;
use algebraeon_rings::parsing::parse_rational_polynomial;
use algebraeon_rings::polynomial::Polynomial;
use gungraun::{library_benchmark, library_benchmark_group, main};
use std::hint::black_box;
use std::str::FromStr;

#[library_benchmark]
#[benches::small(iter = (0u8..=20).map(Natural::from))]
#[bench::large1(Natural::from_str("706000565581575429997696139445280900").unwrap())]
fn bench_factor_natural(value: Natural) {
    black_box(factor(value));
}

library_benchmark_group!(
    name = bench_factor_natural_group;
    benchmarks = bench_factor_natural
);

#[library_benchmark]
#[bench::poly1(parse_rational_polynomial("x^5 - x + 1", "x").unwrap())]
fn bench_count_real_polynomial_roots(polynomial: Polynomial<Rational>) {
    black_box(polynomial.count_real_roots());
}

#[library_benchmark]
#[bench::poly1(parse_rational_polynomial("x^5 - x + 1", "x").unwrap())]
fn bench_isolate_real_polynomial_roots(polynomial: Polynomial<Rational>) {
    black_box(polynomial.all_real_roots());
}

#[library_benchmark]
#[bench::poly1(parse_rational_polynomial("x^5 - x + 1", "x").unwrap())]
fn bench_isolate_complex_polynomial_roots(polynomial: Polynomial<Rational>) {
    black_box(polynomial.all_complex_roots());
}

library_benchmark_group!(
    name = bench_count_polynomial_roots;
    benchmarks =
        bench_count_real_polynomial_roots,
        bench_isolate_real_polynomial_roots,
        bench_isolate_complex_polynomial_roots
);

main!(
    library_benchmark_groups = bench_factor_natural_group,
    bench_count_polynomial_roots
);
