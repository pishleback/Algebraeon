use algebraeon_nzq::{Integer, Natural, Rational};
use algebraeon_rings::natural::NaturalFns;
use algebraeon_rings::parsing::{parse_integer_polynomial, parse_rational_polynomial};
use algebraeon_rings::polynomial::Polynomial;
use algebraeon_rings::structure::MetaFactorableSignature;
use gungraun::{library_benchmark, library_benchmark_group, main};
use std::hint::black_box;
use std::str::FromStr;

#[library_benchmark]
#[benches::small(iter = (100u16..=120).map(Natural::from))]
#[bench::large1(Natural::from_str("706000565581575429997696139445280900").unwrap())]
fn bench_factor_natural(value: Natural) {
    black_box(value.factor());
}

library_benchmark_group!(
    name = bench_factor_natural_group;
    benchmarks = bench_factor_natural
);

#[library_benchmark]
#[bench::quadratic_1(parse_integer_polynomial("(x - 2) * (x + 3)", "x").unwrap())]
#[bench::many_factors_1(parse_integer_polynomial("(x - 1)^5 * (x^2 + 1)", "x").unwrap())]
#[bench::many_factors_2(parse_integer_polynomial("x^12 - 1", "x").unwrap())]
#[bench::many_factors_3(parse_integer_polynomial("x^36 - 1", "x").unwrap())]
#[bench::many_factors_4(parse_integer_polynomial("(x-2)*(x-3)*(x-5)*(x-7)*(x-11)*(x-13)*(x-17)", "x").unwrap())]
#[bench::irreducible_1(parse_integer_polynomial("x^6 + x^5 + x^4 + x^3 + x^2 + x + 1", "x").unwrap())]
#[bench::two_irreducibles_1(parse_integer_polynomial("(x^6 + 2*x^5 - x^4 + x^2 - x + 3) * (x^5 - x^4 + 3*x^3 - 2*x + 1)", "x").unwrap())]
fn bench_factor_integer_polynomial(polynomial: Polynomial<Integer>) {
    black_box(polynomial.factor());
}

library_benchmark_group!(
    name = bench_factor_integer_polynomial_group;
    benchmarks = bench_factor_integer_polynomial
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
    bench_factor_integer_polynomial_group,
    bench_count_polynomial_roots
);
