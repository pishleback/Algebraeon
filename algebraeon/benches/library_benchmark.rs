use algebraeon::nzq::{Integer, IntegerCanonicalStructure, Natural, Rational};
use algebraeon::rings::matrix::{Matrix, RealInnerProduct, StandardInnerProduct};
use algebraeon::rings::parsing::{parse_integer_polynomial, parse_rational_polynomial};
use algebraeon::rings::polynomial::Polynomial;
use algebraeon::rings::structure::MetaFactoringMonoid;
use algebraeon::sets::structure::MetaType;
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
// #[bench::zimmermann_p1(algebraeon_rings::integer::zimmermann_polys::p1())]
// #[bench::zimmermann_p2(algebraeon_rings::integer::zimmermann_polys::p2())]
// fn bench_factor_integer_polynomial(polynomial: Polynomial<Integer>) {
//     black_box(polynomial.factor());
// }

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

#[library_benchmark]
#[bench::lll(
    Matrix::<Integer>::from_rows(vec![
            vec![1, 0, 0, 0, 0, 0, 10000000000i64],
            vec![0, 1, 0, 0, 0, 0, 11614471390],
            vec![0, 0, 1, 0, 0, 0, 13489594567],
            vec![0, 0, 0, 1, 0, 0, 15667451017],
            vec![0, 0, 0, 0, 1, 0, 18196916159],
            vec![0, 0, 0, 0, 0, 1, 21134756212],
    ]),
    &StandardInnerProduct::new(Integer::structure()),
    Rational::from_str("3/4").unwrap()
)]
fn bench_lll_integral_dim6(
    mat: Matrix<Integer>,
    inner_product: &impl RealInnerProduct<IntegerCanonicalStructure>,
    delta: Rational,
) {
    black_box(mat.lll_integral_row_reduction_algorithm(inner_product, &delta));
}

library_benchmark_group!(
    name = lll;
    benchmarks =
        bench_lll_integral_dim6,
);

main!(
    library_benchmark_groups = bench_factor_natural_group,
    bench_factor_integer_polynomial_group,
    bench_count_polynomial_roots,
    lll
);
