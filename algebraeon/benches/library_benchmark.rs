use algebraeon_nzq::Natural;
use algebraeon_rings::natural::factorization::factor;
use gungraun::{library_benchmark, library_benchmark_group, main};
use std::hint::black_box;
use std::str::FromStr;

#[library_benchmark]
#[bench::zero(Natural::from(0u8))]
#[bench::one(Natural::from(1u8))]
#[benches::small(iter = (2u8..=20).map(Natural::from))]
#[bench::large1(Natural::from_str("706000565581575429997696139445280900").unwrap())]
fn bench_factor_natural(value: Natural) {
    black_box(factor(value));
}

library_benchmark_group!(
    name = bench_factor_natural_group;
    benchmarks = bench_factor_natural
);

main!(library_benchmark_groups = bench_factor_natural_group);
