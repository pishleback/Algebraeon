use algebraeon_nzq::Natural;
use algebraeon_rings::natural::factorization::factor;
use gungraun::{library_benchmark, library_benchmark_group, main};
use std::hint::black_box;
use std::str::FromStr;

fn fibonacci(n: u64) -> u64 {
    match n {
        0 => 1,
        1 => 1,
        // 2 => 2,
        // 3 => 3,
        // 4 => 5,
        n => fibonacci(n - 1) + fibonacci(n - 2),
    }
}

#[library_benchmark]
#[bench::short(10)]
#[bench::long(30)]
fn bench_fibonacci(value: u64) -> u64 {
    black_box(fibonacci(value))
}

library_benchmark_group!(
    name = bench_fibonacci_group;
    benchmarks = bench_fibonacci
);


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

main!(
    library_benchmark_groups = bench_fibonacci_group,
    bench_factor_natural_group
);
