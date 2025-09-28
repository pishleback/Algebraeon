use iai_callgrind::{main, library_benchmark, library_benchmark_group};
use std::hint::black_box;

// Your function to benchmark
fn fibonacci(n: u64) -> u64 {
    match n {
        0 => 0,
        1 => 1,
        n => fibonacci(n - 1) + fibonacci(n - 2),
    }
}

// Annotate with library_benchmark
#[library_benchmark]
fn bench_fibonacci() -> u64 {
    // use black_box on inputs/outputs to prevent compiler optimizations
    black_box(fibonacci(black_box(20)))
}

// You can add more with #[library_benchmark] ...

// Group them
library_benchmark_group!(
    name = my_group;
    benchmarks = bench_fibonacci
);

// Use main! with the group
main!(library_benchmark_groups = my_group);