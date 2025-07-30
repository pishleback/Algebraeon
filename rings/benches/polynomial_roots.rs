use algebraeon_nzq::Natural;
use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use std::hint::black_box;

mod sampling;

fn count_real_roots(c: &mut Criterion) {
    let mut group = c.benchmark_group("count_real_roots");

    for degree in [1, 2, 3, 4, 5, 10, 20] {
        let num = 10;
        group.throughput(Throughput::Elements(num));

        group.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, &degree| {
                b.iter(|| {
                    for polynomial in
                        sampling::sample_integer_polynomials(0, num as usize, degree, 10)
                    {
                        black_box(polynomial.count_real_roots());
                    }
                });
            },
        );
    }

    group.finish();
}

fn isolate_real_roots(c: &mut Criterion) {
    let mut group = c.benchmark_group("isolate_real_roots");

    for degree in [1, 2, 3, 4, 5, 10, 20] {
        let num = 10;
        group.throughput(Throughput::Elements(num));

        group.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, &degree| {
                b.iter(|| {
                    for polynomial in
                        sampling::sample_integer_polynomials(0, num as usize, degree, 10)
                    {
                        black_box(polynomial.all_real_roots());
                    }
                });
            },
        );
    }

    group.finish();
}

fn isolate_complex_roots(c: &mut Criterion) {
    let mut group = c.benchmark_group("isolate_complex_roots");

    for degree in [1, 2, 3, 4, 5, 10] {
        let num = 10;
        group.throughput(Throughput::Elements(num));

        group.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, &degree| {
                b.iter(|| {
                    for polynomial in
                        sampling::sample_integer_polynomials(0, num as usize, degree, 10)
                    {
                        black_box(polynomial.all_complex_roots());
                    }
                });
            },
        );
    }

    group.finish();
}

fn isolate_2adic_roots(c: &mut Criterion) {
    let mut group = c.benchmark_group("isolate_2adic_roots");

    for degree in [1, 2, 3, 4, 5, 10, 20] {
        let num = 10;
        group.throughput(Throughput::Elements(num));

        group.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, &degree| {
                b.iter(|| {
                    for polynomial in
                        sampling::sample_integer_polynomials(0, num as usize, degree, 10)
                    {
                        black_box(polynomial.all_padic_roots(&Natural::from(2u32)));
                    }
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    benches,
    count_real_roots,
    isolate_real_roots,
    isolate_complex_roots,
    isolate_2adic_roots
);
criterion_main!(benches);
