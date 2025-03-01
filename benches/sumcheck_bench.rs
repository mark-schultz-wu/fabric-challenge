use criterion::{black_box, criterion_group, criterion_main, Criterion};
use fabric_challenge::{
    Field, GeneralMultivariatePolynomial, MontgomeryFp, MultilinearPolynomial,
    MultivariatePolynomial, Prover, Verifier, VerifierMessage,
};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::HashMap;

// Define our test field
type F = MontgomeryFp<10007>;

// Generate a polynomial with specified complexity
fn gen_poly(num_vars: usize, num_terms: usize, seed: u64) -> GeneralMultivariatePolynomial<F> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut coeffs = HashMap::new();

    // Add a constant term
    coeffs.insert(vec![0; num_vars], F::from(1));

    // Add random terms
    for _ in 1..num_terms {
        let mut exponents = vec![0; num_vars];
        for exp in exponents.iter_mut() {
            *exp = rng.random_range(0..=2); // Max degree 2
        }

        // Ensure at least one variable has non-zero exponent
        if exponents.iter().all(|&e| e == 0) {
            let idx = rng.random_range(0..num_vars);
            exponents[idx] = 1;
        }

        coeffs.insert(exponents, F::random(&mut rng));
    }

    GeneralMultivariatePolynomial::from_coefficients(coeffs)
}

// Generate a random evaluation point
fn gen_point(num_vars: usize, seed: u64) -> Vec<F> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..num_vars).map(|_| F::random(&mut rng)).collect()
}

// Benchmark core operations
fn bench_operations(c: &mut Criterion) {
    // Define the test variables
    let num_vars = 8;
    let num_terms = 20;

    let general_poly = gen_poly(num_vars, num_terms, 42);
    let mle_poly = MultilinearPolynomial::from_general_polynomial(&general_poly);
    let point = gen_point(num_vars, 43);

    // 1. Evaluate
    let mut group = c.benchmark_group("Evaluation");
    group.bench_function("GeneralPoly", |b| {
        b.iter(|| black_box(&general_poly).evaluate(black_box(&point)))
    });
    group.bench_function("MultilinearPoly", |b| {
        b.iter(|| black_box(&mle_poly).evaluate(black_box(&point)))
    });
    group.finish();

    // 2. Boolean Hypercube Sum
    let mut group = c.benchmark_group("BooleanHypercubeSum");
    group.bench_function("GeneralPoly", |b| {
        b.iter(|| black_box(&general_poly).sum_over_boolean_hypercube())
    });
    group.bench_function("MultilinearPoly", |b| {
        b.iter(|| black_box(&mle_poly).sum_over_boolean_hypercube())
    });
    group.finish();

    // 3. Univariate Slice
    let mut group = c.benchmark_group("UnivariateSlice");
    group.bench_function("GeneralPoly", |b| {
        b.iter(|| black_box(&general_poly).univariate_slice_last())
    });
    group.bench_function("MultilinearPoly", |b| {
        b.iter(|| black_box(&mle_poly).univariate_slice_last())
    });
    group.finish();

    // 4. Shrink Last
    let mut group = c.benchmark_group("ShrinkLast");
    group.bench_function("GeneralPoly", |b| {
        b.iter_batched(
            || (general_poly.clone(), F::from(123)),
            |(mut poly, value)| {
                let _ = poly.shrink_last(&value);
                poly
            },
            criterion::BatchSize::SmallInput,
        )
    });
    group.bench_function("MultilinearPoly", |b| {
        b.iter_batched(
            || (mle_poly.clone(), F::from(123)),
            |(mut poly, value)| {
                let _ = poly.shrink_last(&value);
                poly
            },
            criterion::BatchSize::SmallInput,
        )
    });
    group.finish();
}

// Benchmark full protocol execution
fn bench_protocol(c: &mut Criterion) {
    let mut group = c.benchmark_group("FullProtocol");

    // Test with 5 variables
    let num_vars = 5;
    let num_terms = 20;

    // General polynomial
    group.bench_function("GeneralPoly", |b| {
        b.iter_batched(
            || {
                let poly = gen_poly(num_vars, num_terms, 42);
                let prover = Prover::new(poly.clone());
                let verifier = Verifier::new(poly, Some(43));
                (prover, verifier)
            },
            |(mut prover, mut verifier)| {
                let mut verifier_msg = VerifierMessage::Initial;
                while matches!(
                    verifier_msg,
                    VerifierMessage::Challenge(_) | VerifierMessage::Initial
                ) {
                    let prover_msg = prover.message(verifier_msg).unwrap();
                    verifier_msg = verifier.message(prover_msg).unwrap();
                }
                verifier_msg
            },
            criterion::BatchSize::SmallInput,
        )
    });

    // Multilinear polynomial
    group.bench_function("MultilinearPoly", |b| {
        b.iter_batched(
            || {
                let general = gen_poly(num_vars, num_terms, 42);
                let poly = MultilinearPolynomial::from_general_polynomial(&general);
                let prover = Prover::new(poly.clone());
                let verifier = Verifier::new(poly, Some(43));
                (prover, verifier)
            },
            |(mut prover, mut verifier)| {
                let mut verifier_msg = VerifierMessage::Initial;
                while matches!(
                    verifier_msg,
                    VerifierMessage::Challenge(_) | VerifierMessage::Initial
                ) {
                    let prover_msg = prover.message(verifier_msg).unwrap();
                    verifier_msg = verifier.message(prover_msg).unwrap();
                }
                verifier_msg
            },
            criterion::BatchSize::SmallInput,
        )
    });

    group.finish();
}

// Benchmark MLE conversion
fn bench_mle_conversion(c: &mut Criterion) {
    let mut group = c.benchmark_group("MLE Conversion");

    let num_vars = 8;
    let num_terms = 20;
    let general_poly = gen_poly(num_vars, num_terms, 42);

    group.bench_function("GeneralToMLE", |b| {
        b.iter(|| MultilinearPolynomial::from_general_polynomial(black_box(&general_poly)))
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_operations,
    bench_protocol,
    bench_mle_conversion
);
criterion_main!(benches);
