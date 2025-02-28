# Finite Field and Sumcheck Protocol Implementation

A Rust implementation of fundamental components for zero-knowledge proof systems, focusing on the Sumcheck protocol and multilinear extensions.

## Overview

This project implements several fundamental building blocks for modern zero-knowledge proof systems:

1. **Finite Field Arithmetic**: Montgomery arithmetic implementation for prime fields
2. **Polynomial Representations**:
   - Univariate polynomials
   - General multivariate polynomials (sparse representation)
   - Multilinear polynomials (for efficient multilinear extensions)
3. **Sumcheck Protocol**: Complete implementation of prover and verifier components

## Project Structure

```text
src/
├── ff/                  # Finite field implementations
│   ├── mod.rs
│   ├── small_prime_mont.rs  # Montgomery arithmetic for small primes
│   └── traits.rs        # Field trait definitions
├── poly/                # Polynomial implementations
│   ├── general_multivariate.rs
│   ├── mod.rs
│   ├── multilinear.rs   # Multilinear extension implementation
│   ├── traits.rs
│   └── univariate_poly.rs
├── sumcheck/            # Sumcheck protocol implementation
│   ├── mod.rs
│   └── protocol.rs      # Prover and Verifier logic
└── main.rs
```

## Features

### Montgomery Arithmetic

The `MontgomeryFp<P>` type provides an efficient implementation of prime field arithmetic using Montgomery representation. It supports:

- Basic field operations (+, -, *, /)
- Field inversion
- Exponentiation
- Random element generation

Example:
```rust
use fabric_challenge::{Field, MontgomeryFp};

// Define a type for a specific prime field
type F = MontgomeryFp<251>;  // Field with prime modulus 251

// Field operations
let a = F::from(123);
let b = F::from(45);
let c = a * b;  // Multiplication in the field
```

### Polynomial Representations

Multiple polynomial representations are provided for different use cases:

- `UnivariatePolynomial<F>`: Single-variable polynomials
- `GeneralMultivariatePolynomial<F>`: Sparse representation of multivariate polynomials using a HashMap to store non-zero terms
- `MultilinearPolynomial<F>`: Polynomials with degree at most 1 in each variable, using evaluation-based representation

Example:
```rust
use fabric_challenge::{Field, MontgomeryFp, MultilinearPolynomial, MultivariatePolynomial};
use std::collections::HashMap;

type F = MontgomeryFp<251>;

// Create a multilinear polynomial from evaluations on the boolean hypercube
let evaluations = vec![
    F::from(1), F::from(2), F::from(3), F::from(4),
    F::from(5), F::from(6), F::from(7), F::from(8)
];
let poly = MultilinearPolynomial::from_evaluations_on_hypercube(evaluations);

// Evaluate the polynomial at a point
let point = vec![F::from(10), F::from(20), F::from(30)];
let result = poly.evaluate(&point);
```

### Sumcheck Protocol

The `sumcheck` module provides implementations of the Prover and Verifier for the Sumcheck protocol:

- `Prover<F, P>`: Generates proofs that the sum of a polynomial over the boolean hypercube equals a claimed value
- `Verifier<F, P>`: Verifies such proofs with high probability

The implementation supports arbitrary multivariate polynomials that implement the `MultivariatePolynomial` trait.

## Usage

### Sumcheck Protocol Execution

The Sumcheck protocol execution can be observed in the test cases in `src/sumcheck/protocol.rs`. Here's a simplified example:

```rust
use fabric_challenge::{Field, MontgomeryFp, GeneralMultivariatePolynomial, MultivariatePolynomial, Prover, Verifier, VerifierMessage};
use std::collections::HashMap;

// Create a polynomial f(x,y,z) = 5 + 3xy + 2yz^2
let mut coeffs = HashMap::new();
coeffs.insert(vec![0, 0, 0], MontgomeryFp::<251>::from(5)); // Constant term
coeffs.insert(vec![1, 1, 0], MontgomeryFp::<251>::from(3)); // xy term
coeffs.insert(vec![0, 1, 2], MontgomeryFp::<251>::from(2)); // yz^2 term
let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);

// Create prover and verifier
let mut prover = Prover::new(poly.clone());
let mut verifier = Verifier::new(poly, Some(42)); // Use a fixed seed for deterministic testing

// Run the protocol
let mut verifier_msg = VerifierMessage::Initial;
while matches!(verifier_msg, VerifierMessage::Challenge(_) | VerifierMessage::Initial) {
    let prover_msg = prover.message(verifier_msg).unwrap();
    verifier_msg = verifier.message(prover_msg).unwrap();
}

// Check the final verifier message
match verifier_msg {
    VerifierMessage::Accept => println!("Verification succeeded!"),
    VerifierMessage::Reject(reason) => println!("Verification failed: {}", reason),
    _ => unreachable!(),
}
```

## Requirements

- Rust 1.54 or higher
- Standard Rust toolchain (cargo, rustc)

## Building and Testing

```bash
# Build the project
cargo build

# Run tests
cargo test

# Run specific tests
cargo test test_successful_protocol_execution
```

## License

MIT