//! A univariate polynomial over a finite field, stored as a coefficient vector

use crate::Field;

use super::traits::UnivariatePolynomial;

pub struct UnivariatePoly<F: Field>(Vec<F>);

impl<F: Field> UnivariatePolynomial<F> for UnivariatePoly<F> {
    fn new(input: &[F]) -> Self {
        let mut last_nonzero = None;
        for (i, item) in input.iter().enumerate().rev() {
            if *item != F::zero() {
                last_nonzero = Some(i);
                break;
            }
        }
        let vec = match last_nonzero {
            Some(i) => input[0..=i].to_vec(),
            None => Vec::new(), // Slice is all 0, return the empty vector
        };
        Self(vec)
    }

    fn coefficients(&self) -> Vec<F> {
        self.0.clone()
    }
    /// Degree of a polynomial.
    ///
    /// Note: degree of the zero polynomial is None
    fn degree(&self) -> Option<usize> {
        if self.0.is_empty() {
            None
        } else {
            Some(self.0.len() - 1)
        }
    }

    /// Evaluation via Horner's method.
    fn evaluate(&self, point: &F) -> F {
        // Handle empty polynomial case
        if self.0.is_empty() {
            return F::zero();
        }

        // Start with the highest degree coefficient
        let mut result = self.0.last().unwrap().clone();

        // Apply Horner's method, iterating through coefficients in reverse order
        for coeff in self.0.iter().rev().skip(1) {
            result *= point;
            result += coeff;
        }
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // Constants for our test primes
    const SMALL_PRIME: u32 = 17;
    const LARGE_PRIME: u32 = 1_073_741_789; // 2^30 - 35
    use crate::MontgomeryFp;

    fn test_empty_polynomial<const P: u32>() {
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[]);
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(5)),
            MontgomeryFp::<P>::zero()
        );
    }

    fn test_constant_polynomial<const P: u32>() {
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[MontgomeryFp::<P>::new(7)]);
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(3)),
            MontgomeryFp::<P>::new(7)
        );
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(0)),
            MontgomeryFp::<P>::new(7)
        );
    }

    fn test_linear_polynomial<const P: u32>() {
        // Polynomial: 3x + 2
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(2),
            MontgomeryFp::<P>::new(3),
        ]);

        // Test at x = 4: 3*4 + 2 = 14
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(4)),
            MontgomeryFp::<P>::new(14)
        );

        // Test at x = 0: 3*0 + 2 = 2
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(0)),
            MontgomeryFp::<P>::new(2)
        );
    }

    fn test_quadratic_polynomial<const P: u32>() {
        // Polynomial: 2x^2 + 3x + 1
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(1),
            MontgomeryFp::<P>::new(3),
            MontgomeryFp::<P>::new(2),
        ]);

        // Test at x = 5: 2*5^2 + 3*5 + 1 = 2*25 + 15 + 1 = 50 + 16 = 66
        let x = 5;
        let expected = MontgomeryFp::<P>::new(66);
        assert_eq!(poly.evaluate(&MontgomeryFp::<P>::new(x)), expected);
    }

    fn test_higher_degree_polynomial<const P: u32>() {
        // Polynomial: 3x^3 + 2x^2 + 5x + 1
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(1),
            MontgomeryFp::<P>::new(5),
            MontgomeryFp::<P>::new(2),
            MontgomeryFp::<P>::new(3),
        ]);

        // Test at x = 2
        let x = 2;
        // expected is 3*8 + 2*4 + 5*2 + 1 = 24 + 8 + 10 + 1 = 43
        let expected = MontgomeryFp::<P>::new(43);

        assert_eq!(poly.evaluate(&MontgomeryFp::<P>::new(x)), expected);
    }

    fn test_degree<const P: u32>() {
        // Test polynomial degrees
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(1),
            MontgomeryFp::<P>::new(2),
            MontgomeryFp::<P>::new(3),
        ]); // 3x^2 + 2x + 1
        assert_eq!(poly.degree(), Some(2));

        let constant = UnivariatePoly::<MontgomeryFp<P>>::new(&[MontgomeryFp::<P>::new(5)]); // 5
        assert_eq!(constant.degree(), Some(0));

        let empty = UnivariatePoly::<MontgomeryFp<P>>::new(&[]);
        assert_eq!(empty.degree(), None);
    }

    fn test_evaluate_at_large_point<const P: u32>() {
        // Polynomial: 4x^2 + 2x + 9
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(9),
            MontgomeryFp::<P>::new(2),
            MontgomeryFp::<P>::new(4),
        ]);

        // x \equiv -2 \bmod P
        let x = P - 2;
        // Expected result is 4x^2 + 2x +5 = 16 - 4 + 9 = 21
        let expected = MontgomeryFp::<P>::new(21);

        assert_eq!(poly.evaluate(&MontgomeryFp::<P>::new(x)), expected);
    }

    fn test_zero_polynomial<const P: u32>() {
        // Polynomial: 0x^2 + 0x + 0
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(0),
        ]);

        // Should evaluate to 0 at any point
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(3)),
            MontgomeryFp::<P>::new(0)
        );
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(7)),
            MontgomeryFp::<P>::new(0)
        );

        // A polynomial with all zeros should be reduced to empty
        assert_eq!(poly.degree(), None);
    }

    fn test_monic_polynomial<const P: u32>() {
        // Monic polynomial: 1x^3 + 0x^2 + 0x + 0 (just x^3)
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(1),
        ]);

        // Test at x = 2: 2^3 = 8
        let x = 2;
        let expected = MontgomeryFp::<P>::new(8);
        assert_eq!(poly.evaluate(&MontgomeryFp::<P>::new(x)), expected);

        // Test at a larger value
        let x = 123;
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(x)),
            MontgomeryFp::<P>::new(x.pow(3))
        );
    }

    fn test_vanishing_polynomial<const P: u32>() {
        // Create a polynomial that vanishes at a specific point
        // For example, (x - a) vanishes at x = a

        let a = 12345 % P; // Ensure a is within range

        // Polynomial: x - a
        // Coefficients: [-a, 1]
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(P - a),
            MontgomeryFp::<P>::new(1),
        ]);

        // The polynomial should evaluate to zero at x = a
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(a)),
            MontgomeryFp::<P>::zero()
        );

        // And non-zero elsewhere
        assert_ne!(
            poly.evaluate(&MontgomeryFp::<P>::new(a + 1)),
            MontgomeryFp::<P>::zero()
        );
    }

    fn test_leading_zero_coefficients<const P: u32>() {
        // Polynomial with trailing zero coefficients
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(5),
            MontgomeryFp::<P>::new(3),
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(0),
        ]);

        // p(2) should be 3*2 + 5 = 11 in any case
        let x = 2;
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(x)),
            MontgomeryFp::<P>::new(11)
        );

        // Degree check - the constructor should strip trailing zeros
        assert_eq!(poly.degree(), Some(1));

        // Check that the internal representation only has the necessary coefficients
        assert_eq!(poly.coefficients().len(), 2);
    }

    fn test_large_coefficients<const P: u32>() {
        // Polynomial with coefficients close to the field size
        // (P-1)x^2 + (P-2)x + (P-3)
        let poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(P - 3),
            MontgomeryFp::<P>::new(P - 2),
            MontgomeryFp::<P>::new(P - 1),
        ]);

        // Test at x = P - 1 (which is congruent to -1 mod P)
        // This forces the evaluation to use the field arithmetic correctly
        // With x = P-1, we get:
        // (P-1)(P-1)^2 + (P-2)(P-1) + (P-3) mod P
        // \equiv -1 + 2 - 3 \equiv -2
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(P - 1)),
            MontgomeryFp::<P>::new(P - 2)
        );
    }

    fn test_trailing_zeros_constructor<const P: u32>() {
        // Test that constructor correctly strips trailing zeros
        let poly1 = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(1),
            MontgomeryFp::<P>::new(2),
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(0),
        ]);

        let poly2 = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(1),
            MontgomeryFp::<P>::new(2),
        ]);

        // Both should have the same coefficients after construction
        assert_eq!(poly1.coefficients(), poly2.coefficients());
        assert_eq!(poly1.degree(), Some(1));

        // All zeros should create an empty polynomial
        let zero_poly = UnivariatePoly::<MontgomeryFp<P>>::new(&[
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(0),
        ]);

        assert_eq!(zero_poly.coefficients(), &[]);
        assert_eq!(zero_poly.degree(), None);
    }

    // Tests for small prime
    #[test]
    fn small_prime_empty_polynomial() {
        test_empty_polynomial::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_constant_polynomial() {
        test_constant_polynomial::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_linear_polynomial() {
        test_linear_polynomial::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_quadratic_polynomial() {
        test_quadratic_polynomial::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_higher_degree_polynomial() {
        test_higher_degree_polynomial::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_degree() {
        test_degree::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_evaluate_at_large_point() {
        test_evaluate_at_large_point::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_zero_polynomial() {
        test_zero_polynomial::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_monic_polynomial() {
        test_monic_polynomial::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_vanishing_polynomial() {
        test_vanishing_polynomial::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_leading_zero_coefficients() {
        test_leading_zero_coefficients::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_large_coefficients() {
        test_large_coefficients::<SMALL_PRIME>();
    }

    #[test]
    fn small_prime_trailing_zeros_constructor() {
        test_trailing_zeros_constructor::<SMALL_PRIME>();
    }

    // Tests for large prime
    #[test]
    fn large_prime_empty_polynomial() {
        test_empty_polynomial::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_constant_polynomial() {
        test_constant_polynomial::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_linear_polynomial() {
        test_linear_polynomial::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_quadratic_polynomial() {
        test_quadratic_polynomial::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_higher_degree_polynomial() {
        test_higher_degree_polynomial::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_degree() {
        test_degree::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_evaluate_at_large_point() {
        test_evaluate_at_large_point::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_zero_polynomial() {
        test_zero_polynomial::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_monic_polynomial() {
        test_monic_polynomial::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_vanishing_polynomial() {
        test_vanishing_polynomial::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_leading_zero_coefficients() {
        test_leading_zero_coefficients::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_large_coefficients() {
        test_large_coefficients::<LARGE_PRIME>();
    }

    #[test]
    fn large_prime_trailing_zeros_constructor() {
        test_trailing_zeros_constructor::<LARGE_PRIME>();
    }
}
