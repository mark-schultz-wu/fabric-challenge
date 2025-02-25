//! A univariate polynomial over a finite field, stored as a coefficient vector

use crate::Field;

use super::traits::UnivariatePolynomial;

pub struct UnivariatePoly<F: Field>(Vec<F>);

impl<F: Field> UnivariatePolynomial<F> for UnivariatePoly<F> {
    fn coefficients(&self) -> Vec<F> {
        self.0.clone()
    }
    fn degree(&self) -> isize {
        // Note: degree of the 0 polynomial is -1.
        self.0.len() as isize - 1
    }

    /// Evaluation via Horner's method.
    ///
    fn evaluate(&self, point: &F) -> F {
        // Handle empty polynomial case
        if self.0.is_empty() {
            return F::zero();
        }

        // Start with the highest degree coefficient
        let mut result = self.0.last().unwrap().clone();

        // Apply Horner's method, iterating through coefficients in reverse
        // (from second-highest degree to constant term)
        for coeff in self.0.iter().rev().skip(1) {
            // result = result * point + coefficient
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
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![]);
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(5)),
            MontgomeryFp::<P>::zero()
        );
    }

    fn test_constant_polynomial<const P: u32>() {
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![MontgomeryFp::<P>::new(7)]);
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
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![
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
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![
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
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![
            MontgomeryFp::<P>::new(1),
            MontgomeryFp::<P>::new(5),
            MontgomeryFp::<P>::new(2),
            MontgomeryFp::<P>::new(3),
        ]);

        // Test using Horner's method for consistency with the implementation
        let x = 2;
        // expected is 3*8 + 2*4 + 5*2 + 1 = 24 + 8 + 10 + 1 = 43
        let expected = MontgomeryFp::<P>::new(43);

        assert_eq!(poly.evaluate(&MontgomeryFp::<P>::new(x)), expected);
    }

    fn test_degree<const P: u32>() {
        // Test polynomial degrees
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![
            MontgomeryFp::<P>::new(1),
            MontgomeryFp::<P>::new(2),
            MontgomeryFp::<P>::new(3),
        ]); // 3x^2 + 2x + 1
        assert_eq!(poly.degree(), 2);

        let constant = UnivariatePoly::<MontgomeryFp<P>>(vec![MontgomeryFp::<P>::new(5)]); // 5
        assert_eq!(constant.degree(), 0);

        let empty = UnivariatePoly::<MontgomeryFp<P>>(vec![]);
        assert_eq!(empty.degree(), -1);
    }

    fn test_evaluate_at_large_point<const P: u32>() {
        // Polynomial: 4x^2 + 2x + 9
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![
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
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![
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
    }

    fn test_monic_polynomial<const P: u32>() {
        // Monic polynomial: 1x^3 + 0x^2 + 0x + 0 (just x^3)
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(1),
        ]);

        // Test at x = 2: 2^3 = 8
        let x = 2;
        let expected = (x * x * x) % P;
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(x)),
            MontgomeryFp::<P>::new(expected)
        );

        // Test at a larger value
        let x = 123;
        let expected = ((x as u64 * x as u64) % P as u64 * x as u64) % P as u64;
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(x)),
            MontgomeryFp::<P>::new(expected as u32)
        );
    }

    fn test_vanishing_polynomial<const P: u32>() {
        // Create a polynomial that vanishes at a specific point
        // For example, (x - a) vanishes at x = a

        let a = 12345 % P; // Ensure a is within range

        // Polynomial: x - a
        // Coefficients: [-a, 1]
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![
            MontgomeryFp::<P>::new(P - a), // This is -a in the field
            MontgomeryFp::<P>::new(1),
        ]);

        // The polynomial should evaluate to zero at x = a
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(a)),
            MontgomeryFp::<P>::zero()
        );

        // And non-zero elsewhere
        assert_ne!(
            poly.evaluate(&MontgomeryFp::<P>::new((a + 1) % P)),
            MontgomeryFp::<P>::zero()
        );
    }

    fn test_leading_zero_coefficients<const P: u32>() {
        // Polynomial with leading zero coefficients
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![
            MontgomeryFp::<P>::new(5),
            MontgomeryFp::<P>::new(3),
            MontgomeryFp::<P>::new(0),
            MontgomeryFp::<P>::new(0),
        ]);

        // p(2) should be 3*2 + 5 = 11 in any case
        let x = 2;
        let expected = (3 * x + 5) % P;
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(x)),
            MontgomeryFp::<P>::new(expected)
        );

        // Degree check - this depends on the implementation's handling of leading zeros
        assert_eq!(poly.degree(), 3); // Assuming no stripping of leading zeros
    }

    fn test_large_coefficients<const P: u32>() {
        // Polynomial with coefficients close to the field size
        // (P-1)x^2 + (P-2)x + (P-3)
        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![
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

    fn test_interpolation<const P: u32>() {
        // Test polynomial that passes through specific points
        // For 3 points, we need a degree 2 polynomial

        // Let's create a polynomial that passes through:
        // (1, 5), (2, 12), (3, 23)
        // This is the quadratic polynomial: 2x^2 + x + 2

        let poly = UnivariatePoly::<MontgomeryFp<P>>(vec![
            MontgomeryFp::<P>::new(2),
            MontgomeryFp::<P>::new(1),
            MontgomeryFp::<P>::new(2),
        ]);

        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(1)),
            MontgomeryFp::<P>::new(5)
        ); // 2*1^2 + 1*1 + 2 = 5
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(2)),
            MontgomeryFp::<P>::new(12)
        ); // 2*2^2 + 1*2 + 2 = 12
        assert_eq!(
            poly.evaluate(&MontgomeryFp::<P>::new(3)),
            MontgomeryFp::<P>::new(23)
        ); // 2*3^2 + 1*3 + 2 = 23
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
    fn small_prime_interpolation() {
        test_interpolation::<SMALL_PRIME>();
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
    fn large_prime_interpolation() {
        test_interpolation::<LARGE_PRIME>();
    }
}
