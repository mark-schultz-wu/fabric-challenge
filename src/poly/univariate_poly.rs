use crate::Field;

/// Represents a univariate polynomial
#[derive(Debug, Clone)]
pub struct UnivariatePolynomial<F: Field> {
    /// Coefficients in ascending order of degree
    /// e.g., [2, 3, 1] represents 2 + 3x + x²
    pub coefficients: Vec<F>,
}

impl<F: Field> UnivariatePolynomial<F> {
    /// Creates a new univariate polynomial from coefficients
    pub fn new(coefficients: Vec<F>) -> Self {
        let mut coeffs = coefficients;
        // Ensure `coefficients` is non-empty
        if coeffs.is_empty() {
            coeffs.push(F::zero());
        }
        // Trim trailing zeros
        while coeffs.len() > 1 && coeffs.last() == Some(&F::zero()) {
            coeffs.pop();
        }
        Self {
            coefficients: coeffs,
        }
    }

    /// Evaluates the polynomial at a point
    pub fn evaluate(&self, point: F) -> F {
        // Horner's method for efficient evaluation
        let mut result = F::zero();
        for coeff in self.coefficients.iter().rev() {
            result *= &point;
            result += coeff;
        }
        result
    }

    /// Returns the degree of the polynomial
    pub fn degree(&self) -> Option<usize> {
        for (i, v) in self.coefficients.iter().enumerate().rev() {
            if v.has_no_terms() {
                continue;
            } else {
                return Some(i);
            }
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use crate::Field;
    use crate::MontgomeryFp;
    use crate::UnivariatePolynomial;

    // Using a small prime (251) for testing
    type F = MontgomeryFp<251>;

    #[test]
    fn test_univariate_polynomial_creation() {
        // Test empty polynomial creation
        let empty_poly = UnivariatePolynomial::<F>::new(vec![]);
        assert_eq!(empty_poly.coefficients.len(), 1);
        assert_eq!(empty_poly.coefficients[0], F::zero());

        // Test creation with trailing zeros
        let poly_with_zeros =
            UnivariatePolynomial::new(vec![F::from(1), F::from(2), F::from(0), F::from(0)]);
        assert_eq!(poly_with_zeros.coefficients.len(), 2);
        assert_eq!(poly_with_zeros.coefficients[0], F::from(1));
        assert_eq!(poly_with_zeros.coefficients[1], F::from(2));

        // Test normal polynomial creation
        let poly = UnivariatePolynomial::new(vec![F::from(3), F::from(2), F::from(1)]);
        assert_eq!(poly.coefficients.len(), 3);
        assert_eq!(poly.coefficients[0], F::from(3));
        assert_eq!(poly.coefficients[1], F::from(2));
        assert_eq!(poly.coefficients[2], F::from(1));
    }

    #[test]
    fn test_univariate_polynomial_degree() {
        // Test degree of zero polynomial
        let zero_poly = UnivariatePolynomial::<F>::new(vec![F::zero()]);
        assert_eq!(zero_poly.degree(), None);

        // Test degree of constant polynomial
        let constant_poly = UnivariatePolynomial::new(vec![F::from(5)]);
        assert_eq!(constant_poly.degree(), Some(0));

        // Test degree of linear polynomial
        let linear_poly = UnivariatePolynomial::new(vec![F::from(3), F::from(2)]);
        assert_eq!(linear_poly.degree(), Some(1));

        // Test degree of higher-degree polynomial
        let cubic_poly =
            UnivariatePolynomial::new(vec![F::from(1), F::from(2), F::from(3), F::from(4)]);
        assert_eq!(cubic_poly.degree(), Some(3));
    }

    #[test]
    fn test_univariate_polynomial_evaluation() {
        // Test evaluation of constant polynomial
        let constant_poly = UnivariatePolynomial::new(vec![F::from(5)]);
        assert_eq!(constant_poly.evaluate(F::from(10)), F::from(5));

        // Test evaluation of linear polynomial: 3 + 2x
        let linear_poly = UnivariatePolynomial::new(vec![F::from(3), F::from(2)]);
        assert_eq!(linear_poly.evaluate(F::from(0)), F::from(3));
        assert_eq!(linear_poly.evaluate(F::from(1)), F::from(5));
        assert_eq!(linear_poly.evaluate(F::from(2)), F::from(7));

        // Test evaluation of quadratic polynomial: 1 + 2x + 3x²
        let quadratic_poly = UnivariatePolynomial::new(vec![F::from(1), F::from(2), F::from(3)]);
        assert_eq!(quadratic_poly.evaluate(F::from(0)), F::from(1));
        assert_eq!(quadratic_poly.evaluate(F::from(1)), F::from(6));
        assert_eq!(quadratic_poly.evaluate(F::from(2)), F::from(17)); // 1 + 2*2 + 3*4 = 17

        // Test evaluation with field arithmetic: 4x² + 2x + 5 at x=3 mod 251
        // 4*3² + 2*3 + 5 = 4*9 + 6 + 5 = 36 + 11 = 47 mod 251
        let poly = UnivariatePolynomial::new(vec![F::from(5), F::from(2), F::from(4)]);
        assert_eq!(poly.evaluate(F::from(3)), F::from(47));
    }

    #[test]
    fn test_horners_method() {
        // Test Horner's method with larger polynomial
        // p(x) = 10 + 20x + 30x² + 40x³
        // At x=5: 10 + 20*5 + 30*25 + 40*125 = 10 + 100 + 750 + 5000 = 5860
        // But with modulo 251: 5860 % 251 = 93
        let poly =
            UnivariatePolynomial::new(vec![F::from(10), F::from(20), F::from(30), F::from(40)]);

        let result: u32 = poly.evaluate(F::from(5)).into();
        assert_eq!(result, 87); // 5860 % 251 = 87

        // Another example with a different point
        // At x=7: 10 + 20*7 + 30*49 + 40*343 = 10 + 140 + 1470 + 13720 = 15340
        // With modulo 251: 15340 % 251 = 29
        let result2 = poly.evaluate(F::from(7));
        assert_eq!(result2, F::from(29));
    }

    #[test]
    fn test_zero_evaluation() {
        // Test that a zero polynomial evaluates to zero everywhere
        let zero_poly = UnivariatePolynomial::<F>::new(vec![F::zero()]);
        for i in 0..10 {
            assert_eq!(zero_poly.evaluate(F::from(i)), F::zero());
        }
    }

    #[test]
    fn test_edge_cases() {
        // Test evaluation at field elements close to the prime
        let poly = UnivariatePolynomial::new(vec![F::from(250), F::from(1)]);
        // At x=250: 250 + 250 = 500 which is 500 % 251 = 249
        assert_eq!(poly.evaluate(F::from(250)), F::from(249));

        // Test polynomial with coefficients spanning the field
        let large_coef_poly = UnivariatePolynomial::new(vec![F::from(0), F::from(1), F::from(250)]);
        // At x=2: 0 + 2 + 250*4 = 2 + 1000 = 1002 % 251 = 249
        assert_eq!(large_coef_poly.evaluate(F::from(2)), F::from(249));
    }
}
