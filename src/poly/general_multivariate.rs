use crate::poly::MultivariatePolynomial;
use crate::poly::UnivariatePolynomial;
use crate::Field;
use std::collections::HashMap;

use super::traits::ShrinkError;

/// A dense multivariate polynomial represented using a coefficient map
#[derive(Debug, Clone)]
pub struct GeneralMultivariatePolynomial<F: Field> {
    /// Maps exponent vectors to coefficients
    /// Each key is a vector of exponents, one for each variable
    /// e.g., [1, 0, 2] represents x0^1 * x2^2
    coefficients: HashMap<Vec<usize>, F>,

    /// Maximum number of variables in the polynomial (constant after creation)
    max_variables: usize,

    /// Current number of variables
    /// Is decremented after calling `.shrink_last`
    current_variables: usize,
}

impl<F: Field> GeneralMultivariatePolynomial<F> {
    /// Construct a Polynomial from a Hashmap associating exponent vectors with coefficients.
    pub fn from_coefficients(mut coefficients: HashMap<Vec<usize>, F>) -> Self {
        // First, determine the maximum number of variables
        let max_variables = coefficients
            .keys()
            .map(|exps| exps.len())
            .max()
            .unwrap_or(0);

        // Remove keys that have zero value
        coefficients.retain(|_, v| !v.is_zero());

        // Normalize all exponent vectors to have the same length
        let normalized_coefficients = coefficients
            .into_iter()
            .map(|(mut exponents, coefficient)| {
                // Extend the exponent vector if needed
                exponents.resize(max_variables, 0);
                (exponents, coefficient)
            })
            .collect();

        Self {
            coefficients: normalized_coefficients,
            max_variables,
            current_variables: max_variables,
        }
    }

    /// Sets the coefficient for a given exponent vector
    pub fn set_coefficient(&mut self, exponents: Vec<usize>, coefficient: F) {
        if exponents.len() > self.max_variables {
            panic!("Exponent vector too long");
        }

        // Extend exponent vector if needed
        let mut full_exponents = exponents;
        full_exponents.resize(self.max_variables, 0);

        if coefficient.is_zero() {
            self.coefficients.remove(&full_exponents);
        } else {
            self.coefficients.insert(full_exponents, coefficient);
        }
    }

    /// Gets the coefficient for a given exponent vector, if it exists
    pub fn get_coefficient(&self, exponents: &[usize]) -> Option<&F> {
        if exponents.len() > self.max_variables {
            return None;
        }

        // Create full exponent vector
        let mut full_exponents = exponents.to_vec();
        full_exponents.resize(self.max_variables, 0);
        self.coefficients.get(&full_exponents)
    }
    /// Creates a zero polynomial with the specified number of variables
    pub fn zero(num_variables: usize) -> Self {
        Self {
            coefficients: HashMap::new(),
            max_variables: num_variables,
            current_variables: num_variables,
        }
    }
}

impl<F: Field> MultivariatePolynomial<F> for GeneralMultivariatePolynomial<F> {
    /// Returns the current number of variables in the multivariate polynomial
    fn num_variables(&self) -> usize {
        self.current_variables
    }

    /// Evaluates the multivariate polynomial on a point
    fn evaluate(&self, point: &[F]) -> F {
        if point.len() != self.current_variables {
            panic!("Incorrect number of variables provided for evaluation");
        }

        let mut result = F::zero();

        for (exponents, coefficient) in &self.coefficients {
            let mut term = coefficient.clone();

            // Only consider the first current_variables elements
            for (var, &exp) in exponents.iter().take(self.current_variables).enumerate() {
                if exp > 0 {
                    term *= point[var].pow(exp as u32);
                }
            }
            result += term;
        }
        result
    }

    /// Maps the multivariate polynomial g(x1,...,xn) to the polynomial
    /// gi(xn) = \sum_{bi\in\{0,1\}} g(b1,b2,...,b_{n-1},xn).
    ///
    /// Note: reverse order is so we can truncate exponent vectors when responding to verifier challenges, rather than shifting them.
    fn univariate_slice_last(&self) -> UnivariatePolynomial<F> {
        if self.current_variables == 0 {
            return UnivariatePolynomial::new(vec![self.evaluate(&[])]);
        }

        let last_var = self.current_variables - 1;

        // Find the degree of the last variable
        let max_degree = self.degree(last_var).unwrap_or(0);

        // Initialize coefficient vector for the univariate polynomial
        let mut univariate_coeffs = vec![F::zero(); max_degree + 1];

        // For each term in the polynomial
        for (exponents, coefficient) in &self.coefficients {
            let last_var_power = exponents[last_var];

            // Calculate how many points in the boolean hypercube contribute to this term
            // The answer is 2^free_vars, where free_vars is the number
            // of variables with exponent == 0 (e.g. not in the term)
            // among all variables that are NOT the last variable.
            //
            // This is because each of these variables is summed over {0,1}, and
            // 1. If exponent = 0, both settings of Xi contribute to the sum
            // 2. If exponent = 1, only Xi = 1 contributes, and
            // 3. If exponent > 1, then Xi^exponent == Xi, and only Xi = 1 contributes.
            //
            // Skip the last variable as we are not summing over it when taking a univariate slice.

            let free_vars = (0..last_var).filter(|&i| exponents[i] == 0).count();
            let mut contributing = F::from(1 << free_vars as u32);
            // Multiply coefficient by number of contributing points
            contributing *= coefficient;

            // Add to the appropriate coefficient in the univariate polynomial
            univariate_coeffs[last_var_power] += contributing;
        }

        UnivariatePolynomial::new(univariate_coeffs)
    }

    /// Substitutes the `value : F` into the last variable of the multivariate polynomial,
    /// yielding a multivariate polynomial of 1 fewer degree.
    /// Returns `true` if the operation succeeded, or `false` if it failed (due to having no free variables).
    fn shrink_last(&mut self, value: &F) -> Result<(), ShrinkError> {
        if self.current_variables == 0 {
            return Err(ShrinkError::NoVariablesToShrink);
        }

        let last_var = self.current_variables - 1;

        // Identify keys with non-zero exponent for the last variable
        let keys_to_process: Vec<Vec<usize>> = self
            .coefficients
            .keys()
            .filter(|exp| exp.len() > last_var && exp[last_var] > 0)
            .cloned()
            .collect();

        // Process these keys
        for mut key in keys_to_process {
            // Remove the entry with the last variable
            if let Some(coeff) = self.coefficients.remove(&key) {
                let last_var_power = key[last_var];

                // Set the power of the last variable to 0
                key[last_var] = 0;

                // Calculate the substituted value
                let new_value = coeff * value.pow(last_var_power as u32);

                // Only update if the new value is non-zero
                if !new_value.is_zero() {
                    // Update or add the entry with the substituted value
                    *self.coefficients.entry(key).or_insert(F::zero()) += new_value;
                }
            }
        }

        // Clean up any zero coefficients that might have resulted from additions
        self.coefficients.retain(|_, v| !v.is_zero());

        // Decrement the number of current variables
        self.current_variables -= 1;
        Ok(())
    }

    /// Returns the degree of the variable of index `variable_index`, or `None`
    /// if this index is larger than the current number of variables.
    fn degree(&self, variable_index: usize) -> Option<usize> {
        if variable_index >= self.current_variables {
            panic!("Indexing an out-of-bounds variable");
        }
        // Zero polynomial has degree `None`
        if self.has_no_terms() {
            return None;
        }

        let max_degree = self
            .coefficients
            .keys()
            .filter_map(|exponents| {
                if exponents.len() > variable_index {
                    Some(exponents[variable_index])
                } else {
                    None
                }
            })
            .max();

        max_degree
    }

    /// Checks if the polynomial is the zero polynomial
    fn has_no_terms(&self) -> bool {
        self.coefficients.values().all(|v| v.is_zero())
    }

    /// Computes the total degree of the polynomial.
    /// This is the maximum of the sums of the degrees of each term.
    /// Returns `None` for the zero polynomial
    fn total_degree(&self) -> Option<usize> {
        if self.has_no_terms() {
            return None;
        }

        let max_degree = self
            .coefficients
            .keys()
            .map(|exponents| {
                // Only consider current variables in total degree calculation
                exponents.iter().take(self.current_variables).sum()
            })
            .max();

        max_degree
    }

    /// Returns `true` if `self` has no free variables, and otherwise returns `false`.
    fn has_no_variables(&self) -> bool {
        self.current_variables == 0
    }

    /// Returns the maximum degree of any single variable in the polynomial
    /// `None` for the 0 polynomial (degree -\infty)
    fn max_single_degree(&self) -> Option<usize> {
        (0..self.current_variables)
            .flat_map(|i| self.degree(i))
            .max()
    }

    fn sum_over_boolean_hypercube(&self) -> F {
        // For each term, we need to compute how many points in the
        // boolean hypercube contribute to it.
        let mut sum = F::zero();

        for (exponents, coefficient) in &self.coefficients {
            // A variable with exponent > 0 must be 1 to contribute
            // A variable with exponent 0 can be either 0 or 1

            // Each variable with exponent 0 doubles the number of contributing points
            let free_vars = exponents
                .iter()
                .take(self.current_variables)
                .filter(|&&exp| exp == 0)
                .count();

            // Each free variable can be either 0 or 1, so multiply by 2^free_vars
            let multiplier = 1 << free_vars;
            let mut contribution = F::from(multiplier as u32);
            contribution *= coefficient;
            sum += contribution;
        }
        sum
    }
}

#[cfg(test)]
mod tests {
    use crate::Field;
    use crate::GeneralMultivariatePolynomial;
    use crate::MontgomeryFp;
    use crate::MultivariatePolynomial;
    use std::collections::HashMap;

    // Using a small prime (251) for testing
    type F = MontgomeryFp<251>;

    #[test]
    fn test_creation_and_getters() {
        // Create a simple polynomial: f(x,y,z) = 5 + 3xy + 2yz^2
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::from(5)); // Constant term
        coeffs.insert(vec![1, 1, 0], F::from(3)); // xy term
        coeffs.insert(vec![0, 1, 2], F::from(2)); // yz^2 term

        let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);

        // Check number of variables
        assert_eq!(poly.num_variables(), 3);

        // Check coefficients
        assert_eq!(*poly.get_coefficient(&[0, 0, 0]).unwrap(), F::from(5));
        assert_eq!(*poly.get_coefficient(&[1, 1, 0]).unwrap(), F::from(3));
        assert_eq!(*poly.get_coefficient(&[0, 1, 2]).unwrap(), F::from(2));

        // Check non-existent coefficients
        assert!(poly.get_coefficient(&[1, 0, 0]).is_none());
        assert!(poly.get_coefficient(&[0, 0, 1]).is_none());

        // Check oversized exponent vector
        assert!(poly.get_coefficient(&[1, 1, 1, 1]).is_none());
    }

    #[test]
    fn test_zero_polynomial() {
        // Create an empty polynomial (should be the zero polynomial)
        let empty_coeffs: HashMap<Vec<usize>, F> = HashMap::new();
        let zero_poly = GeneralMultivariatePolynomial::from_coefficients(empty_coeffs);

        // Should have 0 variables and be recognized as zero
        assert_eq!(zero_poly.num_variables(), 0);
        assert!(zero_poly.has_no_terms());

        // Create the zero polynomial explicitly with 3 variables
        let zero_poly_3vars = GeneralMultivariatePolynomial::zero(3);

        // Should have 3 variables and be recognized as zero
        assert_eq!(zero_poly_3vars.num_variables(), 3);
        assert!(zero_poly_3vars.has_no_terms());

        // Evaluation should always give zero
        assert_eq!(
            zero_poly_3vars.evaluate(&[F::from(1), F::from(2), F::from(3)]),
            F::zero()
        );
    }

    #[test]
    fn test_zero_coefficient_removal() {
        // Create polynomial with some zero coefficients
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::from(5));
        coeffs.insert(vec![1, 0, 0], F::zero()); // This should be removed
        coeffs.insert(vec![0, 1, 0], F::from(3));

        let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);

        // Should only have the non-zero coefficients
        assert!(poly.get_coefficient(&[1, 0, 0]).is_none());
        assert_eq!(*poly.get_coefficient(&[0, 0, 0]).unwrap(), F::from(5));
        assert_eq!(*poly.get_coefficient(&[0, 1, 0]).unwrap(), F::from(3));
    }

    #[test]
    fn test_polynomial_with_all_zero_coefficients() {
        // Create a polynomial where all coefficients are zero
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::zero());
        coeffs.insert(vec![1, 1, 0], F::zero());
        coeffs.insert(vec![0, 0, 2], F::zero());

        let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);

        // Should be the zero polynomial with 3 variables
        assert_eq!(poly.num_variables(), 3);
        assert!(poly.has_no_terms());
    }

    #[test]
    fn test_set_coefficient() {
        // Create a zero polynomial with 3 variables
        let mut poly = GeneralMultivariatePolynomial::zero(3);

        // Set some coefficients
        poly.set_coefficient(vec![1, 0, 0], F::from(5)); // x term
        poly.set_coefficient(vec![0, 1, 0], F::from(3)); // y term
        poly.set_coefficient(vec![0, 0, 1], F::from(2)); // z term

        // Check coefficients were set correctly
        assert_eq!(*poly.get_coefficient(&[1, 0, 0]).unwrap(), F::from(5));
        assert_eq!(*poly.get_coefficient(&[0, 1, 0]).unwrap(), F::from(3));
        assert_eq!(*poly.get_coefficient(&[0, 0, 1]).unwrap(), F::from(2));

        // Update an existing coefficient
        poly.set_coefficient(vec![1, 0, 0], F::from(7)); // Change x term
        assert_eq!(*poly.get_coefficient(&[1, 0, 0]).unwrap(), F::from(7));

        // Set a coefficient to zero should remove it
        poly.set_coefficient(vec![0, 1, 0], F::zero());
        assert!(poly.get_coefficient(&[0, 1, 0]).is_none());
    }

    #[test]
    fn test_num_variables() {
        // Create polynomial: f(x,y,z) = 5 + 3xy + 2yz^2
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::from(5)); // Constant term
        coeffs.insert(vec![1, 1, 0], F::from(3)); // xy term
        coeffs.insert(vec![0, 1, 2], F::from(2)); // yz^2 term

        let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);
        assert_eq!(poly.num_variables(), 3);

        // Create polynomial with 4 variables: f(x,y,z,w) = 5 + 3xy + 2yz^2 + w
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0, 0], F::from(5)); // Constant term
        coeffs.insert(vec![1, 1, 0, 0], F::from(3)); // xy term
        coeffs.insert(vec![0, 1, 2, 0], F::from(2)); // yz^2 term
        coeffs.insert(vec![0, 0, 0, 1], F::from(1)); // w term

        let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);
        assert_eq!(poly.num_variables(), 4);
    }

    #[test]
    fn test_polynomial_evaluation() {
        // Create polynomial: f(x,y,z) = 5 + 3xy + 2yz^2
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::from(5)); // Constant term
        coeffs.insert(vec![1, 1, 0], F::from(3)); // xy term
        coeffs.insert(vec![0, 1, 2], F::from(2)); // yz^2 term

        let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);

        // Evaluate at (1,1,1): 5 + 3*1*1 + 2*1*1^2 = 5 + 3 + 2 = 10
        let result1 = poly.evaluate(&[F::from(1), F::from(1), F::from(1)]);
        assert_eq!(result1, F::from(10));

        // Evaluate at (2,3,4): 5 + 3*2*3 + 2*3*4^2 = 5 + 18 + 2*3*16 = 5 + 18 + 96 = 119
        let result2 = poly.evaluate(&[F::from(2), F::from(3), F::from(4)]);
        assert_eq!(result2, F::from(119));

        // Evaluate at (0,0,0): 5 + 3*0*0 + 2*0*0^2 = 5
        let result3 = poly.evaluate(&[F::from(0), F::from(0), F::from(0)]);
        assert_eq!(result3, F::from(5));
    }

    #[test]
    #[should_panic(expected = "Incorrect number of variables provided for evaluation")]
    fn test_evaluation_wrong_num_vars() {
        // Create polynomial: f(x,y,z) = 5 + 3xy + 2yz^2
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::from(5));
        coeffs.insert(vec![1, 1, 0], F::from(3));
        coeffs.insert(vec![0, 1, 2], F::from(2));

        let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);

        // Try to evaluate with only 2 variables
        poly.evaluate(&[F::from(1), F::from(2)]);
    }

    #[test]
    #[should_panic]
    fn test_out_of_bounds_access() {
        // Create polynomial: f(x,y,z) = 5 + 3xy + 2yz^2 + 4x^3
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::from(5));
        coeffs.insert(vec![1, 1, 0], F::from(3));
        coeffs.insert(vec![0, 1, 2], F::from(2));
        coeffs.insert(vec![3, 0, 0], F::from(4));

        let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);

        // Check degrees of variables
        assert_eq!(poly.degree(0).unwrap(), 3); // x has max degree 3
        assert_eq!(poly.degree(1).unwrap(), 1); // y has max degree 1
        assert_eq!(poly.degree(2).unwrap(), 2); // z has max degree 2

        // Out of bounds variable
        poly.degree(3);
    }

    #[test]
    fn test_degree() {
        // Create polynomial: f(x,y,z) = 5 + 3xy + 2yz^2 + 4x^3
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::from(5));
        coeffs.insert(vec![1, 1, 0], F::from(3));
        coeffs.insert(vec![0, 1, 2], F::from(2));
        coeffs.insert(vec![3, 0, 0], F::from(4));

        let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);

        // Check degrees of variables
        assert_eq!(poly.degree(0).unwrap(), 3); // x has max degree 3
        assert_eq!(poly.degree(1).unwrap(), 1); // y has max degree 1
        assert_eq!(poly.degree(2).unwrap(), 2); // z has max degree 2

        // Total degree
        assert_eq!(poly.total_degree().unwrap(), 3); // max is x^3 with total degree 3

        // Max single degree
        assert_eq!(poly.max_single_degree().unwrap(), 3); // x has max degree 3
    }

    #[test]
    fn test_univariate_slice_last() {
        // Create polynomial: f(x,y,z) = 5 + 3xy + 2yz^2 + 4z^3
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::from(5)); // Constant term
        coeffs.insert(vec![1, 1, 0], F::from(3)); // xy term
        coeffs.insert(vec![0, 1, 2], F::from(2)); // yz^2 term
        coeffs.insert(vec![0, 0, 3], F::from(4)); // z^3 term

        let poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);

        // Slice the last variable (z)
        let univariate = poly.univariate_slice_last();

        // The correct result should be: g(z) = 23 + 0z + 4z^2 + 16z^3
        // This is because when we sum over the boolean hypercube for x,y:
        // - The constant term 5 appears in all 4 combinations: 5*4 = 20
        // - The xy term 3 appears only when x=1, y=1: 3*1 = 3
        // - The yz^2 term 2z^2 appears when y=1 (and x can be 0 or 1): 2z^2*2 = 4z^2
        // - The z^3 term 4z^3 appears in all 4 combinations: 4z^3*4 = 16z^3
        // So g(z) = 20 + 3 + 4z^2 + 16z^3 = 23 + 4z^2 + 16z^3

        assert_eq!(univariate.coefficients.len(), 4);
        assert_eq!(univariate.coefficients[0], F::from(23)); // Constant term: 5*4 + 3*1 = 23
        assert_eq!(univariate.coefficients[1], F::zero()); // z coefficient: 0
        assert_eq!(univariate.coefficients[2], F::from(4)); // z^2 coefficient: 2*2 = 4
        assert_eq!(univariate.coefficients[3], F::from(16)); // z^3 coefficient: 4*4 = 16
    }

    #[test]
    fn test_shrink_last() {
        // Create polynomial: f(x,y,z) = 5 + 3xy + 2yz^2 + 4z^3
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::from(5)); // Constant term
        coeffs.insert(vec![1, 1, 0], F::from(3)); // xy term
        coeffs.insert(vec![0, 1, 2], F::from(2)); // yz^2 term
        coeffs.insert(vec![0, 0, 3], F::from(4)); // z^3 term

        let mut poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);

        // Substitute z = 2
        let success = poly.shrink_last(&F::from(2));

        assert!(success.is_ok());
        assert_eq!(poly.num_variables(), 2); // Now a polynomial in x,y
                                             // max_variables should still be 3
        assert_eq!(poly.max_variables, 3);

        // New polynomial should be: f(x,y) = 5 + 3xy + 2y*2^2 + 4*2^3 = 5 + 3xy + 8y + 32
        // Simplified: f(x,y) = 37 + 3xy + 8y
        assert_eq!(*poly.get_coefficient(&[0, 0]).unwrap(), F::from(37)); // 5 + 32 = 37
        assert_eq!(*poly.get_coefficient(&[1, 1]).unwrap(), F::from(3));
        assert_eq!(*poly.get_coefficient(&[0, 1]).unwrap(), F::from(8)); // 2*2^2 = 8

        // Substitute y = 3
        let success = poly.shrink_last(&F::from(3));

        assert!(success.is_ok());
        assert_eq!(poly.num_variables(), 1); // Now a polynomial in x
                                             // max_variables should still be 3
        assert_eq!(poly.max_variables, 3);

        // New polynomial should be: f(x) = 37 + 3x*3 + 8*3 = 37 + 9x + 24 = 61 + 9x
        assert_eq!(*poly.get_coefficient(&[0]).unwrap(), F::from(61)); // 37 + 24 = 61
        assert_eq!(*poly.get_coefficient(&[1]).unwrap(), F::from(9)); // 3*3 = 9
    }

    #[test]
    fn test_shrink_last_dimensionality() {
        // Test that shrink_last preserves max_variables while reducing current_variables
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0, 0], F::from(5)); // Constant term in 4-variable poly

        let mut poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);
        assert_eq!(poly.max_variables, 4);
        assert_eq!(poly.num_variables(), 4);

        // Shrink polynomial multiple times
        assert!(poly.shrink_last(&F::from(2)).is_ok());
        assert_eq!(poly.max_variables, 4);
        assert_eq!(poly.num_variables(), 3);

        assert!(poly.shrink_last(&F::from(3)).is_ok());
        assert_eq!(poly.max_variables, 4);
        assert_eq!(poly.num_variables(), 2);

        assert!(poly.shrink_last(&F::from(4)).is_ok());
        assert_eq!(poly.max_variables, 4);
        assert_eq!(poly.num_variables(), 1);

        // After shrinking the last variable, we should have a constant polynomial
        assert_eq!(poly.get_coefficient(&[0]).unwrap().clone(), F::from(5));
    }

    #[test]
    fn test_no_irrelevant_ones_after_shrink() {
        // Create a simple multivariate polynomial with 3 variables
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![0, 0, 0], F::from(5)); // constant term
        coeffs.insert(vec![1, 0, 0], F::from(2)); // x term
        coeffs.insert(vec![0, 1, 0], F::from(3)); // y term
        coeffs.insert(vec![0, 0, 1], F::from(4)); // z term
        coeffs.insert(vec![1, 1, 1], F::from(7)); // xyz term

        let mut poly = GeneralMultivariatePolynomial::from_coefficients(coeffs);
        assert_eq!(poly.num_variables(), 3);

        // Initial state should have the right number of terms
        assert_eq!(poly.coefficients.len(), 5);

        // Shrink the last variable (z)
        assert!(poly.shrink_last(&F::from(2)).is_ok());
        assert_eq!(poly.num_variables(), 2);

        // Check that no exponent vectors have any non-zero values beyond the current_variables
        for exponent in poly.coefficients.keys() {
            assert_eq!(
                exponent.len(),
                3,
                "Exponent vector length should remain the same"
            );

            // Check that no non-zero values appear in positions at or after current_variables
            for degree in exponent.iter().skip(poly.current_variables) {
                assert_eq!(
                    *degree, 0,
                    "No non-zero exponents should appear in positions beyond current_variables"
                );
            }
        }

        // Terms with z=1 should have been merged with their z=0 counterparts
        // So [0,0,1] merges with [0,0,0] -> [0,0,0] with value 5 + 4*2 = 13
        // [1,1,1] gets sent to [1,1,0] with value 7*2 = 14, which was not present before
        // So we have 4 terms: [0,0,0], [1,0,0], [0,1,0], [1,1,0]
        assert_eq!(
            poly.coefficients.len(),
            4,
            "After shrinking z, terms should be merged appropriately"
        );

        // Check specific coefficients
        assert_eq!(
            *poly.get_coefficient(&[0, 0, 0]).unwrap(),
            F::from(13), // 5 + 4*2 = 13
            "Constant term should be updated correctly"
        );

        assert_eq!(
            *poly.get_coefficient(&[1, 1, 0]).unwrap(),
            F::from(14), // 7*2 = 14
            "New xy term should have correct coefficient"
        );
    }

    #[test]
    fn test_has_no_variables_and_has_no_terms() {
        // Zero polynomial
        let zero_poly = GeneralMultivariatePolynomial::<F>::zero(0);
        assert!(zero_poly.has_no_terms());
        assert!(zero_poly.has_no_variables());

        // Constant polynomial
        let mut constant_coeffs = HashMap::new();
        constant_coeffs.insert(vec![0, 0, 0], F::from(5));
        let mut constant_poly = GeneralMultivariatePolynomial::from_coefficients(constant_coeffs);
        assert!(!constant_poly.has_no_terms());
        assert!(!constant_poly.has_no_variables()); // Initially not constant because num_variables == 3

        // Shrink to make it truly constant
        assert!(constant_poly.shrink_last(&F::from(1)).is_ok());
        assert!(constant_poly.shrink_last(&F::from(1)).is_ok());
        assert!(constant_poly.shrink_last(&F::from(1)).is_ok());
        assert!(!constant_poly.has_no_terms());
        assert!(constant_poly.has_no_variables()); // Now it is constant

        // Non-constant polynomial
        let mut non_constant_coeffs = HashMap::new();
        non_constant_coeffs.insert(vec![0, 0], F::from(5));
        non_constant_coeffs.insert(vec![1, 0], F::from(3));
        let non_constant_poly =
            GeneralMultivariatePolynomial::from_coefficients(non_constant_coeffs);
        assert!(!non_constant_poly.has_no_terms());
        assert!(!non_constant_poly.has_no_variables());
    }

    #[test]
    fn test_sum_over_boolean_hypercube() {
        // Test a constant polynomial: 5
        let mut constant_coeffs = HashMap::new();
        constant_coeffs.insert(vec![0, 0, 0], F::from(5));
        let constant_poly = GeneralMultivariatePolynomial::from_coefficients(constant_coeffs);
        // Sum is 5 * 2^3 = 40
        assert_eq!(constant_poly.sum_over_boolean_hypercube(), F::from(40));

        // Test f(x,y) = x + y
        let mut linear_coeffs = HashMap::new();
        linear_coeffs.insert(vec![1, 0], F::from(1)); // x term
        linear_coeffs.insert(vec![0, 1], F::from(1)); // y term
        let linear_poly = GeneralMultivariatePolynomial::from_coefficients(linear_coeffs);
        // Sum is 1 * (0,1) + 1 * (1,0) + 2 * (1,1) = 4
        assert_eq!(linear_poly.sum_over_boolean_hypercube(), F::from(4));

        // Test f(x,y,z) = xyz (only equals 1 at (1,1,1))
        let mut product_coeffs = HashMap::new();
        product_coeffs.insert(vec![1, 1, 1], F::from(1));
        let product_poly = GeneralMultivariatePolynomial::from_coefficients(product_coeffs);
        // Sum is 1 at point (1,1,1), 0 elsewhere
        assert_eq!(product_poly.sum_over_boolean_hypercube(), F::from(1));
    }

    #[test]
    fn test_shrinking_beyond_dimensionality() {
        // Create a zero polynomial with 1 variable
        let mut poly = GeneralMultivariatePolynomial::zero(1);

        // First shrink should succeed
        assert!(poly.shrink_last(&F::from(1)).is_ok());
        assert_eq!(poly.num_variables(), 0);

        // Second shrink should fail since we have no variables left
        assert!(poly.shrink_last(&F::from(2)).is_err());
        assert_eq!(poly.num_variables(), 0);
    }
}
