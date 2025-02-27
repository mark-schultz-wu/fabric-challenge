use crate::poly::MultivariatePolynomial;
use crate::poly::UnivariatePolynomial;
use crate::Field;
use bitvec::prelude::*;
use std::collections::HashMap;

/// A dense multilinear polynomial representation where each variable
/// has an exponent of at most 1 in any term.
#[derive(Clone)]
pub struct DenseMultilinearPolynomial<F: Field> {
    /// Maps exponent vectors (as bitvecs) to coefficients
    /// Each bit indicates whether the corresponding variable appears (1) or not (0)
    coefficients: HashMap<BitVec, F>,

    /// Maximum number of variables in the polynomial (constant after creation)
    max_variables: usize,

    /// Current number of variables
    /// Is decremented after calling `.shrink_last`
    current_variables: usize,
}

impl<F: Field> DenseMultilinearPolynomial<F> {
    /// Create a new empty multilinear polynomial with the specified number of variables
    pub fn new(num_variables: usize) -> Self {
        Self {
            coefficients: HashMap::new(),
            max_variables: num_variables,
            current_variables: num_variables,
        }
    }

    /// Construct a polynomial from a map of exponent vectors to coefficients
    pub fn from_coefficients(mut coefficients: HashMap<BitVec, F>) -> Self {
        // First, determine the maximum number of variables
        let max_variables = coefficients.keys().map(|bv| bv.len()).max().unwrap_or(0);

        // Remove keys that have zero value
        coefficients.retain(|_, v| !v.has_no_terms());

        // Normalize all exponent vectors to have the same length
        let normalized_coefficients = coefficients
            .into_iter()
            .map(|(mut exponents, coefficient)| {
                // Ensure the exponent vector has the correct length
                if exponents.len() < max_variables {
                    exponents.resize(max_variables, false);
                }
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
    pub fn set_coefficient(&mut self, mut exponents: BitVec, coefficient: F) {
        if exponents.len() > self.max_variables {
            panic!("Exponent vector too long");
        }

        // Ensure the exponent vector has the correct length
        if exponents.len() < self.max_variables {
            exponents.resize(self.max_variables, false);
        }

        if coefficient.has_no_terms() {
            self.coefficients.remove(&exponents);
        } else {
            self.coefficients.insert(exponents, coefficient);
        }
    }

    /// Gets the coefficient for a given exponent vector, if it exists
    pub fn get_coefficient(&self, exponents: &BitVec) -> Option<&F> {
        if exponents.len() > self.max_variables {
            return None;
        }

        // Create a properly sized exponent vector if needed
        if exponents.len() < self.max_variables {
            let mut full_exponents = exponents.clone();
            full_exponents.resize(self.max_variables, false);
            self.coefficients.get(&full_exponents)
        } else {
            self.coefficients.get(exponents)
        }
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

impl<F: Field> MultivariatePolynomial<F> for DenseMultilinearPolynomial<F> {
    fn num_variables(&self) -> usize {
        self.current_variables
    }

    fn evaluate(&self, point: &[F]) -> F {
        if point.len() != self.current_variables {
            panic!("Incorrect number of variables provided for evaluation");
        }

        let mut result = F::zero();

        for (exponents, coefficient) in &self.coefficients {
            let mut term_value = coefficient.clone();

            // For each set bit in the exponent vector, multiply by the corresponding point value
            // Only consider the first current_variables elements
            for var_idx in 0..self.current_variables {
                if exponents[var_idx] {
                    term_value *= &point[var_idx];
                }
            }

            result = result + term_value;
        }

        result
    }

    fn univariate_slice_last(&self) -> UnivariatePolynomial<F> {
        if self.current_variables == 0 {
            return UnivariatePolynomial::new(vec![self.evaluate(&[])]);
        }

        // For a multilinear polynomial, we only need coefficients for x^0 and x^1
        let mut const_coeff = F::zero();
        let mut linear_coeff = F::zero();

        let last_var = self.current_variables - 1;

        for (exponents, coefficient) in &self.coefficients {
            // Skip if the exponent vector contains powers for variables beyond current_variables
            if exponents
                .iter()
                .skip(self.current_variables)
                .any(|bit| *bit)
            {
                continue;
            }

            let last_var_present = exponents[last_var];

            // Calculate how many points in the boolean hypercube contribute to this term
            let free_vars = (0..last_var).filter(|&i| !exponents[i]).count();

            // 2^free_vars is the multiplier
            let multiplier = 1 << free_vars;
            let mut contribution = F::from(multiplier as u32);
            contribution *= coefficient;

            if last_var_present {
                linear_coeff = linear_coeff + contribution;
            } else {
                const_coeff = const_coeff + contribution;
            }
        }

        UnivariatePolynomial::new(vec![const_coeff, linear_coeff])
    }

    fn shrink_last(&mut self, value: F) -> bool {
        if self.current_variables == 0 {
            return false;
        }

        let last_var = self.current_variables - 1;

        // Identify keys with the last variable set to 1
        let keys_to_process: Vec<BitVec> = self
            .coefficients
            .keys()
            .filter(|exp| exp[last_var])
            .cloned()
            .collect();

        // Process these keys
        for mut exponents in keys_to_process {
            // Remove the entry with the last variable set
            if let Some(mut coeff) = self.coefficients.remove(&exponents) {
                // Set the last variable to 0
                exponents.set(last_var, false);

                // Calculate the substituted value (for multilinear, always multiply by value once)
                coeff *= &value;

                // Only update if the new value is non-zero
                if !coeff.has_no_terms() {
                    // Update or add the entry with the substituted value
                    *self.coefficients.entry(exponents).or_insert(F::zero()) += coeff;
                }
            }
        }

        // Clean up any zero coefficients that might have resulted from additions
        self.coefficients.retain(|_, v| !v.has_no_terms());

        // Decrement the number of current variables
        self.current_variables -= 1;
        true
    }

    fn degree(&self, variable_index: usize) -> Option<usize> {
        if variable_index >= self.current_variables || self.has_no_terms() {
            return None;
        }

        // For multilinear polynomials, the degree is either 0 or 1
        let has_var = self
            .coefficients
            .iter()
            .filter(|(exponents, _)| {
                !exponents
                    .iter()
                    .skip(self.current_variables)
                    .any(|bit| *bit)
            })
            .any(|(exponents, _)| exponents[variable_index]);

        if has_var {
            Some(1)
        } else {
            Some(0)
        }
    }

    fn has_no_terms(&self) -> bool {
        self.coefficients.is_empty()
    }

    fn total_degree(&self) -> Option<usize> {
        if self.has_no_terms() {
            return None;
        }

        let mut max_degree = 0;

        for exponents in self.coefficients.keys() {
            // Skip terms with variables beyond current_variables
            let has_irrelevant_vars = exponents
                .iter()
                .skip(self.current_variables)
                .any(|bit| *bit);

            if has_irrelevant_vars {
                continue;
            }

            // Count the number of set bits (variables with exponent 1)
            let mut term_degree = 0;
            for i in 0..self.current_variables {
                if exponents[i] {
                    term_degree += 1;
                }
            }

            if term_degree > max_degree {
                max_degree = term_degree;
            }
        }

        Some(max_degree)
    }

    fn sum_over_boolean_hypercube(&self) -> F {
        let mut sum = F::zero();

        for (exponents, coefficient) in &self.coefficients {
            // Skip if the exponent vector contains powers for variables beyond current_variables
            if exponents
                .iter()
                .skip(self.current_variables)
                .any(|bit| *bit)
            {
                continue;
            }

            // Each variable with bit=0 doubles the count of contributing points
            let free_vars = (0..self.current_variables)
                .filter(|&i| !exponents[i])
                .count();

            // 2^free_vars contributing points
            let multiplier = 1 << free_vars;
            let mut contributing = F::from(multiplier as u32);
            contributing *= coefficient;
            sum += contributing;
        }

        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MontgomeryFp;

    // Using a small prime (251) for testing
    type F = MontgomeryFp<251>;

    // Helper function to create a BitVec from a list of 0s and 1s
    fn create_bitvec(values: &[usize]) -> BitVec {
        let mut bv = BitVec::new();
        for &val in values {
            match val {
                0 => bv.push(false),
                1 => bv.push(true),
                _ => panic!("Invalid value in exponent vector: only 0 and 1 are allowed for multilinear polynomials")
            }
        }
        bv
    }

    #[test]
    fn test_creation_and_getters() {
        // Create a simple multilinear polynomial: f(x,y,z) = 5 + 3xy + 2yz
        let mut coeffs = HashMap::new();
        coeffs.insert(create_bitvec(&[0, 0, 0]), F::from(5)); // Constant term
        coeffs.insert(create_bitvec(&[1, 1, 0]), F::from(3)); // xy term
        coeffs.insert(create_bitvec(&[0, 1, 1]), F::from(2)); // yz term

        let poly = DenseMultilinearPolynomial::from_coefficients(coeffs);

        // Check number of variables
        assert_eq!(poly.num_variables(), 3);

        // Check coefficients
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[0, 0, 0])).unwrap(),
            F::from(5)
        );
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[1, 1, 0])).unwrap(),
            F::from(3)
        );
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[0, 1, 1])).unwrap(),
            F::from(2)
        );

        // Check non-existent coefficients
        assert!(poly.get_coefficient(&create_bitvec(&[1, 0, 0])).is_none());
        assert!(poly.get_coefficient(&create_bitvec(&[0, 0, 1])).is_none());

        // Check oversized exponent vector
        assert!(poly
            .get_coefficient(&create_bitvec(&[1, 1, 1, 1]))
            .is_none());
    }

    #[test]
    fn test_zero_polynomial() {
        // Create an empty polynomial (should be the zero polynomial)
        let empty_coeffs: HashMap<BitVec, F> = HashMap::new();
        let zero_poly = DenseMultilinearPolynomial::from_coefficients(empty_coeffs);

        // Should have 0 variables and be recognized as zero
        assert_eq!(zero_poly.num_variables(), 0);
        assert!(zero_poly.has_no_terms());
        assert!(zero_poly.has_no_variables());

        // Create the zero polynomial explicitly with 3 variables
        let zero_poly_3vars = DenseMultilinearPolynomial::zero(3);

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
        coeffs.insert(create_bitvec(&[0, 0, 0]), F::from(5));
        coeffs.insert(create_bitvec(&[1, 0, 0]), F::zero()); // This should be removed
        coeffs.insert(create_bitvec(&[0, 1, 0]), F::from(3));

        let poly = DenseMultilinearPolynomial::from_coefficients(coeffs);

        // Should only have the non-zero coefficients
        assert!(poly.get_coefficient(&create_bitvec(&[1, 0, 0])).is_none());
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[0, 0, 0])).unwrap(),
            F::from(5)
        );
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[0, 1, 0])).unwrap(),
            F::from(3)
        );
    }

    #[test]
    fn test_polynomial_with_all_zero_coefficients() {
        // Create a polynomial where all coefficients are zero
        let mut coeffs = HashMap::new();
        coeffs.insert(create_bitvec(&[0, 0, 0]), F::zero());
        coeffs.insert(create_bitvec(&[1, 1, 0]), F::zero());
        coeffs.insert(create_bitvec(&[0, 0, 1]), F::zero());

        let poly = DenseMultilinearPolynomial::from_coefficients(coeffs);

        // Should be the zero polynomial with 3 variables
        assert_eq!(poly.num_variables(), 3);
        assert!(poly.has_no_terms());
    }

    #[test]
    fn test_set_coefficient() {
        // Create a zero polynomial with 3 variables
        let mut poly = DenseMultilinearPolynomial::zero(3);

        // Set some coefficients
        poly.set_coefficient(create_bitvec(&[1, 0, 0]), F::from(5)); // x term
        poly.set_coefficient(create_bitvec(&[0, 1, 0]), F::from(3)); // y term
        poly.set_coefficient(create_bitvec(&[0, 0, 1]), F::from(2)); // z term

        // Check coefficients were set correctly
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[1, 0, 0])).unwrap(),
            F::from(5)
        );
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[0, 1, 0])).unwrap(),
            F::from(3)
        );
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[0, 0, 1])).unwrap(),
            F::from(2)
        );

        // Update an existing coefficient
        poly.set_coefficient(create_bitvec(&[1, 0, 0]), F::from(7)); // Change x term
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[1, 0, 0])).unwrap(),
            F::from(7)
        );

        // Set a coefficient to zero should remove it
        poly.set_coefficient(create_bitvec(&[0, 1, 0]), F::zero());
        assert!(poly.get_coefficient(&create_bitvec(&[0, 1, 0])).is_none());
    }

    #[test]
    fn test_num_variables() {
        // Create multilinear polynomial: f(x,y,z) = 5 + 3xy + 2yz
        let mut coeffs = HashMap::new();
        coeffs.insert(create_bitvec(&[0, 0, 0]), F::from(5)); // Constant term
        coeffs.insert(create_bitvec(&[1, 1, 0]), F::from(3)); // xy term
        coeffs.insert(create_bitvec(&[0, 1, 1]), F::from(2)); // yz term

        let poly = DenseMultilinearPolynomial::from_coefficients(coeffs);
        assert_eq!(poly.num_variables(), 3);

        // Create multilinear polynomial with 4 variables: f(x,y,z,w) = 5 + 3xy + 2yz + w
        let mut coeffs = HashMap::new();
        coeffs.insert(create_bitvec(&[0, 0, 0, 0]), F::from(5)); // Constant term
        coeffs.insert(create_bitvec(&[1, 1, 0, 0]), F::from(3)); // xy term
        coeffs.insert(create_bitvec(&[0, 1, 1, 0]), F::from(2)); // yz term
        coeffs.insert(create_bitvec(&[0, 0, 0, 1]), F::from(1)); // w term

        let poly = DenseMultilinearPolynomial::from_coefficients(coeffs);
        assert_eq!(poly.num_variables(), 4);
    }

    #[test]
    fn test_polynomial_evaluation() {
        // Create multilinear polynomial: f(x,y,z) = 5 + 3xy + 2yz
        let mut coeffs = HashMap::new();
        coeffs.insert(create_bitvec(&[0, 0, 0]), F::from(5)); // Constant term
        coeffs.insert(create_bitvec(&[1, 1, 0]), F::from(3)); // xy term
        coeffs.insert(create_bitvec(&[0, 1, 1]), F::from(2)); // yz term

        let poly = DenseMultilinearPolynomial::from_coefficients(coeffs);

        // Evaluate at (1,1,1): 5 + 3*1*1 + 2*1*1 = 5 + 3 + 2 = 10
        let result1 = poly.evaluate(&[F::from(1), F::from(1), F::from(1)]);
        assert_eq!(result1, F::from(10));

        // Evaluate at (2,3,4): 5 + 3*2*3 + 2*3*4 = 5 + 18 + 24 = 47
        let result2 = poly.evaluate(&[F::from(2), F::from(3), F::from(4)]);
        assert_eq!(result2, F::from(47));

        // Evaluate at (0,0,0): 5 + 3*0*0 + 2*0*0 = 5
        let result3 = poly.evaluate(&[F::from(0), F::from(0), F::from(0)]);
        assert_eq!(result3, F::from(5));
    }

    #[test]
    #[should_panic(expected = "Incorrect number of variables provided for evaluation")]
    fn test_evaluation_wrong_num_vars() {
        // Create polynomial: f(x,y,z) = 5 + 3xy + 2yz
        let mut coeffs = HashMap::new();
        coeffs.insert(create_bitvec(&[0, 0, 0]), F::from(5));
        coeffs.insert(create_bitvec(&[1, 1, 0]), F::from(3));
        coeffs.insert(create_bitvec(&[0, 1, 1]), F::from(2));

        let poly = DenseMultilinearPolynomial::from_coefficients(coeffs);

        // Try to evaluate with only 2 variables
        poly.evaluate(&[F::from(1), F::from(2)]);
    }

    #[test]
    fn test_degree() {
        // Create multilinear polynomial: f(x,y,z) = 5 + 3xy + 2yz + 4x
        let mut coeffs = HashMap::new();
        coeffs.insert(create_bitvec(&[0, 0, 0]), F::from(5));
        coeffs.insert(create_bitvec(&[1, 1, 0]), F::from(3));
        coeffs.insert(create_bitvec(&[0, 1, 1]), F::from(2));
        coeffs.insert(create_bitvec(&[1, 0, 0]), F::from(4));

        let poly = DenseMultilinearPolynomial::from_coefficients(coeffs);

        // Check degrees of variables - for multilinear, degree is at most 1
        assert_eq!(poly.degree(0).unwrap(), 1); // x has max degree 1
        assert_eq!(poly.degree(1).unwrap(), 1); // y has max degree 1
        assert_eq!(poly.degree(2).unwrap(), 1); // z has max degree 1

        // Out of bounds variable
        assert_eq!(poly.degree(3), None);

        // Total degree - for multilinear, this is the maximum number of variables in any term
        assert_eq!(poly.total_degree().unwrap(), 2); // max is xy with total degree 2

        // Max single degree - for multilinear, this is always 1 if variables are present
        assert_eq!(poly.max_single_degree().unwrap(), 1);
    }

    #[test]
    fn test_univariate_slice_last() {
        // Create multilinear polynomial: f(x,y,z) = 5 + 3xy + 2yz + 4z
        let mut coeffs = HashMap::new();
        coeffs.insert(create_bitvec(&[0, 0, 0]), F::from(5)); // Constant term
        coeffs.insert(create_bitvec(&[1, 1, 0]), F::from(3)); // xy term
        coeffs.insert(create_bitvec(&[0, 1, 1]), F::from(2)); // yz term
        coeffs.insert(create_bitvec(&[0, 0, 1]), F::from(4)); // z term

        let poly = DenseMultilinearPolynomial::from_coefficients(coeffs);

        // Slice the last variable (z)
        let univariate = poly.univariate_slice_last();

        // The correct result should be: g(z) = 23 + 20z
        // This is because:
        // - The constant term 5 appears in all 4 combinations: 5*4 = 20
        // - The xy term 3 appears only when x=1, y=1: 3*1 = 3
        // - The yz term 2z appears when y=1 (and x can be 0 or 1): 2z*2 = 4z
        // - The z term 4z appears in all 4 combinations: 4z*4 = 16z
        // So g(z) = 20 + 3 + 4z + 16z = 23 + 20z

        assert_eq!(univariate.coefficients.len(), 2);
        assert_eq!(univariate.coefficients[0], F::from(23)); // Constant term: 5*4 + 3*1 = 23
        assert_eq!(univariate.coefficients[1], F::from(20)); // z coefficient: 4*4 + 2*2 = 20
    }

    #[test]
    fn test_shrink_last() {
        // Create multilinear polynomial: f(x,y,z) = 5 + 3xy + 2yz + 4z
        let mut coeffs = HashMap::new();
        coeffs.insert(create_bitvec(&[0, 0, 0]), F::from(5)); // Constant term
        coeffs.insert(create_bitvec(&[1, 1, 0]), F::from(3)); // xy term
        coeffs.insert(create_bitvec(&[0, 1, 1]), F::from(2)); // yz term
        coeffs.insert(create_bitvec(&[0, 0, 1]), F::from(4)); // z term

        let mut poly = DenseMultilinearPolynomial::from_coefficients(coeffs);

        // Substitute z = 2
        let success = poly.shrink_last(F::from(2));

        assert!(success);
        assert_eq!(poly.num_variables(), 2); // Now a polynomial in x,y

        // New polynomial should be: f(x,y) = 5 + 3xy + 2y*2 + 4*2 = 5 + 3xy + 4y + 8
        // Simplified: f(x,y) = 13 + 3xy + 4y
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[0, 0])).unwrap(),
            F::from(13)
        ); // 5 + 8 = 13
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[1, 1])).unwrap(),
            F::from(3)
        );
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[0, 1])).unwrap(),
            F::from(4)
        ); // 2*2 = 4

        // Substitute y = 3
        let success = poly.shrink_last(F::from(3));

        assert!(success);
        assert_eq!(poly.num_variables(), 1); // Now a polynomial in x

        // New polynomial should be: f(x) = 13 + 3x*3 + 4*3 = 13 + 9x + 12 = 25 + 9x
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[0])).unwrap(),
            F::from(25)
        ); // 13 + 12 = 25
        assert_eq!(
            *poly.get_coefficient(&create_bitvec(&[1])).unwrap(),
            F::from(9)
        ); // 3*3 = 9
    }

    #[test]
    fn test_no_irrelevant_ones_after_shrink() {
        // Create a simple multilinear polynomial with 3 variables
        let mut coeffs = HashMap::new();
        coeffs.insert(create_bitvec(&[0, 0, 0]), F::from(5)); // constant term
        coeffs.insert(create_bitvec(&[1, 0, 0]), F::from(2)); // x term
        coeffs.insert(create_bitvec(&[0, 1, 0]), F::from(3)); // y term
        coeffs.insert(create_bitvec(&[0, 0, 1]), F::from(4)); // z term
        coeffs.insert(create_bitvec(&[1, 1, 1]), F::from(7)); // xyz term

        let mut poly = DenseMultilinearPolynomial::from_coefficients(coeffs);
        assert_eq!(poly.num_variables(), 3);

        // Initial state should have the right number of terms
        assert_eq!(poly.coefficients.len(), 5);

        // Shrink the last variable (z)
        poly.shrink_last(F::from(2));
        assert_eq!(poly.num_variables(), 2);

        // Check that no exponent vectors have any bits set beyond the current_variables
        for exponent in poly.coefficients.keys() {
            assert_eq!(
                exponent.len(),
                3,
                "Exponent vector length should remain the same"
            );

            // Check that no 1s appear in positions at or after current_variables
            for i in poly.num_variables()..exponent.len() {
                assert!(
                    !exponent[i],
                    "No 1s should appear in positions beyond current_variables"
                );
            }
        }

        // Terms with z=1 should have been merged with their z=0 counterparts
        // So [0,0,1] merges with [0,0,0].
        // [1,1,1] gets sent to [1,1,0], which was not present, so there are still
        // the 4 terms [0,0,0], [1,0,0], [0,1,0], [1,1,0].
        assert_eq!(
            dbg!(poly.coefficients).len(),
            4,
            "After shrinking z, terms should be merged appropriately"
        );
    }

    #[test]
    fn test_shrink_last_dimensionality() {
        // Test that shrink_last decrements current_variables
        let mut coeffs = HashMap::new();
        coeffs.insert(create_bitvec(&[0, 0, 0, 0]), F::from(5)); // Constant term in 4-variable poly

        let mut poly = DenseMultilinearPolynomial::from_coefficients(coeffs);
        assert_eq!(poly.num_variables(), 4);

        // Shrink polynomial multiple times
        poly.shrink_last(F::from(2));
        assert_eq!(poly.num_variables(), 3);

        poly.shrink_last(F::from(3));
        assert_eq!(poly.num_variables(), 2);

        poly.shrink_last(F::from(4));
        assert_eq!(poly.num_variables(), 1);

        // After shrinking the last variable, we should have a constant polynomial
        assert_eq!(
            poly.get_coefficient(&create_bitvec(&[0])).unwrap().clone(),
            F::from(5)
        );
    }

    #[test]
    fn test_has_no_variables_and_has_no_terms() {
        // Zero polynomial
        let zero_poly = DenseMultilinearPolynomial::<F>::zero(0);
        assert!(zero_poly.has_no_terms());
        assert!(zero_poly.has_no_variables());

        // Constant polynomial
        let mut constant_coeffs = HashMap::new();
        constant_coeffs.insert(create_bitvec(&[0, 0, 0]), F::from(5));
        let mut constant_poly = DenseMultilinearPolynomial::from_coefficients(constant_coeffs);
        assert!(!constant_poly.has_no_terms());
        assert!(!constant_poly.has_no_variables()); // Initially not constant because num_variables == 3

        // Shrink to make it truly constant
        constant_poly.shrink_last(F::from(1));
        constant_poly.shrink_last(F::from(1));
        constant_poly.shrink_last(F::from(1));
        assert!(!constant_poly.has_no_terms());
        assert!(constant_poly.has_no_variables()); // Now it is constant

        // Non-constant polynomial
        let mut non_constant_coeffs = HashMap::new();
        non_constant_coeffs.insert(create_bitvec(&[0, 0]), F::from(5));
        non_constant_coeffs.insert(create_bitvec(&[1, 0]), F::from(3));
        let non_constant_poly = DenseMultilinearPolynomial::from_coefficients(non_constant_coeffs);
        assert!(!non_constant_poly.has_no_terms());
        assert!(!non_constant_poly.has_no_variables());
    }

    #[test]
    fn test_sum_over_boolean_hypercube() {
        // Test a constant polynomial: 5
        let mut constant_coeffs = HashMap::new();
        constant_coeffs.insert(create_bitvec(&[0, 0, 0]), F::from(5));
        let constant_poly = DenseMultilinearPolynomial::from_coefficients(constant_coeffs);
        // Sum is 5 * 2^3 = 40
        assert_eq!(constant_poly.sum_over_boolean_hypercube(), F::from(40));

        // Test f(x,y) = x + y
        let mut linear_coeffs = HashMap::new();
        linear_coeffs.insert(create_bitvec(&[1, 0]), F::from(1)); // x term
        linear_coeffs.insert(create_bitvec(&[0, 1]), F::from(1)); // y term
        let linear_poly = DenseMultilinearPolynomial::from_coefficients(linear_coeffs);
        // Sum is 1 * (0,1) + 1 * (1,0) + 2 * (1,1) = 4
        assert_eq!(linear_poly.sum_over_boolean_hypercube(), F::from(4));

        // Test f(x,y,z) = xyz (only equals 1 at (1,1,1))
        let mut product_coeffs = HashMap::new();
        product_coeffs.insert(create_bitvec(&[1, 1, 1]), F::from(1));
        let product_poly = DenseMultilinearPolynomial::from_coefficients(product_coeffs);
        // Sum is 1 at point (1,1,1), 0 elsewhere
        assert_eq!(product_poly.sum_over_boolean_hypercube(), F::from(1));
    }

    #[test]
    fn test_shrinking_beyond_dimensionality() {
        // Create a zero polynomial with 1 variable
        let mut poly = DenseMultilinearPolynomial::zero(1);

        // First shrink should succeed
        assert!(poly.shrink_last(F::from(1)));
        assert_eq!(poly.num_variables(), 0);

        // Second shrink should fail since we have no variables left
        assert!(!poly.shrink_last(F::from(2)));
        assert_eq!(poly.num_variables(), 0);
    }

    #[test]
    #[should_panic(expected = "Invalid value in exponent vector")]
    fn test_bitvec_creation_with_invalid_value() {
        // This should panic because multilinear polynomials only allow exponents of 0 or 1
        create_bitvec(&[0, 2, 0]);
    }
}
