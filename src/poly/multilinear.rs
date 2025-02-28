use super::traits::ShrinkError;
use crate::poly::MultivariatePolynomial;
use crate::poly::UnivariatePolynomial;
use crate::Field;

/// A multilinear polynomial representation using evaluations on the boolean hypercube
#[derive(Debug, Clone)]
pub struct MultilinearPolynomial<F: Field> {
    /// Evaluations on all points of the boolean hypercube, in lexicographic order
    /// e.g.
    /// f(0,0,0), f(0,0,1), f(0,1,0), f(0,1,1), f(1,0,0), f(1,0,1), f(1,1,0), f(1,1,1)
    /// In other words, if [x]_b writes x in (fixed-length) binary, the ordering is
    /// [f([0]_b), f([1]_b), f([2]_b),...]
    evaluations: Vec<F>,

    /// Current number of variables
    current_variables: usize,
}

impl<F: Field> MultilinearPolynomial<F> {
    /// Create a new multilinear polynomial from evaluations on the boolean hypercube
    pub fn from_evaluations_on_hypercube(evaluations: Vec<F>) -> Self {
        let current_variables = evaluations.len().ilog2() as usize;
        assert_eq!(
            1 << current_variables,
            evaluations.len(),
            "Must input a set of evaluations of length 2^k"
        );
        Self {
            evaluations,
            current_variables,
        }
    }

    /// Create a zero polynomial with the specified number of variables
    pub fn zero(num_variables: usize) -> Self {
        let current_variables = num_variables;
        let evaluations = vec![F::zero(); 1 << current_variables];
        Self {
            evaluations,
            current_variables,
        }
    }

    /// Create a constant polynomial with the specified number of variables
    pub fn constant(num_variables: usize, value: F) -> Self {
        let current_variables = num_variables;
        let evaluations = vec![value; 1 << current_variables];
        Self {
            evaluations,
            current_variables,
        }
    }
}

impl<F: Field> MultivariatePolynomial<F> for MultilinearPolynomial<F> {
    /// Number of currently free variables in the Multivariate Polynomial
    fn num_variables(&self) -> usize {
        self.current_variables
    }

    /// Evaluate multilinear polynomial at a point using folding/memoization
    /// See for example Lemma 3.8 of "Proofs, Arguments, and Zero-Knowledge"
    ///
    /// The algorithm evaluates a multilinear polynomial P in variables x_1,...,x_n at point r=(r_1,...,r_n).
    /// For each j from 0 to n, we maintain a table A^j of size 2^(n-j) where:
    /// - A^0 contains P evaluated at all 2^n binary inputs
    /// - A^j contains P (with the first j variables set to r values), and remaining variables as binary.
    /// - A^n is the final evaluation at point r
    ///
    /// The recurrence relation is:
    /// A^j[b_{j+1},...,b_n] = (1-r_j) * A^{j-1}[0,b_{j+1},...,b_n] + r_j * A^{j-1}[1,b_{j+1},...,b_n]
    ///
    /// Note that the points 0,b_{j+1},...,b_n and 1,b_{j+1},...,b_n, when
    /// interpreted as integers, differ by 2^{n-j}.
    fn evaluate(&self, r: &[F]) -> F {
        if r.len() != self.current_variables {
            panic!("Incorrect number of variables provided for evaluation");
        }

        // Start with table A^0: evaluations at all binary inputs
        let mut table = self.evaluations[..(1 << self.current_variables)].to_vec();
        let mut table_size = table.len();

        // Process each variable r_j, building table A^j from A^{j-1}
        for r_j in r.iter().take(self.current_variables) {
            table_size /= 2;

            // For each binary assignment to remaining variables (b_{j+1},...,b_n)
            for idx in 0..table_size {
                // Access A^{j-1} entries with w_j=0 and w_j=1, keeping other bits the same
                let val_with_bit_0 = &table[idx];
                let val_with_bit_1 = &table[idx + table_size];

                let term0 = (F::one() - r_j) * val_with_bit_0;
                let term1 = r_j.clone() * val_with_bit_1;
                table[idx] = term0 + term1;
            }
        }

        // Final result is A^n (a single value) = P(r_1,...,r_n)
        table.swap_remove(0)
    }

    /// Extract the univariate polynomial in the last variable
    fn univariate_slice_last(&self) -> UnivariatePolynomial<F> {
        let n = self.current_variables;

        // For a multilinear polynomial, the resulting univariate is always
        // degree 1, h(x) := a + b*x. Note that
        //     * h(x) = \sum_{bi} p(b1,b2,...,bk,x), where
        //     * a := h(0), b := h(1) - h(0)
        let mut coefficients = vec![F::zero(); 2];

        // Sum over all boolean assignments to the first n-1 variables
        for prefix in 0..(1 << (n - 1)) {
            // idx_0: prefix followed by 0
            let idx_0 = prefix << 1;
            // idx_1: prefix followed by 1
            let idx_1 = idx_0 | 1;

            // a = \sum_i p(b1,...,bk,0)
            coefficients[0] += &self.evaluations[idx_0];
            // b = \sum_i p(b1,...,bk,1) - p(b1,...,bk,0)
            coefficients[1] += &self.evaluations[idx_1];
            coefficients[1] -= &self.evaluations[idx_0];
        }
        UnivariatePolynomial::new(coefficients)
    }

    /// Substitute the last variable with a field element
    ///
    /// Reduces an n-variable multilinear polynomial to an (n-1)-variable polynomial
    /// by evaluating the last variable at the given value.
    ///
    /// Uses the same technique as the evaluate function,
    /// but applies it only to the last variable.
    fn shrink_last(&mut self, value: &F) -> Result<(), ShrinkError> {
        // Check if we have any variables to remove
        if self.current_variables == 0 {
            return Err(ShrinkError::NoVariablesToRemove);
        }

        // We'll use the folding technique, but only for the last variable
        let new_size = 1 << (self.current_variables - 1);
        for idx in 0..new_size {
            // last bit 0
            let eval_with_0 = &self.evaluations[idx << 1];
            // last bit 1
            let eval_with_1 = &self.evaluations[(idx << 1) | 1];
            // computing (1 - value) * eval_with_0 + value * eval_with_1
            let mut one_minus_value = F::one();
            one_minus_value -= value;

            let term0 = one_minus_value * eval_with_0;
            let term1 = value.clone() * eval_with_1;

            self.evaluations[idx] = term0 + term1;
        }
        // Truncate the evaluations array to the new size
        self.evaluations.truncate(new_size);
        self.current_variables -= 1;
        Ok(())
    }

    /// For multilinear polynomials, the degree of any variable is at most 1.
    fn degree(&self, variable_index: usize) -> usize {
        if variable_index >= self.current_variables {
            panic!("Out of bounds variable access");
        }
        1
    }

    /// Sums the polynomial over all possible evaluations
    fn sum_over_boolean_hypercube(&self) -> F {
        self.evaluations.iter().fold(F::zero(), |acc, x| acc + x)
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::clone_on_copy)]
    use super::*;
    use crate::{Field, MontgomeryFp};
    use std::collections::HashMap;

    // Using a small prime (251) for testing
    type F = MontgomeryFp<251>;

    // Helper function to convert coefficient representation to evaluations on boolean hypercube
    fn coeffs_to_evaluations<F: Field>(coeffs: &HashMap<Vec<usize>, F>) -> Vec<F> {
        // Infer number of variables from the coefficient map
        let num_vars = if coeffs.is_empty() {
            0
        } else {
            coeffs.keys().map(|exps| exps.len()).max().unwrap_or(0)
        };

        // If no variables, just return the constant term or zero
        if num_vars == 0 {
            return vec![coeffs.get(&vec![]).cloned().unwrap_or_else(F::zero)];
        }

        let mut evals = vec![F::zero(); 1 << num_vars];

        // Iterate through all possible evaluation points in the boolean hypercube
        for point_idx in 0..(1 << num_vars) {
            let mut eval = F::zero();

            // For each monomial with its coefficient
            for (exponents, coeff) in coeffs.iter() {
                let mut padded_exponents = exponents.clone();
                // Pad the exponents vector if needed
                if padded_exponents.len() < num_vars {
                    padded_exponents.resize(num_vars, 0);
                }

                // Calculate monomial value for this point on the boolean cube
                let mut term_value = coeff.clone();

                // Check each variable's contribution
                for i in 0..num_vars {
                    if padded_exponents[i] > 0 {
                        // Variable has a non-zero exponent
                        // Check if this bit is set in point_idx
                        // The MSB corresponds to variable 0, LSB to variable num_vars-1
                        let bit_pos = num_vars - 1 - i;
                        let var_bit = (point_idx >> bit_pos) & 1;

                        if var_bit == 0 {
                            // If any variable with non-zero exponent is zero, the whole term is zero
                            term_value = F::zero();
                            break;
                        }
                        // If var_bit is 1, it doesn't matter what the exponent is on the boolean cube
                    }
                    // If exponent is 0, the variable doesn't contribute (x^0 = 1)
                }

                // Add the term's value to the evaluation
                eval += term_value;
            }

            evals[point_idx] = eval;
        }

        evals
    }

    #[test]
    fn test_from_evaluations() {
        // f(x,y) = 1 + 2x + 3y + 4xy
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![], F::from(1)); // constant term
        coeffs.insert(vec![1], F::from(2)); // x term
        coeffs.insert(vec![0, 1], F::from(3)); // y term
        coeffs.insert(vec![1, 1], F::from(4)); // xy term

        let evaluations = coeffs_to_evaluations(&coeffs);
        // Evaluations should be:
        // f(0,0) = 1
        // f(0,1) = 1 + 3 = 4
        // f(1,0) = 1 + 2 = 3
        // f(1,1) = 1 + 2 + 3 + 4 = 10
        assert_eq!(
            evaluations,
            vec![F::from(1), F::from(4), F::from(3), F::from(10)]
        );

        let poly = MultilinearPolynomial::from_evaluations_on_hypercube(evaluations);
        assert_eq!(poly.num_variables(), 2);
    }

    #[test]
    fn test_zero_and_constant() {
        // Test zero polynomial
        let zero = MultilinearPolynomial::<F>::zero(2);
        assert_eq!(zero.num_variables(), 2);
        assert_eq!(zero.evaluations.len(), 4);
        for eval in &zero.evaluations {
            assert_eq!(*eval, F::zero());
        }

        // Test constant polynomial
        let constant = MultilinearPolynomial::constant(3, F::from(5));
        assert_eq!(constant.num_variables(), 3);
        assert_eq!(constant.evaluations.len(), 8);
        for eval in &constant.evaluations {
            assert_eq!(*eval, F::from(5));
        }
    }

    #[test]
    fn test_evaluate() {
        // f(x,y) = 1 + 2x + 3y + 4xy
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![], F::from(1)); // constant term
        coeffs.insert(vec![1], F::from(2)); // x term
        coeffs.insert(vec![0, 1], F::from(3)); // y term
        coeffs.insert(vec![1, 1], F::from(4)); // xy term

        let evaluations = coeffs_to_evaluations(&coeffs);
        let poly = MultilinearPolynomial::from_evaluations_on_hypercube(evaluations);

        // Evaluate at (0,0): 1
        let result1 = poly.evaluate(&[F::from(0), F::from(0)]);
        assert_eq!(result1, F::from(1));

        // Evaluate at (1,0): 1 + 2 = 3
        let result2 = poly.evaluate(&[F::from(1), F::from(0)]);
        assert_eq!(result2, F::from(3));

        // Evaluate at (0,1): 1 + 3 = 4
        let result3 = poly.evaluate(&[F::from(0), F::from(1)]);
        assert_eq!(result3, F::from(4));

        // Evaluate at (1,1): 1 + 2 + 3 + 4 = 10
        let result4 = poly.evaluate(&[F::from(1), F::from(1)]);
        assert_eq!(result4, F::from(10));

        // Test at non-binary points - using field values that aren't just 0 or 1
        let x = F::from(123);
        let y = F::from(45);

        // Direct calculation: 1 + 2x + 3y + 4xy
        let expected = F::from(1)
            + F::from(2) * x.clone()
            + F::from(3) * y.clone()
            + F::from(4) * x.clone() * y.clone();

        let result5 = poly.evaluate(&[x, y]);
        assert_eq!(result5, expected);
    }

    #[test]
    fn test_univariate_slice_last() {
        // f(x,y) = 1 + 2x + 3y + 4xy
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![], F::from(1)); // constant term
        coeffs.insert(vec![1], F::from(2)); // x term
        coeffs.insert(vec![0, 1], F::from(3)); // y term
        coeffs.insert(vec![1, 1], F::from(4)); // xy term

        let evaluations = coeffs_to_evaluations(&coeffs);
        let poly = MultilinearPolynomial::from_evaluations_on_hypercube(evaluations);

        // The univariate slice in y should be:
        // g(y) = sum_{x in {0,1}} f(x,y)
        // g(y) = f(0,y) + f(1,y)
        // g(y) = (1 + 3y) + (1 + 2 + 3y + 4y) = 4 + 6y + 4y = 4 + 10y

        let univariate = poly.univariate_slice_last();
        assert_eq!(univariate.coefficients.len(), 2);
        assert_eq!(univariate.coefficients[0], F::from(4)); // Constant term: 1 + 1 + 2 = 4
        assert_eq!(univariate.coefficients[1], F::from(10)); // y coefficient: 3 + 3 + 4 = 10
    }

    #[test]
    fn test_shrink_last() {
        // f(x,y,z) = 1 + 2x + 3y + 4z + 5xy + 6xz + 7yz + 8xyz
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![], F::from(1)); // constant term
        coeffs.insert(vec![1], F::from(2)); // x term
        coeffs.insert(vec![0, 1], F::from(3)); // y term
        coeffs.insert(vec![0, 0, 1], F::from(4)); // z term
        coeffs.insert(vec![1, 1], F::from(5)); // xy term
        coeffs.insert(vec![1, 0, 1], F::from(6)); // xz term
        coeffs.insert(vec![0, 1, 1], F::from(7)); // yz term
        coeffs.insert(vec![1, 1, 1], F::from(8)); // xyz term

        let evaluations = coeffs_to_evaluations(&coeffs);
        let mut poly = MultilinearPolynomial::from_evaluations_on_hypercube(evaluations);

        // Substitute z = 2
        assert!(poly.shrink_last(&F::from(2)).is_ok());
        assert_eq!(poly.num_variables(), 2);

        // The resulting polynomial should be:
        // f(x,y) = 1 + 2x + 3y + 4*2 + 5xy + 6x*2 + 7y*2 + 8xy*2
        //        = 1 + 2x + 3y + 8 + 5xy + 12x + 14y + 16xy
        //        = 9 + 14x + 17y + 21xy

        // Evaluate at some test points to verify
        assert_eq!(poly.evaluate(&[F::from(0), F::from(0)]), F::from(9));
        assert_eq!(poly.evaluate(&[F::from(1), F::from(0)]), F::from(23)); // 9 + 14
        assert_eq!(poly.evaluate(&[F::from(0), F::from(1)]), F::from(26)); // 9 + 17
        assert_eq!(poly.evaluate(&[F::from(1), F::from(1)]), F::from(61)); // 9 + 14 + 17 + 21
    }

    // TODO: rewrite
    #[test]
    fn test_degree() {
        // For multilinear polynomials, the degree of each variable is always at most 1
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![], F::from(1)); // constant term
        coeffs.insert(vec![1], F::from(2)); // x term
        coeffs.insert(vec![0, 1], F::from(3)); // y term
        coeffs.insert(vec![1, 1, 1], F::from(8)); // xyz term

        let evaluations = coeffs_to_evaluations(&coeffs);
        let poly = MultilinearPolynomial::from_evaluations_on_hypercube(evaluations);

        assert_eq!(poly.degree(0), 1); // Degree of x is 1
        assert_eq!(poly.degree(1), 1); // Degree of y is 1
        assert_eq!(poly.degree(2), 1); // Degree of z is 1
    }

    #[test]
    #[should_panic(expected = "Out of bounds variable access")]
    fn test_degree_out_of_bounds() {
        let poly = MultilinearPolynomial::<F>::zero(2);
        poly.degree(2); // This should panic - only 2 variables, indices 0 and 1
    }

    #[test]
    fn test_sum_over_boolean_hypercube() {
        // f(x,y) = 1 + 2x + 3y + 4xy
        let mut coeffs = HashMap::new();
        coeffs.insert(vec![], F::from(1)); // constant term
        coeffs.insert(vec![1], F::from(2)); // x term
        coeffs.insert(vec![0, 1], F::from(3)); // y term
        coeffs.insert(vec![1, 1], F::from(4)); // xy term

        let evaluations = coeffs_to_evaluations(&coeffs);
        let poly = MultilinearPolynomial::from_evaluations_on_hypercube(evaluations);

        // Sum over boolean hypercube:
        // f(0,0) + f(0,1) + f(1,0) + f(1,1) = 1 + 4 + 3 + 10 = 18
        assert_eq!(poly.sum_over_boolean_hypercube(), F::from(18));

        // Test with a constant polynomial
        let constant = MultilinearPolynomial::constant(3, F::from(5));
        // Sum is 5 * 2^3 = 5 * 8 = 40
        assert_eq!(constant.sum_over_boolean_hypercube(), F::from(40));
    }

    #[test]
    fn test_shrink_last_error() {
        // Create a zero polynomial with 0 variables
        let mut poly = MultilinearPolynomial::<F>::zero(0);

        // Should fail since we have no variables to remove
        assert!(poly.shrink_last(&F::from(2)).is_err());
    }
}
