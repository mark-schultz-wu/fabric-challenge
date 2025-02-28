use crate::poly::univariate_poly::UnivariatePolynomial;
use crate::Field;

/// Error types for polynomial operations
#[derive(Debug, Clone, Copy)]
pub enum ShrinkError {
    NoVariablesToShrink,
}

/// Trait for multivariate polynomials
pub trait MultivariatePolynomial<F: Field>: Clone {
    /// Returns the current number of variables in this polynomial
    fn num_variables(&self) -> usize;

    /// Evaluates the polynomial at the given point
    fn evaluate(&self, point: &[F]) -> F;

    /// Returns the univariate polynomial in the last variable,
    /// with all other variables summed over their domains
    fn univariate_slice_last(&self) -> UnivariatePolynomial<F>;

    /// Substitutes the given value for the last variable,
    /// modifying `&mut self` in-place
    ///
    /// Returns `true` if the modification is successful.
    /// returns `false` if it is unsuccessful, say because
    /// `self` is already a constant polynomial.
    fn shrink_last(&mut self, value: &F) -> Result<(), ShrinkError>;

    /// Returns the degree of the polynomial in the i-th variable
    /// (The highest power of the i-th variable that appears in any term)
    /// A non-zero polynomial has degree `Some(d)`, while the zero polynomial has degree `None`.
    fn degree(&self, variable_index: usize) -> Option<usize>;

    /// Checks if the polynomial has no terms, e.g. is the empty polynomial
    fn has_no_terms(&self) -> bool;

    /// Returns the maximum degree of any single variable in the polynomial
    /// `None` is the 0 polynomial (degree -\infty)
    fn max_single_degree(&self) -> Option<usize> {
        if self.has_no_terms() {
            None
        } else {
            (0..self.num_variables())
                .filter_map(|i| self.degree(i))
                .max()
        }
    }

    /// Returns the total degree of the polynomial
    /// (The highest sum of exponents across all terms)
    /// Returns `None` for the zero polynomial
    fn total_degree(&self) -> Option<usize>;

    /// Returns true if this is a constant polynomial (no variables)
    fn has_no_variables(&self) -> bool {
        self.num_variables() == 0
    }

    /// If this is a constant polynomial, return its value
    fn constant_value(&self) -> Option<F> {
        if self.has_no_variables() {
            Some(self.evaluate(&[]))
        } else {
            None
        }
    }

    /// Sums the polynomial over the boolean hypercube
    fn sum_over_boolean_hypercube(&self) -> F;
}
