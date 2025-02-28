use crate::poly::univariate_poly::UnivariatePolynomial;
use crate::Field;

/// Error types for polynomial operations
#[derive(Debug, Clone, Copy)]
pub enum ShrinkError {
    NoVariablesToRemove,
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

    /// Returns a degree bound for the degree of the polynomial in the i-th variable
    fn degree(&self, variable_index: usize) -> usize;

    /// Sums the polynomial over the boolean hypercube
    fn sum_over_boolean_hypercube(&self) -> F;
}
