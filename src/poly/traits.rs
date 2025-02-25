use crate::Field;

// TODO: remove when used
#[allow(dead_code)]
/// A multivariate polynomial over a finite field.
pub trait MultivariatePolynomial<F: Field> {
    /// The concrete type used to represent univariate polynomials derived from this polynomial.
    type UnivariatePoly: UnivariatePolynomial<F>;

    /// Returns the number of variables in this polynomial.
    fn num_variables(&self) -> usize;

    /// Returns the maximum degree of any variable in the polynomial.
    /// For multilinear polynomials, this will always be 1.
    fn max_degree(&self) -> usize;

    /// Evaluates the polynomial at the given point.
    fn evaluate(&self, point: &[F]) -> F;

    /// Takes a univariate slice of the polynomial by fixing all but one variable.
    fn univariate_slice(
        &self,
        variable_index: usize,
        fixed_variables: &[F],
    ) -> Self::UnivariatePoly;

    /// Computes the sum of this polynomial over the boolean hypercube for a given variable.
    fn sum_over_boolean_hypercube(&self, variable_index: usize, partial_point: &[Option<F>]) -> F;
}

// TODO: remove when used
#[allow(dead_code)]
/// A univariate polynomial over a finite field.
pub trait UnivariatePolynomial<F: Field> {
    /// Returns the degree of this polynomial.
    /// degree of the zero polynomial is -1
    fn degree(&self) -> isize;

    /// Evaluates the polynomial at the given point.
    fn evaluate(&self, point: &F) -> F;

    /// Get the coefficients associated with the polynomial
    fn coefficients(&self) -> Vec<F>;
}
