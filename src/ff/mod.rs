/// An implementation of a Finite (Prime) Field using Montgomery Arithmetic
mod small_prime_mont;
/// A trait for defining a Finite Field
mod traits;
pub use small_prime_mont::MontgomeryFp;
pub use traits::Field;
