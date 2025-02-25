use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

/// A finite field, supporting arbitrary field operations
#[allow(dead_code)]
pub trait Field:
    Sized
    + Clone
    + Add
    + for<'a> Add<&'a Self>
    + AddAssign
    + for<'a> AddAssign<&'a Self>
    + Sub
    + for<'a> Sub<&'a Self>
    + SubAssign
    + for<'a> SubAssign<&'a Self>
    + Mul
    + for<'a> Mul<&'a Self>
    + MulAssign
    + for<'a> MulAssign<&'a Self>
    + Div
    + Neg
{
    /// Additive identity of the field
    fn zero() -> Self;
    /// Multiplicative identity of the field
    fn one() -> Self;
}
