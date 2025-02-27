use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

/// A finite field, supporting arbitrary field operations
pub trait Field:
    Sized
    + Clone
    + Add<Output = Self>
    + for<'a> Add<&'a Self, Output = Self>
    + AddAssign
    + for<'a> AddAssign<&'a Self>
    + Sub<Output = Self>
    + for<'a> Sub<&'a Self, Output = Self>
    + SubAssign
    + for<'a> SubAssign<&'a Self>
    + Mul<Output = Self>
    + for<'a> Mul<&'a Self, Output = Self>
    + MulAssign
    + for<'a> MulAssign<&'a Self>
    + Div<Output = Option<Self>>
    + Neg
    + PartialEq
    + Eq
    + From<u32>
{
    /// Additive identity of the field
    fn zero() -> Self;
    /// Multiplicative identity of the field
    fn one() -> Self;
    /// Powering
    fn pow(&self, exp: u32) -> Self;
    /// Checking if zero
    fn is_zero(&self) -> bool;
}
