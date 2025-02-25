use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

/// A finite field, supporting arbitrary field operations
#[allow(dead_code)]
pub trait Field: Sized + Add + AddAssign + Neg + Sub + SubAssign + Mul + MulAssign + Div {
    /// Additive identity of the field
    fn zero() -> Self;
    /// Multiplicative identity of the field
    fn one() -> Self;
}
