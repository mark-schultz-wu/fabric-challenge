//! An implementation of Montgomery arithmetic for a Prime Field of at most 31 bits.
//! Not currently constant time, though it would not be particularly difficult to make it constant time.

use crate::ff::traits::Field;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use std::convert::From;
use std::ops::Div;

/// Computes -a^(-1) mod 2^32 using the Newton-Raphson method
pub const fn mont_neg_inv(a: u32) -> u32 {
    debug_assert!(a % 2 == 1); // Must be odd
    let mut inv = 1u32;
    // 5 iterations is enough for 32 bits
    let mut i = 0;
    while i < 5 {
        inv = inv.wrapping_mul(2u32.wrapping_sub(a.wrapping_mul(inv)));
        i += 1;
    }
    inv.wrapping_neg()
}

/// A prime finite field relative to a prime p <= 2^31, implemented with Montgomery arithmetic
///
/// Note: Not currently fully constant-time, though this would not be hard to achieve.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MontgomeryFp<const P: u32>(u32);

impl<const P: u32> Field for MontgomeryFp<P> {
    fn zero() -> Self {
        Self(0)
    }
    fn one() -> Self {
        Self(Self::R)
    }
    /// Exponentiation by squaring
    fn pow(&self, exp: u32) -> Self {
        let mut result = Self::one();
        let mut base = *self;
        let mut e = exp;

        while e > 0 {
            if e & 1 == 1 {
                result *= base;
            }
            base *= base;
            e >>= 1;
        }

        result
    }
    /// Checking if zero
    fn has_no_terms(&self) -> bool {
        self.0 == 0
    }
}

/// Convert from Standard Representation to Montgomery Representation
impl<const P: u32> From<u32> for MontgomeryFp<P> {
    fn from(value: u32) -> Self {
        Self(Self::montgomery_multiply(value % P, Self::R_SQUARED))
    }
}

/// Convert from Montgomery Representation to Standard Representation
impl<const P: u32> From<MontgomeryFp<P>> for u32 {
    fn from(value: MontgomeryFp<P>) -> Self {
        MontgomeryFp::<P>::montgomery_multiply(value.0, 1)
    }
}

#[allow(dead_code)]
impl<const P: u32> MontgomeryFp<P> {
    // Constants for Montgomery arithmetic
    const R: u32 = ((1u64 << 32) % (P as u64)) as u32;
    const R_SQUARED: u32 = (((1u64 << 32) % (P as u64)).pow(2) % (P as u64)) as u32;
    const N_PRIME: u32 = mont_neg_inv(P);

    /// Constant-time conditional subtraction
    /// is `a` if condition = false, otherwise is `a - b`
    const fn const_sub(a: u32, b: u32, condition: bool) -> u32 {
        let mask = if condition { u32::MAX } else { 0 };
        a.wrapping_sub(b & mask)
    }

    /// Montgomery multiplication
    const fn montgomery_multiply(a: u32, b: u32) -> u32 {
        let ab = (a as u64) * (b as u64);
        let m = (ab as u32).wrapping_mul(Self::N_PRIME);
        let t = ((ab + (m as u64 * P as u64)) >> 32) as u32;
        Self::const_sub(t, P, t >= P)
    }

    /// Constructor
    pub const fn new(value: u32) -> Self {
        // Can't use from_u32 since it's not const-compatible
        let v = value % P;
        let r_squared = ((1u64 << 32) % (P as u64)).pow(2) as u32 % P;
        Self(Self::montgomery_multiply(v, r_squared))
    }

    /// Check if element is zero
    ///
    /// TODO: Change to be constant-time, and change each callsite.
    pub const fn has_no_terms(&self) -> bool {
        self.0 == 0
    }

    /// Multiplicative inverse using Fermat's Little Theorem
    pub fn inv(&self) -> Option<Self> {
        // TODO: change to be constant time
        if self.has_no_terms() {
            None
        } else {
            // x^(p-2) ≡ x^(-1) mod p when p is prime
            Some(self.pow(P - 2))
        }
    }

    /// Negation in-place
    pub fn neg_inplace(&mut self) {
        self.0 = MontgomeryFp::<P>::const_sub(P, self.0, !self.has_no_terms());
    }
}

// Add implementations
impl<'a, const P: u32> AddAssign<&'a MontgomeryFp<P>> for MontgomeryFp<P> {
    fn add_assign(&mut self, other: &'a Self) {
        let sum = self.0 + other.0;
        self.0 = Self::const_sub(sum, P, sum >= P);
    }
}

impl<const P: u32> AddAssign for MontgomeryFp<P> {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<const P: u32> Add for MontgomeryFp<P> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut result = self;
        result += &other;
        result
    }
}

impl<'a, const P: u32> Add<&'a MontgomeryFp<P>> for MontgomeryFp<P> {
    type Output = Self;

    fn add(self, other: &'a Self) -> Self {
        let mut result = self;
        result += other;
        result
    }
}

impl<const P: u32> Neg for &mut MontgomeryFp<P> {
    type Output = Self;
    fn neg(self) -> Self::Output {
        self.neg_inplace();
        self
    }
}

#[allow(clippy::suspicious_op_assign_impl)]
impl<'a, const P: u32> SubAssign<&'a MontgomeryFp<P>> for MontgomeryFp<P> {
    fn sub_assign(&mut self, other: &'a Self) {
        // Use negation to implement subtraction
        // a - b = a + (-b) = -((-a) + b).
        // Two negations so we do not consume b
        self.neg_inplace();
        *self += other;
        self.neg_inplace();
    }
}

impl<const P: u32> SubAssign for MontgomeryFp<P> {
    fn sub_assign(&mut self, other: Self) {
        *self -= &other;
    }
}

impl<'a, const P: u32> Sub<&'a MontgomeryFp<P>> for MontgomeryFp<P> {
    type Output = Self;
    fn sub(self, rhs: &'a Self) -> Self::Output {
        let mut output = self;
        output -= rhs;
        output
    }
}

impl<const P: u32> Sub for MontgomeryFp<P> {
    type Output = Self;
    #[allow(clippy::op_ref)]
    fn sub(self, rhs: Self) -> Self::Output {
        self - &rhs
    }
}

impl<'a, const P: u32> MulAssign<&'a MontgomeryFp<P>> for MontgomeryFp<P> {
    fn mul_assign(&mut self, other: &'a Self) {
        self.0 = Self::montgomery_multiply(self.0, other.0);
    }
}

impl<const P: u32> MulAssign for MontgomeryFp<P> {
    fn mul_assign(&mut self, other: Self) {
        *self *= &other;
    }
}

impl<'a, const P: u32> Mul<&'a MontgomeryFp<P>> for MontgomeryFp<P> {
    type Output = Self;

    fn mul(self, other: &'a Self) -> Self {
        let mut output = self;
        output *= other;
        output
    }
}

impl<const P: u32> Mul for MontgomeryFp<P> {
    type Output = Self;

    #[allow(clippy::op_ref)]
    fn mul(self, other: Self) -> Self {
        self * &other
    }
}

impl<const P: u32> Neg for MontgomeryFp<P> {
    type Output = Self;

    fn neg(self) -> Self {
        Self(Self::const_sub(P, self.0, !self.has_no_terms()))
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<const P: u32> Div for MontgomeryFp<P> {
    type Output = Option<Self>;
    fn div(self, rhs: Self) -> Self::Output {
        Some(self * rhs.inv()?)
    }
}

#[cfg(test)]
mod exhaustive_tests {
    use super::*;

    // Prime constants
    const SMALL_PRIME: u32 = 13; // 4-bit prime
    const MEDIUM_PRIME: u32 = 251; // 8-bit prime (largest 8-bit prime)
    const LARGE_PRIME: u32 = 4093; // 12-bit prime (2^12 - 3)

    // Exhaustive test functions
    fn exhaustive_test_add<const P: u32>() {
        for x in 0..P {
            for y in 0..P {
                let x_fp = MontgomeryFp::<P>::from(x);
                let y_fp = MontgomeryFp::<P>::from(y);
                let corr_x_plus_y = (x + y) % P;
                let our_x_plus_y = (x_fp + y_fp).into();
                assert_eq!(corr_x_plus_y, our_x_plus_y);
            }
        }
    }

    fn exhaustive_test_sub<const P: u32>() {
        for x in 0..P {
            for y in 0..P {
                let x_fp = MontgomeryFp::<P>::from(x);
                let y_fp = MontgomeryFp::<P>::from(y);
                let corr_x_minus_y = (x + (P - y)) % P;
                let our_x_minus_y = (x_fp - y_fp).into();
                assert_eq!(corr_x_minus_y, our_x_minus_y);
            }
        }
    }

    fn exhaustive_test_neg<const P: u32>() {
        for x in 0..P {
            let x_fp = MontgomeryFp::<P>::from(x);
            let corr_x_neg = (P - x) % P;
            let our_x_neg = x_fp.neg().into();
            assert_eq!(corr_x_neg, our_x_neg);
        }
    }

    fn exhaustive_test_mul<const P: u32>() {
        for x in 0..P {
            for y in 0..P {
                let x_fp = MontgomeryFp::<P>::from(x);
                let y_fp = MontgomeryFp::<P>::from(y);
                let corr_x_mul_y = (((x as u64) * (y as u64)) % (P as u64)) as u32;
                let our_x_mul_y = (x_fp * y_fp).into();
                assert_eq!(corr_x_mul_y, our_x_mul_y);
            }
        }
    }

    fn exhaustive_test_inv<const P: u32>() {
        assert_eq!(MontgomeryFp::<P>::from(0).inv(), None);
        for x in 1..P {
            let x_fp = MontgomeryFp::<P>::from(x);
            let x_fp_inv = x_fp.inv().unwrap();
            let one: u32 = (x_fp * x_fp_inv).into();
            assert_eq!(one, 1);
        }
    }

    fn exhaustive_test_div<const P: u32>() {
        for x in 0..P {
            for y in 1..P {
                let x_fp = MontgomeryFp::<P>::from(x);
                let y_fp = MontgomeryFp::<P>::from(y);
                let our_x_div_y = (x_fp / y_fp).unwrap();
                let possible_x_fp = our_x_div_y * y_fp;
                assert_eq!(possible_x_fp, x_fp);
            }
        }
    }

    // Test for Montgomery constants
    fn test_montgomery_constant<const P: u32>() {
        // Verify R is correctly calculated as 2^32 mod P
        let r_expected = ((1u64 << 32) % (P as u64)) as u32;
        assert_eq!(
            MontgomeryFp::<P>::R,
            r_expected,
            "Incorrect R value for prime {}",
            P
        );

        // Verify one() returns correct Montgomery form
        let one: u32 = MontgomeryFp::<P>::one().into();
        assert_eq!(one, 1, "one() doesn't convert back to 1 for prime {}", P);
    }

    // Small prime tests
    #[test]
    fn test_small_prime_all() {
        test_montgomery_constant::<SMALL_PRIME>();
        exhaustive_test_add::<SMALL_PRIME>();
        exhaustive_test_sub::<SMALL_PRIME>();
        exhaustive_test_neg::<SMALL_PRIME>();
        exhaustive_test_mul::<SMALL_PRIME>();
        exhaustive_test_inv::<SMALL_PRIME>();
        exhaustive_test_div::<SMALL_PRIME>();
    }

    // Medium prime tests (release mode only)
    #[test]
    #[cfg(not(debug_assertions))]
    fn test_medium_prime_all() {
        test_montgomery_constant::<MEDIUM_PRIME>();
        exhaustive_test_add::<MEDIUM_PRIME>();
        exhaustive_test_sub::<MEDIUM_PRIME>();
        exhaustive_test_neg::<MEDIUM_PRIME>();
        exhaustive_test_mul::<MEDIUM_PRIME>();
        exhaustive_test_inv::<MEDIUM_PRIME>();
        exhaustive_test_div::<MEDIUM_PRIME>();
    }

    // Large prime tests (release mode only)
    #[test]
    #[cfg(not(debug_assertions))]
    fn test_large_prime_all() {
        test_montgomery_constant::<LARGE_PRIME>();
        exhaustive_test_add::<LARGE_PRIME>();
        exhaustive_test_sub::<LARGE_PRIME>();
        exhaustive_test_neg::<LARGE_PRIME>();
        exhaustive_test_mul::<LARGE_PRIME>();
        exhaustive_test_inv::<LARGE_PRIME>();
        exhaustive_test_div::<LARGE_PRIME>();
    }

    // Special cases test
    fn test_div_special_cases<const P: u32>() {
        // Test 0/y = 0
        for y in 1..P {
            let zero = MontgomeryFp::<P>::from(0);
            let y_fp = MontgomeryFp::<P>::from(y);
            let result: u32 = (zero / y_fp).unwrap().into();
            assert_eq!(result, 0);
        }

        // Test x/1 = x
        for x in 0..P {
            let x_fp = MontgomeryFp::<P>::from(x);
            let one = MontgomeryFp::<P>::from(1);
            let result: u32 = (x_fp / one).unwrap().into();
            assert_eq!(result, x);
        }

        // Test div by zero returns None
        let one = MontgomeryFp::<P>::from(1);
        let zero = MontgomeryFp::<P>::from(0);
        assert_eq!(one / zero, None);
    }
    #[test]
    fn test_div_special_cases_all_primes() {
        test_div_special_cases::<SMALL_PRIME>();
        test_div_special_cases::<MEDIUM_PRIME>();
        test_div_special_cases::<LARGE_PRIME>();
    }

    // Can't exhaustively test, so solely edge cases
    #[test]
    fn test_30bit_prime_specific_values() {
        // 30-bit prime (2^30 - 35)
        const PRIME_30BIT: u32 = 1_073_741_789;

        // Test specific values and edge cases
        let test_values = [
            0,
            1,
            2,
            3,
            1_000_000,
            PRIME_30BIT / 2,
            PRIME_30BIT - 2,
            PRIME_30BIT - 1,
        ];

        for &x in &test_values {
            for &y in &test_values {
                if y == 0 {
                    continue;
                } // Skip division by zero

                // Test that a/b * b = a
                let a = MontgomeryFp::<PRIME_30BIT>::from(x);
                let b = MontgomeryFp::<PRIME_30BIT>::from(y);

                if let Some(div_result) = a / b {
                    assert_eq!(div_result * b, x.into());
                }
            }
        }
    }
}

#[cfg(test)]
mod large_prime_tests {
    use super::*;

    // 30-bit prime (2^30 - 35)
    const PRIME_30BIT: u32 = 1_073_741_789;
    const TEST_VALUES: [u32; 8] = [
        0,
        1,
        2,
        3,
        1_000_000,
        PRIME_30BIT / 2,
        PRIME_30BIT - 2,
        PRIME_30BIT - 1,
    ];

    // Test specific values instead of exhaustive testing
    #[test]
    fn targeted_test_add() {
        for &x in &TEST_VALUES {
            for &y in &TEST_VALUES {
                let x_fp = MontgomeryFp::<PRIME_30BIT>::from(x);
                let y_fp = MontgomeryFp::<PRIME_30BIT>::from(y);
                let corr_x_plus_y = (x + y) % PRIME_30BIT;
                let our_x_plus_y = (x_fp + y_fp).into();
                assert_eq!(corr_x_plus_y, our_x_plus_y);
            }
        }
    }

    #[test]
    fn targeted_test_sub() {
        for &x in &TEST_VALUES {
            for &y in &TEST_VALUES {
                let x_fp = MontgomeryFp::<PRIME_30BIT>::from(x);
                let y_fp = MontgomeryFp::<PRIME_30BIT>::from(y);
                let corr_x_minus_y = (x + (PRIME_30BIT - y)) % PRIME_30BIT;
                let our_x_minus_y = (x_fp - y_fp).into();
                assert_eq!(corr_x_minus_y, our_x_minus_y);
            }
        }
    }

    #[test]
    fn targeted_test_neg() {
        for &x in &TEST_VALUES {
            let x_fp = MontgomeryFp::<PRIME_30BIT>::from(x);
            let corr_x_neg = (PRIME_30BIT - x) % PRIME_30BIT;
            let our_x_neg = x_fp.neg().into();
            assert_eq!(corr_x_neg, our_x_neg);
        }
    }

    #[test]
    fn targeted_test_mul() {
        for &x in &TEST_VALUES {
            for &y in &TEST_VALUES {
                let x_fp = MontgomeryFp::<PRIME_30BIT>::from(x);
                let y_fp = MontgomeryFp::<PRIME_30BIT>::from(y);
                let corr_x_mul_y = (((x as u64) * (y as u64)) % (PRIME_30BIT as u64)) as u32;
                let our_x_mul_y = (x_fp * y_fp).into();
                assert_eq!(corr_x_mul_y, our_x_mul_y);
            }
        }
    }

    #[test]
    fn targeted_test_inv() {
        assert_eq!(MontgomeryFp::<PRIME_30BIT>::from(0).inv(), None);
        for &x in TEST_VALUES.iter().skip(1) {
            let x_fp = MontgomeryFp::<PRIME_30BIT>::from(x);
            let x_fp_inv = x_fp.inv().unwrap();
            let one = x_fp * x_fp_inv;
            assert_eq!(one, 1.into());
        }
    }

    #[test]
    fn targeted_test_div() {
        for &x in &TEST_VALUES {
            for &y in &TEST_VALUES {
                if y == 0 {
                    continue;
                } // Skip division by zero

                let x_fp = MontgomeryFp::<PRIME_30BIT>::from(x);
                let y_fp = MontgomeryFp::<PRIME_30BIT>::from(y);
                let our_x_div_y = (x_fp / y_fp).unwrap();
                let reconstructed_x = our_x_div_y * y_fp;
                assert_eq!(reconstructed_x, x.into());
            }
        }

        // Also verify division by zero returns None
        let one = MontgomeryFp::<PRIME_30BIT>::from(1);
        let zero = MontgomeryFp::<PRIME_30BIT>::from(0);
        assert_eq!(one / zero, None);
    }

    #[test]
    fn test_montgomery_constant_30bit() {
        // Verify R is correctly calculated
        let r_expected = ((1u64 << 32) % (PRIME_30BIT as u64)) as u32;
        assert_eq!(MontgomeryFp::<PRIME_30BIT>::R, r_expected);

        // Verify one() returns correct Montgomery form
        let one = MontgomeryFp::<PRIME_30BIT>::one();
        assert_eq!(one, 1.into());
    }
}
