use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MontgomeryFp<const P: u32>(u32);

#[allow(dead_code)]
impl<const P: u32> MontgomeryFp<P> {
    // Constants for Montgomery arithmetic
    const R: u32 = (1u64 << 32) as u32 % P;
    const R_SQUARED: u32 = ((1u64 << 32) % (P as u64)).pow(2) as u32 % P;
    const N_PRIME: u32 = mont_neg_inv(P);

    // Constant-time conditional subtraction
    // is `a` if condition = false, otherwise is `a - b`
    const fn const_sub(a: u32, b: u32, condition: bool) -> u32 {
        let mask = if condition { u32::MAX } else { 0 };
        a.wrapping_sub(b & mask)
    }

    // Convert from standard to Montgomery representation
    pub fn from_u32(a: u32) -> Self {
        Self(Self::montgomery_multiply(a % P, Self::R_SQUARED))
    }

    // Convert from Montgomery to standard representation
    pub fn to_u32(self) -> u32 {
        Self::montgomery_multiply(self.0, 1)
    }

    // Core Montgomery multiplication
    const fn montgomery_multiply(a: u32, b: u32) -> u32 {
        let ab = (a as u64) * (b as u64);
        let m = (ab as u32).wrapping_mul(Self::N_PRIME);
        let t = ((ab + (m as u64 * P as u64)) >> 32) as u32;
        Self::const_sub(t, P, t >= P)
    }

    // Constructor
    pub const fn new(value: u32) -> Self {
        // Can't use from_u32 since it's not const-compatible
        let v = value % P;
        let r_squared = ((1u64 << 32) % (P as u64)).pow(2) as u32 % P;
        Self(Self::montgomery_multiply(v, r_squared))
    }

    // Constants
    pub const fn zero() -> Self {
        Self(0)
    }
    pub const fn one() -> Self {
        Self(Self::R)
    }

    // Check if element is zero
    pub const fn is_zero(&self) -> bool {
        self.0 == 0
    }

    // Multiplicative inverse using Fermat's Little Theorem
    pub fn inv(&self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            // x^(p-2) ≡ x^(-1) mod p when p is prime
            Some(self.pow(P - 2))
        }
    }

    // Exponentiation by squaring
    pub fn pow(&self, exp: u32) -> Self {
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
}

// Operator implementations
impl<const P: u32> Add for MontgomeryFp<P> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        let mut output = self;
        output += other;
        output
    }
}

impl<const P: u32> AddAssign for MontgomeryFp<P> {
    fn add_assign(&mut self, other: Self) {
        let sum = self.0 + other.0;
        self.0 = Self::const_sub(sum, P, sum >= P);
    }
}

impl<const P: u32> Sub for MontgomeryFp<P> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        let mut output = self;
        output -= other;
        output
    }
}

#[allow(clippy::suspicious_op_assign_impl)]
impl<const P: u32> SubAssign for MontgomeryFp<P> {
    fn sub_assign(&mut self, other: Self) {
        // Use negation to implement subtraction
        // a - b = a + (-b)
        *self += -other;
    }
}

impl<const P: u32> Mul for MontgomeryFp<P> {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Self(Self::montgomery_multiply(self.0, other.0))
    }
}

impl<const P: u32> MulAssign for MontgomeryFp<P> {
    fn mul_assign(&mut self, other: Self) {
        self.0 = Self::montgomery_multiply(self.0, other.0);
    }
}

impl<const P: u32> Neg for MontgomeryFp<P> {
    type Output = Self;

    fn neg(self) -> Self {
        Self(Self::const_sub(P, self.0, !self.is_zero()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn exhaustive_test_add<const P: u32>() {
        for x in 0..P {
            for y in 0..P {
                let x_fp = MontgomeryFp::<P>::from_u32(x);
                let y_fp = MontgomeryFp::<P>::from_u32(y);
                let corr_x_plus_y = (x + y) % P;
                let our_x_plus_y = (x_fp + y_fp).to_u32();
                assert_eq!(corr_x_plus_y, our_x_plus_y);
            }
        }
    }
    fn exhaustive_test_sub<const P: u32>() {
        for x in 0..P {
            for y in 0..P {
                let x_fp = MontgomeryFp::<P>::from_u32(x);
                let y_fp = MontgomeryFp::<P>::from_u32(y);
                let corr_x_minus_y = (x + (P - y)) % P;
                let our_x_minus_y = (x_fp - y_fp).to_u32();
                assert_eq!(corr_x_minus_y, our_x_minus_y);
            }
        }
    }
    fn exhaustive_test_neg<const P: u32>() {
        for x in 0..P {
            let x_fp = MontgomeryFp::<P>::from_u32(x);
            let corr_x_neg = (P - x) % P;
            let our_x_neg = x_fp.neg().to_u32();
            assert_eq!(corr_x_neg, our_x_neg);
        }
    }
    fn exhaustive_test_mul<const P: u32>() {
        for x in 0..P {
            for y in 0..P {
                let x_fp = MontgomeryFp::<P>::from_u32(x);
                let y_fp = MontgomeryFp::<P>::from_u32(y);
                let corr_x_mul_y = (((x as u64) * (y as u64)) % (P as u64)) as u32;
                let our_x_mul_y = (x_fp * y_fp).to_u32();
                assert_eq!(corr_x_mul_y, our_x_mul_y);
            }
        }
    }
    // Arbitrary 14 bit prime 2^14-3 to exhaustively test
    const TEST_PRIME: u32 = 16381;
    #[test]
    fn test_14_bit_add() {
        exhaustive_test_add::<TEST_PRIME>();
    }
    #[test]
    fn test_14_bit_sub() {
        exhaustive_test_sub::<TEST_PRIME>();
    }
    #[test]
    fn test_14_bit_neg() {
        exhaustive_test_neg::<TEST_PRIME>();
    }
    #[test]
    fn test_14_bit_mul() {
        exhaustive_test_mul::<TEST_PRIME>();
    }
}
