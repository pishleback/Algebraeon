//! Traits for operations on numbers not covered by the standard library.

pub trait DivMod<T> {
    type DivOutput;
    type ModOutput;

    /// For `a.div_mod(b)` return the unique `(q, r)` such that
    /// - `a` = `q`*`b` + `r`
    /// - 0 ≤ |`r`| < |`b`|
    /// - `r` has the same sign as `b`
    ///
    /// # Panics
    /// When `b` = 0
    ///
    /// # Examples
    /// ```
    /// use algebraeon_nzq::{Integer, traits::DivMod};
    ///
    /// let a = Integer::from(13);
    /// let b = Integer::from(5);
    /// let (q, r) = a.div_mod(b);
    /// assert_eq!(q, Integer::from(2));
    /// assert_eq!(r, Integer::from(3));
    ///
    /// let a = Integer::from(13);
    /// let b = Integer::from(-5);
    /// let (q, r) = a.div_mod(b);
    /// assert_eq!(q, Integer::from(-3));
    /// assert_eq!(r, Integer::from(-2));
    /// ```
    fn div_mod(self, other: T) -> (Self::DivOutput, Self::ModOutput);
}

pub trait Abs {
    type Output;

    /// The absolute value of `self`
    fn abs(self) -> Self::Output;
}

pub trait AbsDiff<T> {
    type Output;

    /// Computes |`self` - `other`|
    fn abs_diff(self, other: T) -> Self::Output;
}

pub trait Fraction {
    type NumeratorOutput;
    type DenominatorOutput;
    /// The numerator of `self`.
    /// - Has the same sign as `self`.
    fn numerator(self) -> Self::NumeratorOutput;
    /// The denominator of `self`.
    /// - Is always positive.
    fn denominator(self) -> Self::DenominatorOutput;
}

pub trait Floor {
    type Output;
    /// The largest integer `n` such that `n` ≤ `self`.
    fn floor(self) -> Self::Output;
}

pub trait Ceil {
    type Output;
    /// The smallest integer `n` such that `self` ≤ `n`.
    fn ceil(self) -> Self::Output;
}

pub trait ModPow<Exponent, Modulus> {
    type Output;

    /// Raise `self` the power of `exp` and reduce modulo `m`.
    /// ```
    /// use algebraeon_nzq::{Natural, traits::ModPow};
    /// assert_eq!(
    ///     Natural::from(4u32),
    ///     Natural::from(3u32).mod_pow(Natural::from(100u32), Natural::from(7u32))
    /// );
    /// ```
    fn mod_pow(self, exp: Exponent, m: Modulus) -> Self::Output;
}

pub trait ModInv<Modulus> {
    type Output;

    /// The modular inverse of `self` modulo `m`.
    /// ```
    /// use algebraeon_nzq::{Natural, traits::ModInv};
    /// assert_eq!(
    ///     Natural::from(4u32),
    ///     Natural::from(3u32).mod_inv(Natural::from(11u32)).unwrap(),
    /// );
    /// ```
    fn mod_inv(self, m: Modulus) -> Self::Output;
}
