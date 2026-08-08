use std::error::Error;
use std::fmt;

/// An error from checked Gaussian-integer arithmetic or a finite experiment.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum GaussianError {
    /// A checked arithmetic operation could not be represented with `i64` coordinates.
    ArithmeticOverflow { operation: &'static str },
    /// Gaussian division was attempted with a zero divisor.
    DivisionByZero,
    /// A support/topology operation was requested for zero, which is outside the paper's space.
    ZeroOutsideTopology,
    /// Exact factorization was requested beyond the documented finite range.
    FactorizationLimitExceeded { norm: u128, max_norm: u128 },
    /// A value supplied as a Gaussian-prime class representative was not prime.
    NotGaussianPrime { value: GaussianInt },
}

impl fmt::Display for GaussianError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ArithmeticOverflow { operation } => {
                write!(f, "Gaussian integer overflow during {operation}")
            }
            Self::DivisionByZero => f.write_str("division by zero Gaussian integer"),
            Self::ZeroOutsideTopology => {
                f.write_str("zero is outside the coprimality space Z[i] \\ {0}")
            }
            Self::FactorizationLimitExceeded { norm, max_norm } => write!(
                f,
                "norm {norm} exceeds the exact factorization limit {max_norm}"
            ),
            Self::NotGaussianPrime { value } => {
                write!(f, "{value} is not a Gaussian prime in the supported range")
            }
        }
    }
}

impl Error for GaussianError {}

/// A Gaussian integer `re + im*i` with checked `i64` coordinate arithmetic.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct GaussianInt {
    /// Real coordinate.
    pub re: i64,
    /// Imaginary coordinate.
    pub im: i64,
}

impl GaussianInt {
    /// The additive identity.
    pub const ZERO: Self = Self::new(0, 0);
    /// The multiplicative identity.
    pub const ONE: Self = Self::new(1, 0);
    /// The four Gaussian units, in deterministic order.
    pub const UNITS: [Self; 4] = [
        Self::new(1, 0),
        Self::new(-1, 0),
        Self::new(0, 1),
        Self::new(0, -1),
    ];

    /// Creates `re + im*i`.
    #[must_use]
    pub const fn new(re: i64, im: i64) -> Self {
        Self { re, im }
    }

    /// Returns the exact norm `re^2 + im^2`.
    ///
    /// The norm of every pair of `i64` coordinates fits in `u128`.
    #[must_use]
    pub fn norm(self) -> u128 {
        let re = u128::from(self.re.unsigned_abs());
        let im = u128::from(self.im.unsigned_abs());
        re * re + im * im
    }

    /// Returns whether this value is zero.
    #[must_use]
    pub const fn is_zero(self) -> bool {
        self.re == 0 && self.im == 0
    }

    /// Returns whether this value is one of `1, -1, i, -i`.
    #[must_use]
    pub fn is_unit(self) -> bool {
        self.norm() == 1
    }

    /// Returns the checked additive inverse.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::ArithmeticOverflow`] if either coordinate is `i64::MIN`.
    pub fn checked_neg(self) -> Result<Self, GaussianError> {
        Ok(Self::new(
            self.re
                .checked_neg()
                .ok_or(GaussianError::ArithmeticOverflow {
                    operation: "negation",
                })?,
            self.im
                .checked_neg()
                .ok_or(GaussianError::ArithmeticOverflow {
                    operation: "negation",
                })?,
        ))
    }

    /// Returns the checked complex conjugate.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::ArithmeticOverflow`] if the imaginary coordinate is
    /// `i64::MIN`.
    pub fn checked_conjugate(self) -> Result<Self, GaussianError> {
        Ok(Self::new(
            self.re,
            self.im
                .checked_neg()
                .ok_or(GaussianError::ArithmeticOverflow {
                    operation: "conjugation",
                })?,
        ))
    }

    /// Returns the checked sum.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::ArithmeticOverflow`] if a coordinate is unrepresentable.
    pub fn checked_add(self, rhs: Self) -> Result<Self, GaussianError> {
        Ok(Self::new(
            self.re
                .checked_add(rhs.re)
                .ok_or(GaussianError::ArithmeticOverflow {
                    operation: "addition",
                })?,
            self.im
                .checked_add(rhs.im)
                .ok_or(GaussianError::ArithmeticOverflow {
                    operation: "addition",
                })?,
        ))
    }

    /// Returns the checked difference.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::ArithmeticOverflow`] if a coordinate is unrepresentable.
    pub fn checked_sub(self, rhs: Self) -> Result<Self, GaussianError> {
        Ok(Self::new(
            self.re
                .checked_sub(rhs.re)
                .ok_or(GaussianError::ArithmeticOverflow {
                    operation: "subtraction",
                })?,
            self.im
                .checked_sub(rhs.im)
                .ok_or(GaussianError::ArithmeticOverflow {
                    operation: "subtraction",
                })?,
        ))
    }

    /// Returns the checked Gaussian product.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::ArithmeticOverflow`] if an intermediate or result
    /// coordinate is unrepresentable.
    pub fn checked_mul(self, rhs: Self) -> Result<Self, GaussianError> {
        let ar = i128::from(self.re);
        let ai = i128::from(self.im);
        let br = i128::from(rhs.re);
        let bi = i128::from(rhs.im);
        let real = checked_i128_sub(ar.checked_mul(br), ai.checked_mul(bi), "multiplication")?;
        let imaginary = checked_i128_add(ar.checked_mul(bi), ai.checked_mul(br), "multiplication")?;
        Ok(Self::new(
            i64::try_from(real).map_err(|_| GaussianError::ArithmeticOverflow {
                operation: "multiplication",
            })?,
            i64::try_from(imaginary).map_err(|_| GaussianError::ArithmeticOverflow {
                operation: "multiplication",
            })?,
        ))
    }

    /// Returns the unique canonical associate in the sector `re > 0, im >= 0`.
    ///
    /// Zero is returned unchanged. An associate may be unrepresentable for coordinates
    /// involving `i64::MIN`; in that case this method returns an overflow error.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::ArithmeticOverflow`] if the canonical associate is not
    /// representable with `i64` coordinates.
    pub fn canonical_associate(self) -> Result<Self, GaussianError> {
        if self.is_zero() {
            return Ok(self);
        }
        if is_canonical(self) {
            return Ok(self);
        }

        let rotations = [
            self.im.checked_neg().map(|re| Self::new(re, self.re)),
            match (self.re.checked_neg(), self.im.checked_neg()) {
                (Some(re), Some(im)) => Some(Self::new(re, im)),
                _ => None,
            },
            self.re.checked_neg().map(|im| Self::new(self.im, im)),
        ];
        rotations
            .into_iter()
            .flatten()
            .find(|value| is_canonical(*value))
            .ok_or(GaussianError::ArithmeticOverflow {
                operation: "canonical associate normalization",
            })
    }

    /// Returns `[z, -z, i*z, -i*z]`, the four unit multiples of this value.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::ArithmeticOverflow`] if negating a coordinate is
    /// unrepresentable.
    pub fn checked_associates(self) -> Result<[Self; 4], GaussianError> {
        let negative = self.checked_neg()?;
        let i_times = Self::new(
            self.im
                .checked_neg()
                .ok_or(GaussianError::ArithmeticOverflow {
                    operation: "associate enumeration",
                })?,
            self.re,
        );
        let negative_i_times = i_times.checked_neg()?;
        Ok([self, negative, i_times, negative_i_times])
    }

    /// Returns whether two Gaussian integers differ by multiplication by a unit.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::ArithmeticOverflow`] if associate enumeration needs an
    /// unrepresentable negation.
    pub fn is_associate_of(self, other: Self) -> Result<bool, GaussianError> {
        Ok(self.checked_associates()?.contains(&other))
    }

    /// Computes Euclidean quotient and remainder, with each quotient coordinate rounded
    /// to a nearest integer (ties away from zero).
    ///
    /// For a nonzero divisor, the result satisfies `self = divisor*q + r` and
    /// `N(r) < N(divisor)`.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::DivisionByZero`] for a zero divisor, or
    /// [`GaussianError::ArithmeticOverflow`] if checked intermediates are unrepresentable.
    pub fn div_rem(self, divisor: Self) -> Result<(Self, Self), GaussianError> {
        let quotient = self.division_quotient(divisor)?;
        let remainder = self.checked_sub(divisor.checked_mul(quotient)?)?;
        Ok((quotient, remainder))
    }

    /// Returns an exact quotient, or `None` when `divisor` does not divide this value.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::DivisionByZero`] for a zero divisor, or
    /// [`GaussianError::ArithmeticOverflow`] if checked intermediates are unrepresentable.
    pub fn checked_exact_div(self, divisor: Self) -> Result<Option<Self>, GaussianError> {
        if divisor.is_zero() {
            return Err(GaussianError::DivisionByZero);
        }
        let (real, imaginary, denominator) = division_numerators(self, divisor)?;
        if real.unsigned_abs() % denominator != 0 || imaginary.unsigned_abs() % denominator != 0 {
            return Ok(None);
        }
        Ok(Some(Self::new(
            exact_ratio_to_i64(real, denominator, "exact division")?,
            exact_ratio_to_i64(imaginary, denominator, "exact division")?,
        )))
    }

    /// Returns whether `self` divides `value` exactly.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::DivisionByZero`] if `self` is zero, or
    /// [`GaussianError::ArithmeticOverflow`] if checked intermediates are unrepresentable.
    pub fn divides(self, value: Self) -> Result<bool, GaussianError> {
        Ok(value.checked_exact_div(self)?.is_some())
    }

    /// Computes a Euclidean gcd, normalized to its canonical associate.
    ///
    /// `gcd(0, 0)` is defined as zero.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::ArithmeticOverflow`] if a checked Euclidean-division
    /// intermediate or the normalized result is unrepresentable.
    pub fn gcd(self, rhs: Self) -> Result<Self, GaussianError> {
        let mut a = self;
        let mut b = rhs;
        while !b.is_zero() {
            let (_, remainder) = a.div_rem(b)?;
            a = b;
            b = remainder;
        }
        a.canonical_associate()
    }

    fn division_quotient(self, divisor: Self) -> Result<Self, GaussianError> {
        if divisor.is_zero() {
            return Err(GaussianError::DivisionByZero);
        }
        let (real, imaginary, denominator) = division_numerators(self, divisor)?;
        Ok(Self::new(
            rounded_ratio_to_i64(real, denominator)?,
            rounded_ratio_to_i64(imaginary, denominator)?,
        ))
    }
}

impl fmt::Display for GaussianInt {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (self.re, self.im) {
            (0, 0) => f.write_str("0"),
            (re, 0) => write!(f, "{re}"),
            (0, 1) => f.write_str("i"),
            (0, -1) => f.write_str("-i"),
            (0, im) => write!(f, "{im}i"),
            (re, 1) => write!(f, "{re}+i"),
            (re, -1) => write!(f, "{re}-i"),
            (re, im) if im > 0 => write!(f, "{re}+{im}i"),
            (re, im) => write!(f, "{re}{im}i"),
        }
    }
}

fn is_canonical(value: GaussianInt) -> bool {
    value.re > 0 && value.im >= 0
}

fn checked_i128_add(
    lhs: Option<i128>,
    rhs: Option<i128>,
    operation: &'static str,
) -> Result<i128, GaussianError> {
    lhs.and_then(|left| rhs.and_then(|right| left.checked_add(right)))
        .ok_or(GaussianError::ArithmeticOverflow { operation })
}

fn checked_i128_sub(
    lhs: Option<i128>,
    rhs: Option<i128>,
    operation: &'static str,
) -> Result<i128, GaussianError> {
    lhs.and_then(|left| rhs.and_then(|right| left.checked_sub(right)))
        .ok_or(GaussianError::ArithmeticOverflow { operation })
}

fn division_numerators(
    dividend: GaussianInt,
    divisor: GaussianInt,
) -> Result<(i128, i128, u128), GaussianError> {
    let ar = i128::from(dividend.re);
    let ai = i128::from(dividend.im);
    let br = i128::from(divisor.re);
    let bi = i128::from(divisor.im);
    let real = checked_i128_add(ar.checked_mul(br), ai.checked_mul(bi), "division numerator")?;
    let imaginary = checked_i128_sub(ai.checked_mul(br), ar.checked_mul(bi), "division numerator")?;
    Ok((real, imaginary, divisor.norm()))
}

fn rounded_ratio_to_i64(numerator: i128, denominator: u128) -> Result<i64, GaussianError> {
    let magnitude = numerator.unsigned_abs();
    let quotient = magnitude / denominator;
    let remainder = magnitude % denominator;
    let rounded = if remainder >= denominator - remainder {
        quotient.checked_add(1)
    } else {
        Some(quotient)
    }
    .ok_or(GaussianError::ArithmeticOverflow {
        operation: "Euclidean division",
    })?;
    signed_magnitude_to_i64(rounded, numerator.is_negative(), "Euclidean division")
}

fn exact_ratio_to_i64(
    numerator: i128,
    denominator: u128,
    operation: &'static str,
) -> Result<i64, GaussianError> {
    signed_magnitude_to_i64(
        numerator.unsigned_abs() / denominator,
        numerator.is_negative(),
        operation,
    )
}

fn signed_magnitude_to_i64(
    magnitude: u128,
    negative: bool,
    operation: &'static str,
) -> Result<i64, GaussianError> {
    if negative {
        if magnitude == (1_u128 << 63) {
            return Ok(i64::MIN);
        }
        let value = i64::try_from(magnitude)
            .map_err(|_| GaussianError::ArithmeticOverflow { operation })?;
        value
            .checked_neg()
            .ok_or(GaussianError::ArithmeticOverflow { operation })
    } else {
        i64::try_from(magnitude).map_err(|_| GaussianError::ArithmeticOverflow { operation })
    }
}
