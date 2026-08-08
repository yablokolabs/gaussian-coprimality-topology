use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::{GaussianError, GaussianInt};

/// Largest norm accepted by exact primality, factorization, and support routines.
///
/// Arithmetic and gcd operations are not subject to this limit. The factorizer uses
/// deterministic trial division, so this deliberately modest ceiling keeps runtime
/// predictable without probabilistic tests or large dependencies.
pub const MAX_FACTOR_NORM: u128 = 1_000_000_000;

/// A canonical associate class of a Gaussian prime.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct PrimeClass {
    representative: GaussianInt,
}

impl PrimeClass {
    /// Validates a Gaussian prime and returns its canonical associate class.
    ///
    /// # Errors
    ///
    /// Returns a range error above [`MAX_FACTOR_NORM`],
    /// [`GaussianError::NotGaussianPrime`] for a composite/unit/zero, or a checked
    /// arithmetic error if normalization is unrepresentable.
    pub fn new(value: GaussianInt) -> Result<Self, GaussianError> {
        if !is_gaussian_prime(value)? {
            return Err(GaussianError::NotGaussianPrime { value });
        }
        Ok(Self {
            representative: value.canonical_associate()?,
        })
    }

    /// Returns the representative with positive real and nonnegative imaginary part.
    #[must_use]
    pub const fn representative(self) -> GaussianInt {
        self.representative
    }
}

impl fmt::Display for PrimeClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}]", self.representative)
    }
}

/// One canonical prime class and its positive exponent in a factorization.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct PrimePower {
    prime: PrimeClass,
    exponent: u32,
}

impl PrimePower {
    /// Returns the prime associate class.
    #[must_use]
    pub const fn prime(self) -> PrimeClass {
        self.prime
    }

    /// Returns the positive exponent.
    #[must_use]
    pub const fn exponent(self) -> u32 {
        self.exponent
    }
}

/// A unit times sorted powers of canonical Gaussian-prime representatives.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Factorization {
    unit: GaussianInt,
    factors: Vec<PrimePower>,
}

impl Factorization {
    /// Returns the unit multiplier.
    #[must_use]
    pub const fn unit(&self) -> GaussianInt {
        self.unit
    }

    /// Returns prime powers sorted by canonical representative.
    #[must_use]
    pub fn factors(&self) -> &[PrimePower] {
        &self.factors
    }

    /// Reconstructs the original value using checked arithmetic.
    ///
    /// # Errors
    ///
    /// Returns [`GaussianError::ArithmeticOverflow`] if the represented product has
    /// unrepresentable `i64` coordinates.
    pub fn reconstruct(&self) -> Result<GaussianInt, GaussianError> {
        let mut value = self.unit;
        for factor in &self.factors {
            for _ in 0..factor.exponent {
                value = value.checked_mul(factor.prime.representative)?;
            }
        }
        Ok(value)
    }

    /// Returns the factorization's finite prime support.
    #[must_use]
    pub fn support(&self) -> BTreeSet<PrimeClass> {
        self.factors.iter().map(|power| power.prime).collect()
    }
}

/// Tests Gaussian primality exactly within [`MAX_FACTOR_NORM`].
///
/// A non-axis value is prime exactly when its norm is a rational prime. An axis
/// value is prime exactly when its nonzero coordinate is a rational prime congruent
/// to `3 mod 4`.
///
/// # Errors
///
/// Returns [`GaussianError::FactorizationLimitExceeded`] when the norm is above
/// [`MAX_FACTOR_NORM`].
pub fn is_gaussian_prime(value: GaussianInt) -> Result<bool, GaussianError> {
    ensure_supported_norm(value.norm())?;
    if value.is_zero() || value.is_unit() {
        return Ok(false);
    }
    if value.re == 0 || value.im == 0 {
        let coordinate = if value.re == 0 {
            u128::from(value.im.unsigned_abs())
        } else {
            u128::from(value.re.unsigned_abs())
        };
        return Ok(coordinate % 4 == 3 && is_rational_prime(coordinate));
    }
    Ok(is_rational_prime(value.norm()))
}

/// Factors a nonzero Gaussian integer exactly within [`MAX_FACTOR_NORM`].
///
/// The result uses canonical prime representatives and a unit, and reconstructs the
/// input exactly. Trial division is deterministic; its worst-case work is roughly
/// proportional to the number of lattice points of norm at most `sqrt(N(value))`.
///
/// # Errors
///
/// Returns [`GaussianError::ZeroOutsideTopology`] for zero,
/// [`GaussianError::FactorizationLimitExceeded`] above [`MAX_FACTOR_NORM`], or a
/// checked arithmetic error if an intermediate is unrepresentable.
pub fn factor(value: GaussianInt) -> Result<Factorization, GaussianError> {
    if value.is_zero() {
        return Err(GaussianError::ZeroOutsideTopology);
    }
    ensure_supported_norm(value.norm())?;

    let mut remainder = value;
    let mut exponents = BTreeMap::<PrimeClass, u32>::new();
    let trial_norm_limit = integer_sqrt(value.norm());

    for prime in trial_primes(trial_norm_limit)? {
        while let Some(quotient) = remainder.checked_exact_div(prime.representative())? {
            remainder = quotient;
            let exponent = exponents.entry(prime).or_default();
            *exponent = exponent
                .checked_add(1)
                .ok_or(GaussianError::ArithmeticOverflow {
                    operation: "factor exponent",
                })?;
        }
        if remainder.is_unit() {
            break;
        }
    }

    if !remainder.is_unit() {
        let prime = PrimeClass::new(remainder)?;
        let canonical = prime.representative();
        remainder =
            remainder
                .checked_exact_div(canonical)?
                .ok_or(GaussianError::ArithmeticOverflow {
                    operation: "factorization normalization",
                })?;
        let exponent = exponents.entry(prime).or_default();
        *exponent = exponent
            .checked_add(1)
            .ok_or(GaussianError::ArithmeticOverflow {
                operation: "factor exponent",
            })?;
    }

    let factors = exponents
        .into_iter()
        .map(|(prime, exponent)| PrimePower { prime, exponent })
        .collect();
    Ok(Factorization {
        unit: remainder,
        factors,
    })
}

/// Computes the finite set of canonical Gaussian-prime classes dividing `value`.
///
/// # Errors
///
/// Returns the same domain, range, and checked-arithmetic errors as [`factor`].
pub fn prime_support(value: GaussianInt) -> Result<BTreeSet<PrimeClass>, GaussianError> {
    Ok(factor(value)?.support())
}

fn ensure_supported_norm(norm: u128) -> Result<(), GaussianError> {
    if norm > MAX_FACTOR_NORM {
        Err(GaussianError::FactorizationLimitExceeded {
            norm,
            max_norm: MAX_FACTOR_NORM,
        })
    } else {
        Ok(())
    }
}

fn trial_primes(norm_limit: u128) -> Result<Vec<PrimeClass>, GaussianError> {
    let coordinate_limit =
        i64::try_from(integer_sqrt(norm_limit)).map_err(|_| GaussianError::ArithmeticOverflow {
            operation: "trial-prime bound",
        })?;
    let mut primes = Vec::new();
    for re in 1..=coordinate_limit {
        for im in 0..=coordinate_limit {
            let candidate = GaussianInt::new(re, im);
            if candidate.norm() <= norm_limit && is_gaussian_prime(candidate)? {
                primes.push(PrimeClass {
                    representative: candidate,
                });
            }
        }
    }
    primes.sort_unstable();
    Ok(primes)
}

fn is_rational_prime(value: u128) -> bool {
    if value < 2 {
        return false;
    }
    if value % 2 == 0 {
        return value == 2;
    }
    let mut divisor = 3_u128;
    while divisor <= value / divisor {
        if value % divisor == 0 {
            return false;
        }
        divisor += 2;
    }
    true
}

fn integer_sqrt(value: u128) -> u128 {
    if value < 2 {
        return value;
    }
    let mut low = 1_u128;
    let mut high = value / 2 + 1;
    while low <= high {
        let midpoint = low + (high - low) / 2;
        if midpoint <= value / midpoint {
            low = midpoint + 1;
        } else {
            high = midpoint - 1;
        }
    }
    high
}
