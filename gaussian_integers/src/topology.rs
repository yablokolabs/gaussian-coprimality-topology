use std::collections::BTreeSet;

use crate::{prime_support, GaussianError, GaussianInt, PrimeClass};

/// Returns whether two Gaussian integers are coprime (their normalized gcd is a unit).
///
/// This algebraic predicate accepts zero; topology-specific functions reject zero.
///
/// # Errors
///
/// Returns a checked arithmetic error if Euclidean division is unrepresentable.
pub fn coprime(left: GaussianInt, right: GaussianInt) -> Result<bool, GaussianError> {
    Ok(left.gcd(right)?.is_unit())
}

/// Tests membership `point in sigma_index` in the paper's coprimality basis.
///
/// # Errors
///
/// Returns [`GaussianError::ZeroOutsideTopology`] if either argument is zero, or a
/// checked arithmetic error from the coprimality computation.
pub fn sigma_contains(index: GaussianInt, point: GaussianInt) -> Result<bool, GaussianError> {
    require_nonzero(index)?;
    require_nonzero(point)?;
    coprime(index, point)
}

/// Tests membership in `sigma_left_index intersect sigma_right_index` through the
/// paper's identity `sigma_a intersect sigma_b = sigma_(a*b)`.
///
/// # Errors
///
/// Returns [`GaussianError::ZeroOutsideTopology`] for a zero argument, or a checked
/// arithmetic error if the product or coprimality computation is unrepresentable.
pub fn sigma_intersection_contains(
    left_index: GaussianInt,
    right_index: GaussianInt,
    point: GaussianInt,
) -> Result<bool, GaussianError> {
    require_nonzero(left_index)?;
    require_nonzero(right_index)?;
    sigma_contains(left_index.checked_mul(right_index)?, point)
}

/// Tests topological indistinguishability by exact equality of prime supports.
///
/// # Errors
///
/// Returns the domain, finite-range, and checked-arithmetic errors of support computation.
pub fn topologically_indistinguishable(
    left: GaussianInt,
    right: GaussianInt,
) -> Result<bool, GaussianError> {
    Ok(prime_support(left)? == prime_support(right)?)
}

/// Tests the support-inclusion specialization direction `source -> target`.
///
/// This returns true exactly when `Supp(source) <= Supp(target)`, equivalently when
/// `target` belongs to the closure of the singleton `{source}`. The method name avoids
/// relying on conventions that reverse the specialization-order symbol.
///
/// # Errors
///
/// Returns the domain, finite-range, and checked-arithmetic errors of support computation.
pub fn specializes_to(source: GaussianInt, target: GaussianInt) -> Result<bool, GaussianError> {
    let source_support = prime_support(source)?;
    let target_support = prime_support(target)?;
    Ok(source_support.is_subset(&target_support))
}

/// Tests whether `point` belongs to the closure of the singleton `{source}`.
///
/// # Errors
///
/// Returns the domain, finite-range, and checked-arithmetic errors of support computation.
pub fn in_closure_of(point: GaussianInt, source: GaussianInt) -> Result<bool, GaussianError> {
    specializes_to(source, point)
}

/// A point of the finite-support model of the Kolmogorov quotient.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct KolmogorovPoint {
    support: BTreeSet<PrimeClass>,
}

impl KolmogorovPoint {
    /// Constructs a quotient point from a validated finite set of prime classes.
    #[must_use]
    pub const fn new(support: BTreeSet<PrimeClass>) -> Self {
        Self { support }
    }

    /// Maps a supported nonzero Gaussian integer to its quotient point.
    ///
    /// # Errors
    ///
    /// Returns the domain, finite-range, and checked-arithmetic errors of support computation.
    pub fn from_gaussian(value: GaussianInt) -> Result<Self, GaussianError> {
        Ok(Self::new(prime_support(value)?))
    }

    /// Returns the finite set of prime associate classes.
    #[must_use]
    pub const fn support(&self) -> &BTreeSet<PrimeClass> {
        &self.support
    }
}

/// The quotient basic open `O_F = {S : S intersect F is empty}`.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct QuotientBasicOpen {
    forbidden: BTreeSet<PrimeClass>,
}

impl QuotientBasicOpen {
    /// Constructs `O_F` from the finite forbidden support `F`.
    #[must_use]
    pub const fn new(forbidden: BTreeSet<PrimeClass>) -> Self {
        Self { forbidden }
    }

    /// Returns the finite forbidden set `F`.
    #[must_use]
    pub const fn forbidden(&self) -> &BTreeSet<PrimeClass> {
        &self.forbidden
    }

    /// Tests whether a quotient point belongs to this basic open.
    #[must_use]
    pub fn contains(&self, point: &KolmogorovPoint) -> bool {
        self.forbidden.is_disjoint(point.support())
    }
}

fn require_nonzero(value: GaussianInt) -> Result<(), GaussianError> {
    if value.is_zero() {
        Err(GaussianError::ZeroOutsideTopology)
    } else {
        Ok(())
    }
}
