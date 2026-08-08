//! Exact Gaussian-integer arithmetic and finite coprimality-topology experiments.

mod arithmetic;
mod factor;
mod topology;

pub use arithmetic::{GaussianError, GaussianInt};
pub use factor::{
    factor, is_gaussian_prime, prime_support, Factorization, PrimeClass, PrimePower,
    MAX_FACTOR_NORM,
};
pub use topology::{
    coprime, in_closure_of, sigma_contains, sigma_intersection_contains, specializes_to,
    topologically_indistinguishable, KolmogorovPoint, QuotientBasicOpen,
};
