use std::collections::BTreeSet;

use gaussian_integers::{
    factor, in_closure_of, is_gaussian_prime, prime_support, sigma_contains,
    sigma_intersection_contains, specializes_to, topologically_indistinguishable, GaussianError,
    GaussianInt, KolmogorovPoint, QuotientBasicOpen, MAX_FACTOR_NORM,
};

const fn g(re: i64, im: i64) -> GaussianInt {
    GaussianInt::new(re, im)
}

#[test]
fn gaussian_prime_classification_is_exact_in_supported_range() {
    for prime in [g(1, 1), g(2, 1), g(1, 2), g(3, 0), g(0, -3)] {
        assert!(is_gaussian_prime(prime).unwrap(), "{prime}");
    }
    for composite in [g(0, 0), g(1, 0), g(2, 0), g(5, 0), g(3, 3)] {
        assert!(!is_gaussian_prime(composite).unwrap(), "{composite}");
    }
}

#[test]
fn factorization_reconstructs_with_canonical_prime_classes() {
    for z in [g(2, 0), g(5, 0), g(-7, 11), g(0, -27), g(1, 0)] {
        let f = factor(z).unwrap();
        assert_eq!(f.reconstruct().unwrap(), z, "{f:?}");
        assert!(f.factors().windows(2).all(|w| w[0].prime() < w[1].prime()));
    }

    let two = factor(g(2, 0)).unwrap();
    assert_eq!(two.factors().len(), 1);
    assert_eq!(two.factors()[0].prime().representative(), g(1, 1));
    assert_eq!(two.factors()[0].exponent(), 2);
}

#[test]
fn property_factorization_reconstructs_a_bounded_lattice() {
    for re in -10..=10 {
        for im in -10..=10 {
            let value = g(re, im);
            if value.is_zero() {
                continue;
            }
            let result = factor(value).unwrap();
            assert_eq!(result.reconstruct().unwrap(), value);
            for power in result.factors() {
                assert!(is_gaussian_prime(power.prime().representative()).unwrap());
                assert!(power.prime().representative().divides(value).unwrap());
            }
        }
    }
}

#[test]
fn factor_and_support_reject_zero_and_out_of_range_norms() {
    assert!(matches!(
        factor(g(0, 0)),
        Err(GaussianError::ZeroOutsideTopology)
    ));
    let outside = g(31_623, 0);
    assert!(outside.norm() > MAX_FACTOR_NORM);
    assert!(matches!(
        prime_support(outside),
        Err(GaussianError::FactorizationLimitExceeded { .. })
    ));
}

#[test]
fn paper_examples_units_and_two_are_indistinguishable() {
    let empty = BTreeSet::new();
    for unit in [g(1, 0), g(-1, 0), g(0, 1), g(0, -1)] {
        assert_eq!(prime_support(unit).unwrap(), empty);
        assert!(topologically_indistinguishable(g(1, 0), unit).unwrap());
    }
    assert_eq!(
        prime_support(g(2, 0)).unwrap(),
        prime_support(g(1, 1)).unwrap()
    );
    assert!(topologically_indistinguishable(g(2, 0), g(1, 1)).unwrap());
}

#[test]
fn sigma_one_plus_i_is_exactly_mixed_parity() {
    for re in -8_i64..=8 {
        for im in -8_i64..=8 {
            if re == 0 && im == 0 {
                continue;
            }
            let expected = re.rem_euclid(2) != im.rem_euclid(2);
            assert_eq!(sigma_contains(g(1, 1), g(re, im)).unwrap(), expected);
        }
    }
}

#[test]
fn property_sigma_intersection_law_holds_on_a_bounded_domain() {
    let nonzero = (-2..=2)
        .flat_map(|re| (-2..=2).map(move |im| g(re, im)))
        .filter(|z| !z.is_zero())
        .collect::<Vec<_>>();
    for &alpha in &nonzero {
        for &beta in &nonzero {
            let product = alpha.checked_mul(beta).unwrap();
            for &point in &nonzero {
                let intersection =
                    sigma_contains(alpha, point).unwrap() && sigma_contains(beta, point).unwrap();
                assert_eq!(intersection, sigma_contains(product, point).unwrap());
                assert_eq!(
                    intersection,
                    sigma_intersection_contains(alpha, beta, point).unwrap()
                );
            }
        }
    }
}

#[test]
fn closure_and_specialization_follow_support_inclusion() {
    let x = g(1, 1);
    let y = g(3, 3);
    assert!(specializes_to(x, y).unwrap());
    assert!(in_closure_of(y, x).unwrap());
    assert!(!specializes_to(y, x).unwrap());
    assert!(!in_closure_of(x, y).unwrap());
}

#[test]
fn kolmogorov_points_and_quotient_basic_opens_use_finite_supports() {
    let point = KolmogorovPoint::from_gaussian(g(6, 0)).unwrap();
    let forbidden = KolmogorovPoint::from_gaussian(g(5, 0)).unwrap();
    let open = QuotientBasicOpen::new(forbidden.support().clone());
    assert!(open.contains(&point));

    let forbidden = KolmogorovPoint::from_gaussian(g(3, 0)).unwrap();
    let open = QuotientBasicOpen::new(forbidden.support().clone());
    assert!(!open.contains(&point));
}
