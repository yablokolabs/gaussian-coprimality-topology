use gaussian_integers::{GaussianError, GaussianInt};

const fn g(re: i64, im: i64) -> GaussianInt {
    GaussianInt::new(re, im)
}

#[test]
fn exact_checked_arithmetic_and_norm() {
    let z = g(3, -4);
    assert_eq!(z.norm(), 25);
    assert_eq!(z.checked_conjugate().unwrap(), g(3, 4));
    assert_eq!(z.checked_add(g(-2, 7)).unwrap(), g(1, 3));
    assert_eq!(z.checked_mul(g(1, 2)).unwrap(), g(11, 2));
}

#[test]
fn overflow_and_division_errors_are_structured() {
    assert!(matches!(
        g(i64::MAX, 0).checked_add(g(1, 0)),
        Err(GaussianError::ArithmeticOverflow { .. })
    ));
    assert!(matches!(
        g(1, 0).div_rem(g(0, 0)),
        Err(GaussianError::DivisionByZero)
    ));
}

#[test]
fn canonical_associate_is_shared_by_all_associates() {
    let associates = [g(-2, -3), g(3, -2), g(2, 3), g(-3, 2)];
    for z in associates {
        assert_eq!(z.canonical_associate().unwrap(), g(2, 3));
    }
    assert_eq!(g(0, -7).canonical_associate().unwrap(), g(7, 0));
    assert_eq!(
        g(2, 3).checked_associates().unwrap(),
        [g(2, 3), g(-2, -3), g(-3, 2), g(3, -2)]
    );
    assert!(g(2, 3).is_associate_of(g(-3, 2)).unwrap());
    assert!(!g(2, 3).is_associate_of(g(2, -3)).unwrap());
    assert!(GaussianInt::UNITS.iter().all(|unit| unit.is_unit()));
}

#[test]
fn euclidean_division_reconstructs_and_reduces_norm() {
    for ar in -5..=5 {
        for ai in -5..=5 {
            let a = g(ar, ai);
            for br in -3..=3 {
                for bi in -3..=3 {
                    let b = g(br, bi);
                    if b.is_zero() {
                        continue;
                    }
                    let (q, r) = a.div_rem(b).unwrap();
                    let reconstructed = b.checked_mul(q).unwrap().checked_add(r).unwrap();
                    assert_eq!(reconstructed, a, "a={a}, b={b}, q={q}, r={r}");
                    assert!(r.norm() < b.norm(), "a={a}, b={b}, q={q}, r={r}");
                }
            }
        }
    }
}

#[test]
fn exact_divisibility_and_normalized_gcd() {
    assert_eq!(g(5, 0).checked_exact_div(g(2, 1)).unwrap(), Some(g(2, -1)));
    assert_eq!(g(1, 0).checked_exact_div(g(1, 1)).unwrap(), None);
    assert!(g(2, 1).divides(g(5, 0)).unwrap());
    assert_eq!(g(5, 0).gcd(g(2, 1)).unwrap(), g(2, 1));
    assert_eq!(g(-5, 0).gcd(g(-2, -1)).unwrap(), g(2, 1));
    assert_eq!(g(0, 0).gcd(g(0, 0)).unwrap(), g(0, 0));
}

#[test]
fn property_norm_is_multiplicative_on_a_small_exhaustive_domain() {
    for ar in -4..=4 {
        for ai in -4..=4 {
            for br in -4..=4 {
                for bi in -4..=4 {
                    let a = g(ar, ai);
                    let b = g(br, bi);
                    assert_eq!(a.checked_mul(b).unwrap().norm(), a.norm() * b.norm());
                }
            }
        }
    }
}
