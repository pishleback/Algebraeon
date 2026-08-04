use crate::finite_fields::conway_finite_fields::ConwayFiniteFieldStructure;
use crate::parsing::parse_quaternion;
use crate::quaternion_algebra::{QuaternionAlgebraElement, QuaternionAlgebraStructure};
use crate::structure::*;
use algebraeon_structures::*;
use std::sync::Arc;

/// Hamilton quaternion algebra: H = (-1, -1 / QQ)
fn hamilton() -> Arc<QuaternionAlgebraStructure<RationalCanonicalStructure>> {
    QuaternionAlgebraStructure::new(Rational::structure(), -Rational::ONE, -Rational::ONE)
}

fn quaternion(
    h: &Arc<QuaternionAlgebraStructure<RationalCanonicalStructure>>,
    x: i32,
    y: i32,
    z: i32,
    w: i32,
) -> QuaternionAlgebraElement<Rational> {
    h.from_components(
        Rational::from(x),
        Rational::from(y),
        Rational::from(z),
        Rational::from(w),
    )
}

#[test]
fn test_parse_quaternion_implicit_multiplication() {
    let h = hamilton();
    let q = parse_quaternion("2+3i+5j-2k", &h).unwrap();
    assert!(h.equal(&q, &quaternion(&h, 2, 3, 5, -2)));
}

#[test]
fn test_parse_quaternion_explicit_multiplication() {
    let h = hamilton();
    let q = parse_quaternion("2 + 3*i + 5*j - 2*k", &h).unwrap();
    assert!(h.equal(&q, &quaternion(&h, 2, 3, 5, -2)));
}

#[test]
fn test_parse_quaternion_basis_elements() {
    let h = hamilton();
    assert!(h.equal(&parse_quaternion("1", &h).unwrap(), &h.one()));
    assert!(h.equal(&parse_quaternion("i", &h).unwrap(), &h.i()));
    assert!(h.equal(&parse_quaternion("j", &h).unwrap(), &h.j()));
    assert!(h.equal(&parse_quaternion("k", &h).unwrap(), &h.k()));
    assert!(h.equal(&parse_quaternion("-i", &h).unwrap(), &h.neg(&h.i())));
}

#[test]
fn test_parse_quaternion_rational_coefficients() {
    let h = hamilton();
    let q = parse_quaternion("1/3 - 1/5*i + 1/7j - k", &h).unwrap();
    let expected = h.from_components(
        Rational::from_integers(1, 3),
        -Rational::from_integers(1, 5),
        Rational::from_integers(1, 7),
        -Rational::ONE,
    );
    assert!(h.equal(&q, &expected));
}

#[test]
fn test_parse_quaternion_collects_repeated_terms() {
    let h = hamilton();
    let q = parse_quaternion("i + 2 + 3i - 2", &h).unwrap();
    assert!(h.equal(&q, &quaternion(&h, 0, 4, 0, 0)));
}

#[test]
fn test_parse_quaternion_expands_brackets() {
    let h = hamilton();
    let q = parse_quaternion("2*(1 + i) - 3(j - k)", &h).unwrap();
    assert!(h.equal(&q, &quaternion(&h, 2, 2, -3, 3)));
}

#[test]
fn test_parse_quaternion_over_finite_field() {
    let f5 = ConwayFiniteFieldStructure::new(5, 1).unwrap();
    let algebra = QuaternionAlgebraStructure::new(f5.clone(), f5.from_int(2), f5.from_int(3));
    // 7 = 2 and 1/2 = 3 modulo 5
    let q = parse_quaternion("7 + 1/2*i", &algebra).unwrap();
    let expected = algebra.from_components(f5.from_int(2), f5.from_int(3), f5.zero(), f5.zero());
    assert!(algebra.equal(&q, &expected));
}

#[test]
fn test_parse_quaternion_nonlinear() {
    let h = hamilton();
    // i*j is not determined by the commutative polynomial i*j, so it is rejected
    assert!(parse_quaternion("i*j", &h).is_err());
    assert!(parse_quaternion("i^2", &h).is_err());
    assert!(parse_quaternion("(1 + i)*(1 + j)", &h).is_err());
}

#[test]
fn test_parse_quaternion_invalid() {
    let h = hamilton();
    // only i, j and k are valid variables
    assert!(parse_quaternion("1 + x", &h).is_err());
    // implicit multiplication of two variables is still not allowed
    assert!(parse_quaternion("ij", &h).is_err());
    assert!(parse_quaternion("1 +", &h).is_err());
}
