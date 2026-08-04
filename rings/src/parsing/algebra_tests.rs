use crate::finite_fields::conway_finite_fields::ConwayFiniteFieldStructure;
use crate::parsing::parse_algebra_expression;
use crate::polynomial::{Polynomial, ToPolynomialSignature};
use algebraeon_structures::*;
use std::collections::HashMap;

#[test]
fn test_parse_algebra_expression_in_a_polynomial_ring() {
    let polynomials = Rational::structure().polynomials();
    let variables = HashMap::from([("x", polynomials.var())]);

    let f = parse_algebra_expression("(x + 2)*(x - 2) + 1/2", &polynomials, &variables).unwrap();
    let expected = Polynomial::from_coeffs(vec![
        Rational::from_integers(-7, 2),
        Rational::ZERO,
        Rational::ONE,
    ]);
    assert!(polynomials.equal(&f, &expected));
}

#[test]
fn test_parse_algebra_expression_unknown_variable() {
    let polynomials = Rational::structure().polynomials();
    let variables = HashMap::from([("x", polynomials.var())]);

    assert!(parse_algebra_expression("x + y", &polynomials, &variables).is_err());
}

#[test]
fn test_parse_algebra_expression_coefficient_not_in_base_field() {
    let f5 = ConwayFiniteFieldStructure::new(5, 1).unwrap().polynomials();
    let variables = HashMap::from([("x", f5.var())]);

    // 1/5 does not exist in characteristic 5
    assert!(parse_algebra_expression("1/5*x", &f5, &variables).is_err());
}
