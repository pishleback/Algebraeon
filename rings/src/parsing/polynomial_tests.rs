use crate::parsing::ast::*;
use crate::parsing::polynomial::*;
use crate::polynomial::{MultiPolynomial, Polynomial, Variable};
use crate::structure::{MetaRing, MetaSemiRing};
use algebraeon_nzq::*;
use std::collections::HashMap;
use std::str::FromStr;

#[test]
fn test_integer_polynomial_basic() {
    let result = parse_integer_polynomial("3*x^2 + 2*x - 1", "x").unwrap();
    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![Integer::from(-1), Integer::from(2), Integer::from(3)])
    );
}

#[test]
fn test_integer_polynomial_different_var() {
    let result = parse_integer_polynomial("y^3 + 2*y - 7", "y").unwrap();
    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![
            Integer::from(-7),
            Integer::from(2),
            Integer::from(0),
            Integer::from(1)
        ])
    );
}

#[test]
fn test_integer_polynomial_constant() {
    let result = parse_integer_polynomial("42", "x").unwrap();
    assert_eq!(result, Polynomial::from_coeffs(vec![Integer::from(42)]));
}

#[test]
fn test_rational_polynomial_basic() {
    let result = parse_rational_polynomial("3*x^2 + 2*x - 1", "x").unwrap();
    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![
            Rational::from(-1),
            Rational::from(2),
            Rational::from(3)
        ])
    );
}

#[test]
fn test_rational_polynomial_fractions() {
    let result = parse_rational_polynomial("2/3*x^3 + 1/4*x - 7/6", "x").unwrap();

    // Create expected rational coefficients
    let neg_seven_sixths = Rational::from_integers(Integer::from(-7), Integer::from(6));
    let one_fourth = Rational::from_integers(Integer::from(1), Integer::from(4));
    let zero = Rational::from(0);
    let two_thirds = Rational::from_integers(Integer::from(2), Integer::from(3));

    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![neg_seven_sixths, one_fourth, zero, two_thirds])
    );
}

#[test]
fn test_rational_polynomial_expansion() {
    let result = parse_rational_polynomial("(x + 2)^3", "x").unwrap();

    // (x + 2)^3 = x^3 + 6x^2 + 12x + 8
    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![
            Rational::from(8),
            Rational::from(12),
            Rational::from(6),
            Rational::from(1)
        ])
    );
}

#[test]
fn test_rational_polynomial_nested_subexpression() {
    let result = parse_rational_polynomial("(((x + 2)))", "x").unwrap();

    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![Rational::from(2), Rational::from(1)])
    );
}

#[test]
fn test_integer_polynomial_expansion() {
    let result = parse_integer_polynomial("(x + 2)^3", "x").unwrap();

    // (x + 2)^3 = x^3 + 6x^2 + 12x + 8
    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![
            Integer::from(8),
            Integer::from(12),
            Integer::from(6),
            Integer::from(1)
        ])
    );
}

// New test cases for multivariate polynomials
#[test]
fn test_multivariate_integer_polynomial_basic() {
    // Create variables first
    let x_var = Variable::new("x");
    let y_var = Variable::new("y");
    let variable_mapping = [("x", x_var.clone()), ("y", y_var.clone())].into();

    let result = parse_multivariate_integer_polynomial("x + y", variable_mapping).unwrap();

    // Create expected polynomial using the same variables
    let x = MultiPolynomial::<Integer>::var(x_var);
    let y = MultiPolynomial::<Integer>::var(y_var);
    let expected = MultiPolynomial::add(&x, &y);

    assert_eq!(result, expected);
}

#[test]
fn test_multivariate_integer_polynomial_product() {
    // Create variables first
    let x_var = Variable::new("x");
    let y_var = Variable::new("y");
    let variable_mapping = [("x", x_var.clone()), ("y", y_var.clone())].into();

    let result = parse_multivariate_integer_polynomial("x * y", variable_mapping).unwrap();

    // Create expected polynomial using the same variables
    let x = MultiPolynomial::<Integer>::var(x_var);
    let y = MultiPolynomial::<Integer>::var(y_var);
    let expected = MultiPolynomial::mul(&x, &y);

    assert_eq!(result, expected);
}

#[test]
fn test_multivariate_integer_polynomial_power() {
    // Create variables first
    let x_var = Variable::new("x");
    let y_var = Variable::new("y");
    let variable_mapping = [("x", x_var.clone()), ("y", y_var.clone())].into();

    let result = parse_multivariate_integer_polynomial("x^2 + y^3", variable_mapping).unwrap();

    // Create expected polynomial using the same variables
    let x = MultiPolynomial::<Integer>::var(x_var);
    let y = MultiPolynomial::<Integer>::var(y_var);
    let x_squared = x.nat_pow(&Natural::from(2u32));
    let y_cubed = y.nat_pow(&Natural::from(3u32));
    let expected = MultiPolynomial::add(&x_squared, &y_cubed);

    assert_eq!(result, expected);
}

#[test]
fn test_multivariate_integer_polynomial_complex() {
    // Create variables first
    let x_var = Variable::new("x");
    let y_var = Variable::new("y");
    let variable_mapping = [("x", x_var.clone()), ("y", y_var.clone())].into();

    let result =
        parse_multivariate_integer_polynomial("3*x^2*y + 2*x*y^2 - x + 5", variable_mapping)
            .unwrap();

    // Create expected polynomial using the same variables
    let x = MultiPolynomial::<Integer>::var(x_var);
    let y = MultiPolynomial::<Integer>::var(y_var);
    let constant_3 = MultiPolynomial::<Integer>::constant(Integer::from(3));
    let constant_2 = MultiPolynomial::<Integer>::constant(Integer::from(2));
    let constant_5 = MultiPolynomial::<Integer>::constant(Integer::from(5));

    let x_squared = x.nat_pow(&Natural::from(2u32));
    let y_squared = y.nat_pow(&Natural::from(2u32));

    let term1 = MultiPolynomial::mul(&MultiPolynomial::mul(&constant_3, &x_squared), &y);
    let term2 = MultiPolynomial::mul(&MultiPolynomial::mul(&constant_2, &x), &y_squared);
    let neg_x = MultiPolynomial::neg(&x);

    let sum1 = MultiPolynomial::add(&term1, &term2);
    let sum2 = MultiPolynomial::add(&sum1, &neg_x);
    let expected = MultiPolynomial::add(&sum2, &constant_5);

    assert_eq!(result, expected);
}

#[test]
fn test_multivariate_rational_polynomial_basic() {
    // Create variables first
    let x_var = Variable::new("x");
    let y_var = Variable::new("y");
    let variable_mapping = [("x", x_var.clone()), ("y", y_var.clone())].into();

    let result = parse_multivariate_rational_polynomial("x + y", variable_mapping).unwrap();

    // Create expected polynomial using the same variables
    let x = MultiPolynomial::<Rational>::var(x_var);
    let y = MultiPolynomial::<Rational>::var(y_var);
    let expected = MultiPolynomial::add(&x, &y);

    assert_eq!(result, expected);
}

#[test]
fn test_multivariate_rational_polynomial_fractions() {
    // Create variables first
    let x_var = Variable::new("x");
    let y_var = Variable::new("y");
    let variable_mapping = [("x", x_var.clone()), ("y", y_var.clone())].into();

    let result = parse_multivariate_rational_polynomial("1/2*x + 3/4*y", variable_mapping).unwrap();

    // Create expected polynomial using the same variables
    let x = MultiPolynomial::<Rational>::var(x_var);
    let y = MultiPolynomial::<Rational>::var(y_var);
    let half = MultiPolynomial::<Rational>::constant(Rational::from_integers(
        Integer::from(1),
        Integer::from(2),
    ));
    let three_fourths = MultiPolynomial::<Rational>::constant(Rational::from_integers(
        Integer::from(3),
        Integer::from(4),
    ));

    let term1 = MultiPolynomial::mul(&half, &x);
    let term2 = MultiPolynomial::mul(&three_fourths, &y);
    let expected = MultiPolynomial::add(&term1, &term2);

    assert_eq!(result, expected);
}

#[test]
fn test_multivariate_polynomial_expansion() {
    // Create variables first
    let x_var = Variable::new("x");
    let y_var = Variable::new("y");
    let variable_mapping = [("x", x_var.clone()), ("y", y_var.clone())].into();

    let result = parse_multivariate_integer_polynomial("(x + y)^2", variable_mapping).unwrap();

    // Create expected polynomial: x^2 + 2*x*y + y^2
    let x = MultiPolynomial::<Integer>::var(x_var);
    let y = MultiPolynomial::<Integer>::var(y_var);
    let sum_xy = MultiPolynomial::add(&x, &y);
    let expected = sum_xy.nat_pow(&Natural::from(2u32));

    assert_eq!(result, expected);
}

#[test]
fn test_multivariate_polynomial_large_expression() {
    // Create variables first
    let x_var = Variable::new("x");
    let y_var = Variable::new("y");
    let z_var = Variable::new("z");
    let variable_mapping = [
        ("x", x_var.clone()),
        ("y", y_var.clone()),
        ("z", z_var.clone()),
    ]
    .into();

    let result = parse_multivariate_integer_polynomial("(x + y + z)^30", variable_mapping).unwrap();

    // Create expected polynomial using the same variables
    let x = MultiPolynomial::<Integer>::var(x_var);
    let y = MultiPolynomial::<Integer>::var(y_var);
    let z = MultiPolynomial::<Integer>::var(z_var);
    let sum_xy = MultiPolynomial::add(&x, &y);
    let sum_xyz = MultiPolynomial::add(&sum_xy, &z);
    let expected = sum_xyz.nat_pow(&Natural::from(30u32));

    assert_eq!(result, expected);
}

// New test cases for multi-character variables with braces
#[test]
fn test_multi_character_variables() {
    let result = parse_rational_polynomial("{var}^2 + 3*{var} + 1", "var").unwrap();

    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![
            Rational::from(1),
            Rational::from(3),
            Rational::from(1)
        ])
    );
}

#[test]
fn test_multivariate_multi_character_variables() {
    // Create variables first
    let var1 = Variable::new("var1");
    let var2 = Variable::new("var2");
    let variable_mapping = [("var1", var1.clone()), ("var2", var2.clone())].into();

    let result =
        parse_multivariate_integer_polynomial("{var1} + {var2}", variable_mapping).unwrap();

    // Create expected polynomial using the same variables
    let x = MultiPolynomial::<Integer>::var(var1);
    let y = MultiPolynomial::<Integer>::var(var2);
    let expected = MultiPolynomial::add(&x, &y);

    assert_eq!(result, expected);
}

#[test]
fn test_mixed_variable_formats() {
    // This should work because we're still using "x" as the variable name
    // even though it's represented as {x} in the expression
    let result = parse_rational_polynomial("{x}^2 + 3*{x} + 1", "x").unwrap();

    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![
            Rational::from(1),
            Rational::from(3),
            Rational::from(1)
        ])
    );
}

#[test]
fn test_mixed_variable_formats_multivariate() {
    // Create variables first
    let x_var = Variable::new("x");
    let y_var = Variable::new("y");
    let variable_mapping = [("x", x_var.clone()), ("y", y_var.clone())].into();

    let result = parse_multivariate_integer_polynomial("{x} + y", variable_mapping).unwrap();

    // Create expected polynomial using the same variables
    let x = MultiPolynomial::<Integer>::var(x_var);
    let y = MultiPolynomial::<Integer>::var(y_var);
    let expected = MultiPolynomial::add(&x, &y);

    assert_eq!(result, expected);
}

#[test]
fn test_multiple_multi_character_variables() {
    let result = parse_rational_polynomial("{foo} + {foo}*{bar}", "foo");

    // This should fail because it contains multiple variables
    assert!(result.is_err());
}

#[test]
fn test_invalid_polynomial() {
    let result = parse_rational_polynomial("x^-2", "x");
    assert!(result.is_err());

    let result = parse_rational_polynomial("x^(1/2)", "x");
    assert!(result.is_err());

    let result = parse_rational_polynomial("1/(x+1)", "x");
    assert!(result.is_err());

    // Test multivariate versions too
    let x_var = Variable::new("x");
    let variable_mapping = [("x", x_var)].into();

    let result = parse_multivariate_rational_polynomial("x^-2", variable_mapping);
    assert!(result.is_err());

    let x_var = Variable::new("x");
    let variable_mapping = [("x", x_var)].into();

    let result = parse_multivariate_rational_polynomial("x^(1/2)", variable_mapping);
    assert!(result.is_err());
}

#[test]
fn test_invalid_implicit_multiplication() {
    // These should fail because implicit multiplication is not allowed
    let result = parse_rational_polynomial("2x", "x");
    assert!(result.is_err());

    let result = parse_rational_polynomial("xy", "x");
    assert!(result.is_err());

    let result = parse_rational_polynomial("x{foo}", "x");
    assert!(result.is_err());

    // Test multivariate versions too - these should fail at parsing level
    // so we can use empty variable mappings
    let empty_mapping = HashMap::new();

    let result = parse_multivariate_rational_polynomial("2x", empty_mapping);
    assert!(result.is_err());

    let empty_mapping = HashMap::new();
    let result = parse_multivariate_rational_polynomial("x{foo}", empty_mapping);
    assert!(result.is_err());
}

#[test]
fn test_integer_polynomial_big_constant() {
    let result = parse_integer_polynomial(
        "123456789123456789123456789123456789123456789123456789123456789",
        "x",
    )
    .unwrap();
    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![
            Integer::from_str("123456789123456789123456789123456789123456789123456789123456789")
                .unwrap()
        ])
    );
}

#[test]
fn test_rational_polynomial_big_constant() {
    let result = parse_rational_polynomial(
        "123456789123456789123456789123456789123456789123456789123456789/123456789123456789123456789123456789123456789123456789123456789",
        "x",
    )
        .unwrap();
    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![
            Rational::from_str(
                "123456789123456789123456789123456789123456789123456789123456789/123456789123456789123456789123456789123456789123456789123456789"
            )
                .unwrap()
        ])
    );
}

#[test]
fn test_rational_polynomial_big_coefficients() {
    let result = parse_rational_polynomial("123456789123456789123456789123456789123456789123456789123456789*x^2 + 123456789123456789123456789123456789123456789123456789123456789*x - 123456789123456789123456789123456789123456789123456789123456789", "x").unwrap();
    assert_eq!(
        result,
        Polynomial::from_coeffs(vec![
            Rational::from_str("-123456789123456789123456789123456789123456789123456789123456789")
                .unwrap(),
            Rational::from_str("123456789123456789123456789123456789123456789123456789123456789")
                .unwrap(),
            Rational::from_str("123456789123456789123456789123456789123456789123456789123456789")
                .unwrap()
        ])
    );
}
