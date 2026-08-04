use crate::parsing::ast::*;
use crate::parsing::polynomial::parse_expression;
use crate::structure::{FieldSignature, RingSignature, SemiModuleSignature};
use algebraeon_structures::*;
use std::collections::HashMap;
use std::sync::Arc;

/// Evaluate a string such as `"2 + 3i + 5j - 2k"` in an algebra over a field, given the value of
/// each variable appearing in it.
///
/// Unlike parsing into a polynomial, the order of the factors of each product is preserved, so
/// this works for algebras with non-commutative multiplication such as quaternion algebras.
///
/// ```
/// use algebraeon_rings::parsing::parse_algebra_expression;
/// use algebraeon_rings::quaternion_algebra::QuaternionAlgebraStructure;
/// use algebraeon_rings::structure::*;
/// use algebraeon_structures::*;
/// use std::collections::HashMap;
///
/// // Hamilton quaternion algebra: H = (-1, -1 / QQ)
/// let h = QuaternionAlgebraStructure::new(Rational::structure(), -Rational::ONE, -Rational::ONE);
/// let variables = HashMap::from([("u", h.i()), ("v", h.j())]);
///
/// // uv = -vu, so uv - vu = 2ij = 2k
/// let value = parse_algebra_expression("u*v - v*u", &h, &variables).unwrap();
/// assert!(h.equal(&value, &h.scalar_mul(&h.k(), &Rational::TWO)));
/// ```
pub fn parse_algebra_expression<
    Field: FieldSignature,
    Algebra: RingSignature + SemiModuleSignature<Field>,
>(
    expression_str: &str,
    algebra: &Arc<Algebra>,
    variables: &HashMap<&str, Algebra::Elem>,
) -> Result<Algebra::Elem, String> {
    let expression = parse_expression(expression_str)?;
    expression.validate()?;
    evaluate(&expression, algebra, variables)
}

fn evaluate<Field: FieldSignature, Algebra: RingSignature + SemiModuleSignature<Field>>(
    expression: &Expr,
    algebra: &Arc<Algebra>,
    variables: &HashMap<&str, Algebra::Elem>,
) -> Result<Algebra::Elem, String> {
    match expression {
        Expr::Var(v) => match variables.get(v.name.as_str()) {
            Some(value) => Ok(value.clone()),
            None => Err(format!("Variable '{}' has no value", v.name)),
        },
        Expr::Num(n) => {
            let scalar = scalar(n, algebra)?;
            Ok(algebra.scalar_mul(&algebra.one(), &scalar))
        }
        Expr::Sum(s) => {
            let mut result = algebra.zero();
            for term in &s.terms {
                let value = evaluate(&term.term, algebra, variables)?;
                result = if term.sign {
                    algebra.add(&result, &value)
                } else {
                    algebra.sub(&result, &value)
                };
            }
            Ok(result)
        }
        Expr::Product(p) => {
            // Multiplication by a numeric coefficient is scalar multiplication. As well as being
            // faster, this means that expressions which are linear in the variables can be
            // evaluated in algebras whose multiplication is only partially implemented.
            match (p.left.as_ref(), p.right.as_ref()) {
                (Expr::Num(n), other) | (other, Expr::Num(n)) => {
                    let scalar = scalar(n, algebra)?;
                    Ok(algebra.scalar_mul(&evaluate(other, algebra, variables)?, &scalar))
                }
                (left, right) => Ok(algebra.mul(
                    &evaluate(left, algebra, variables)?,
                    &evaluate(right, algebra, variables)?,
                )),
            }
        }
        Expr::Power(p) => {
            let base = evaluate(&p.base, algebra, variables)?;
            match p.exponent.as_ref() {
                // Validation has already checked that the exponent is a non-negative integer
                Expr::Num(n) => {
                    let exponent: usize = (&n.numerator)
                        .try_into()
                        .map_err(|_| format!("Exponent {} is too large", n.numerator))?;
                    Ok(algebra.nat_pow(&base, &Natural::from(exponent)))
                }
                _ => Err("Exponents must be integer constants".to_string()),
            }
        }
        Expr::Grouped(e) => evaluate(e, algebra, variables),
    }
}

/// The value of a number in the base field of the algebra.
fn scalar<Field: FieldSignature, Algebra: RingSignature + SemiModuleSignature<Field>>(
    number: &Number,
    algebra: &Arc<Algebra>,
) -> Result<Field::Elem, String> {
    let value = Rational::from_integers(number.numerator.clone(), number.denominator.clone());
    algebra
        .ring()
        .try_from_rat(&value)
        .ok_or_else(|| format!("No element of the base field is equal to {}", value))
}
