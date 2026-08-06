use crate::parsing::algebra::parse_algebra_expression;
use crate::quaternion_algebra::{QuaternionAlgebraElement, QuaternionAlgebraStructure};
use crate::structure::FieldSignature;
use std::collections::HashMap;
use std::sync::Arc;

/// Create an element of a quaternion algebra from a string such as `"2 + 3i + 5j - 2k"`.
///
/// The order of the factors of each product is respected, so `"i*j"` and `"j*i"` give different
/// elements.
///
/// ```
/// use algebraeon_rings::parsing::parse_quaternion;
/// use algebraeon_rings::quaternion_algebra::QuaternionAlgebraStructure;
/// use algebraeon_rings::structure::*;
/// use algebraeon_structures::*;
///
/// // Hamilton quaternion algebra: H = (-1, -1 / QQ)
/// let h = QuaternionAlgebraStructure::new(Rational::structure(), -Rational::ONE, -Rational::ONE);
/// let q = parse_quaternion("2+3i+5j-2k", &h).unwrap();
/// let expected = h.from_components(
///     Rational::from(2),
///     Rational::from(3),
///     Rational::from(5),
///     Rational::from(-2),
/// );
/// assert!(h.equal(&q, &expected));
///
/// // ij = k and ji = -k
/// let ij = parse_quaternion("i*j", &h).unwrap();
/// assert!(h.equal(&ij, &h.k()));
/// assert!(h.equal(&parse_quaternion("j*i", &h).unwrap(), &h.neg(&ij)));
/// ```
pub fn parse_quaternion<Field: FieldSignature>(
    quaternion_str: &str,
    quaternion_algebra: &Arc<QuaternionAlgebraStructure<Field>>,
) -> Result<QuaternionAlgebraElement<Field::Elem>, String> {
    let variables = HashMap::from([
        ("i", quaternion_algebra.i()),
        ("j", quaternion_algebra.j()),
        ("k", quaternion_algebra.k()),
    ]);
    parse_algebra_expression(quaternion_str, quaternion_algebra, &variables)
}
