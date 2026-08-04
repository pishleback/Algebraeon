use crate::parsing::polynomial::parse_multivariate_rational_polynomial;
use crate::polynomial::Variable;
use crate::quaternion_algebra::{QuaternionAlgebraElement, QuaternionAlgebraStructure};
use crate::structure::{FieldSignature, SemiModuleSignature};
use algebraeon_structures::*;
use std::collections::HashMap;
use std::sync::Arc;

/// Create an element of a quaternion algebra from a string such as `"2 + 3i + 5j - 2k"`.
///
/// The expression must be linear in `i`, `j` and `k`, since the quaternion algebra is not
/// commutative and so, for example, `i*j` is not determined by the parsed expression alone.
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
/// ```
pub fn parse_quaternion<Field: FieldSignature>(
    quaternion_str: &str,
    quaternion_algebra: &Arc<QuaternionAlgebraStructure<Field>>,
) -> Result<QuaternionAlgebraElement<Field::Elem>, String> {
    let i_var = Variable::new("i");
    let j_var = Variable::new("j");
    let k_var = Variable::new("k");
    let variable_mapping = HashMap::from([
        ("i", i_var.clone()),
        ("j", j_var.clone()),
        ("k", k_var.clone()),
    ]);

    let polynomial = parse_multivariate_rational_polynomial(quaternion_str, variable_mapping)?;

    // The coefficients of 1, i, j and k.
    let mut components = [
        Rational::ZERO,
        Rational::ZERO,
        Rational::ZERO,
        Rational::ZERO,
    ];
    for term in &polynomial.terms {
        let monomial = &term.monomial;
        let component = match monomial.degree() {
            0 => 0,
            1 => {
                if monomial.get_var_pow(&i_var) == 1 {
                    1
                } else if monomial.get_var_pow(&j_var) == 1 {
                    2
                } else {
                    debug_assert_eq!(monomial.get_var_pow(&k_var), 1);
                    3
                }
            }
            _ => {
                return Err(format!(
                    "Quaternion expressions must be linear in i, j and k, but found the term '{}'",
                    monomial
                ));
            }
        };
        components[component] += term.coeff.clone();
    }

    let base_field = quaternion_algebra.ring();
    let [x, y, z, w] = components.map(|c| {
        base_field
            .try_from_rat(&c)
            .ok_or_else(|| format!("No element of the base field is equal to {}", c))
    });

    Ok(quaternion_algebra.from_components(x?, y?, z?, w?))
}
