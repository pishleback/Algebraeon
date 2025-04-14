use crate::polynomial::Polynomial;
use algebraeon_nzq::*;

use std::collections::BTreeMap;
use std::str::FromStr;


#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Monomial {
    degree: usize, // total degree
}

impl Monomial {
    fn new(degree: usize) -> Self {
        Self { degree }
    }
}

/// Parses a univariate polynomial in x with Integer coefficients.
pub fn parse_univariate_polynomial(input: &str) -> Result<Polynomial<Integer>, String> {
    let input = input.replace(" ", "").replace("*", "");

    let mut terms = BTreeMap::new();

    let re = regex::Regex::new(r"([+-]?[^+-]+)").unwrap();
    for cap in re.captures_iter(&input) {
        let term = &cap[1];

        let (coeff_str, deg) = if let Some(idx) = term.find('x') {
            let coeff_part = &term[..idx];
            let power = if let Some(pow_idx) = term.find("^") {
                term[pow_idx + 1..].parse::<usize>().map_err(|e| e.to_string())?
            } else {
                1
            };
            (coeff_part, power)
        } else {
            (term, 0)
        };

        let coeff = if coeff_str == "+" || coeff_str.is_empty() {
            Integer::from(1)
        } else if coeff_str == "-" {
            Integer::from(-1)
        } else {
            Integer::from_str(coeff_str).map_err(|e| format!("Bad coeff '{}': {:?}", coeff_str, e))?
        };

        let monomial = Monomial::new(deg);
        terms
            .entry(monomial)
            .and_modify(|c: &mut Integer| *c += &coeff)
            .or_insert(coeff);
    }

    let max_degree = terms.keys().map(|m| m.degree).max().unwrap_or(0);
    let mut coeffs = vec![Integer::from(0); max_degree + 1];
    for (mono, coeff) in terms {
        coeffs[mono.degree] = coeff;
    }

    Ok(Polynomial::from_coeffs(coeffs))
}


#[test]
fn test_parse_univariate_explicit_mult() {
    let input = "3*x^2 - 4*x + 5";
    let poly = parse_univariate_polynomial(input).unwrap();
    assert_eq!(poly, Polynomial::from_coeffs(vec![5, -4, 3]));
}

#[test]
fn test_parse_univariate_implicit_mult() {
    let input = "3x^2 - 4x + 5";
    let poly = parse_univariate_polynomial(input).unwrap();
    assert_eq!(poly, Polynomial::from_coeffs(vec![5, -4, 3]));
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     #[test]
//     fn test_parser() {
//         // Simple constant polynomial
//         {
//             let poly = parse_univariate_polynomial("42").unwrap();
//             let expected = Polynomial::from_coeffs(vec![Integer::from(42)]);
//             assert_eq!(poly, expected);
//         }
//         
//         // Test with whitespace variations
//         {
//             let poly = parse_univariate_polynomial("  1 + x^2  - 5x^3  ").unwrap();
//             let expected = Polynomial::from_coeffs(vec![
//                 Integer::from(1),
//                 Integer::from(0),
//                 Integer::from(1),
//                 Integer::from(-5),
//             ]);
//             assert_eq!(poly, expected);
//         }
//     }
// }
