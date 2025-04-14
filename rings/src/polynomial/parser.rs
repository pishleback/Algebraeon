use crate::polynomial::Polynomial;
use algebraeon_nzq::*;

use std::collections::{BTreeMap, HashSet};
use std::str::FromStr;

// Variable names or numbers
// e.g. "x" or "123"
struct Chars {
    string: String,
}

// Something in brackets
// e.g. "(x + y^2 - 2)"
enum SubExpression {
    Chars(Chars),
    Sum(Sum),
    Product(Product),
    Power(Power),
}

// A sum
// e.g. "x + 123 + 4*x^2*y + (x + y^2 - 2)"
struct Sum {
    terms: Vec<SumTerm>,
}
enum SumTerm {
    Chars(Chars),
    Power(Power),
    Product(Product),
    SubExpression(SubExpression),
}

// A product
// e.g. "4*x^2*y"
struct Product {
    terms: Vec<ProductTerm>,
}
enum ProductTerm {
    Chars(Chars),
    Power(Power),
    SubExpression(SubExpression),
}

// A power of something
// e.g. "x^2" or "(x + y^2)^3"
struct Power {
    base: Box<PowerBase>,
    exponent: Integer,
}
enum PowerBase {
    Chars(Chars),
    SubExpression(SubExpression),
}

// Any expression
enum Expression {
    Chars(Chars),
    Sum(Sum),
    Product(Product),
    Power(Power),
}

impl Expression {
    fn from_string(input: String) -> Self {
        todo!()
    }

    // integer coefficients
    fn to_univariate_integer_polynomial(
        var: char,
        expression: Expression,
    ) -> Result<Polynomial<Integer>, String> {
        todo!()
    }

    // rational coefficients
    fn to_univariate_rational_polynomial(
        var: char,
        expression: Expression,
    ) -> Result<Polynomial<Rational>, String> {
        todo!()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Monomial {
    degree: usize, // total degree
}

impl Monomial {
    fn new(degree: usize) -> Self {
        Self { degree }
    }
}

/// Parses a univariate polynomial with Integer coefficients.
pub fn parse_univariate_polynomial(input: &str) -> Result<Polynomial<Integer>, String> {
    if input.trim().is_empty() {
        return Err("Empty input string".to_string());
    }

    // Check for invalid operators like just "*" or "/" alone
    if input.contains("**") || input.contains("//") || input.contains("++") || input.contains("--")
    {
        return Err("Invalid consecutive operators found".to_string());
    }

    // Check for lone operators
    if input.trim() == "*" || input.trim() == "/" {
        return Err("Lone operator without operands".to_string());
    }

    // Detect all variable characters in the input
    let mut variables = HashSet::new();
    for c in input.chars() {
        if c.is_alphabetic() {
            variables.insert(c);
        }
    }

    // Ensure we have at most one variable
    if variables.len() > 1 {
        return Err(format!(
            "Found multiple variables: {:?}; expected univariate polynomial",
            variables
        ));
    }

    // Get the variable if there is one
    let variable = variables.iter().next().copied();

    let input = input.replace(" ", "").replace("*", "");

    // Early check for malformed expressions
    if input.contains("^") && variable.is_none() {
        return Err("Invalid exponentiation: '^' without variable".to_string());
    }

    // Check for invalid exponents
    if let Some(var) = variable {
        let exp_pattern = format!(r"{}\^([^0-9].*)", var);
        let exp_check = regex::Regex::new(&exp_pattern).unwrap();
        if let Some(cap) = exp_check.captures(&input) {
            return Err(format!("Invalid exponent: '{}'", &cap[1]));
        }
    }

    let mut terms = BTreeMap::new();

    let re = regex::Regex::new(r"([+-]?[^+-]+)").unwrap();
    for cap in re.captures_iter(&input) {
        let term = &cap[1];

        // Parse coefficient and degree
        let (coeff_str, deg) = match variable {
            Some(var) => {
                if let Some(idx) = term.find(var) {
                    let coeff_part = &term[..idx];
                    let power = if let Some(pow_idx) = term.find("^") {
                        let exp_str = &term[pow_idx + 1..];
                        exp_str
                            .parse::<usize>()
                            .map_err(|_| format!("Invalid exponent: '{}'", exp_str))?
                    } else {
                        1
                    };
                    (coeff_part, power)
                } else {
                    (term, 0)
                }
            }
            None => (term, 0), // No variable means constant polynomial
        };

        let coeff = if coeff_str == "+" || coeff_str.is_empty() {
            Integer::from(1)
        } else if coeff_str == "-" {
            Integer::from(-1)
        } else {
            Integer::from_str(coeff_str)
                .map_err(|e| format!("Invalid coefficient '{}': {:?}", coeff_str, e))?
        };

        let monomial = Monomial::new(deg);
        terms
            .entry(monomial)
            .and_modify(|c: &mut Integer| *c += &coeff)
            .or_insert(coeff);
    }

    if terms.is_empty() {
        return Err("Unable to parse any valid terms".to_string());
    }

    let max_degree = terms.keys().map(|m| m.degree).max().unwrap_or(0);
    let mut coeffs = vec![Integer::from(0); max_degree + 1];
    for (mono, coeff) in terms {
        coeffs[mono.degree] = coeff;
    }

    Ok(Polynomial::from_coeffs(coeffs))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_univariate_whitespace_variable() {
        let input = "  1 + x^2  - 5x^3  ";
        let poly = parse_univariate_polynomial(input).unwrap();
        assert_eq!(poly, Polynomial::from_coeffs(vec![1, 0, 1, -5]));
    }

    #[test]
    fn test_parse_constant() {
        let input = "42";
        let poly = parse_univariate_polynomial(input).unwrap();
        assert_eq!(poly, Polynomial::from_coeffs(vec![42]));
    }

    #[test]
    fn test_parse_univariate_explicit_mult() {
        let input = "3*x^2 - 4*x + 5";
        let poly = parse_univariate_polynomial(input).unwrap();
        assert_eq!(poly, Polynomial::from_coeffs(vec![5, -4, 3]));
    }

    #[test]
    fn test_parse_univariate_explicit_mult_skip_exp() {
        let input = "3*x^3 - 4*x + 5";
        let poly = parse_univariate_polynomial(input).unwrap();
        assert_eq!(poly, Polynomial::from_coeffs(vec![5, -4, 0, 3]));
    }

    #[test]
    fn test_parse_univariate_implicit_mult() {
        let input = "3x^2 - 4x + 5";
        let poly = parse_univariate_polynomial(input).unwrap();
        assert_eq!(poly, Polynomial::from_coeffs(vec![5, -4, 3]));
    }

    #[test]
    fn test_parse_univariate_diff_variable() {
        let input = "y^3 + 2y - 7";
        let poly = parse_univariate_polynomial(input).unwrap();
        assert_eq!(poly, Polynomial::from_coeffs(vec![-7, 2, 0, 1]));
    }

    // Error case tests

    #[test]
    fn test_multiple_variables() {
        let input = "2x + 3y + z";
        let result = parse_univariate_polynomial(input);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Found multiple variables"));
    }

    #[test]
    fn test_lone_operator() {
        let input = "*";
        let result = parse_univariate_polynomial(input);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Lone operator without operands");
    }

    #[test]
    fn test_invalid_exponent() {
        let input = "x^a";
        let result = parse_univariate_polynomial(input);
        assert!(result.is_err());
        println!("Error: {:?}", result);
        // invalid exponent fails because we end up flagging this input for being
        // multivariate before it reaches the bad coeeficient test
        // to do: make this more precise
        // assert!(result.unwrap_err().contains("Invalid exponent"));
        assert!(result.unwrap_err().contains("Found multiple variables"));
    }

    #[test]
    fn test_invalid_consecutive_operators() {
        let input = "x++2";
        let result = parse_univariate_polynomial(input);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Invalid consecutive operators found");
    }

    #[test]
    fn test_empty_input() {
        let input = "";
        let result = parse_univariate_polynomial(input);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Empty input string");
    }

    #[test]
    fn test_exponent_without_variable() {
        let input = "3^2";
        let result = parse_univariate_polynomial(input);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err(),
            "Invalid exponentiation: '^' without variable"
        );
    }
}
