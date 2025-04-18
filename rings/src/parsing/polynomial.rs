use crate::polynomial::Polynomial;
use algebraeon_nzq::*;
use lalrpop_util::lalrpop_mod;
use std::collections::{BTreeMap, HashSet};
use std::str::FromStr;
use std::fmt;

lalrpop_mod!(polynomial_parser, "/parsing/polynomial_grammar.rs"); // synthesized by LALRPOP


#[derive(Debug, Clone, PartialEq)]
pub struct Variable {
    pub name: String,
}

impl fmt::Display for Variable {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.name)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Number {
    pub value: i64,
}

impl fmt::Display for Number {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    Var(Variable),
    Num(Number),
    Sum(Sum),
    Product(Product),
    Power(Power),
    Grouped(Box<Expr>),
    Neg(Box<Expr>), // Added for unary negation
}

impl Expr {
    pub fn validate_polynomial(&self) -> Result<(), String> {
        match self {
            Expr::Power(p) => p.validate_polynomial(),
            Expr::Product(p) => p.validate_polynomial(),
            Expr::Sum(s) => s.validate_polynomial(),
            Expr::Neg(e) => e.validate_polynomial(),
            Expr::Grouped(e) => e.validate_polynomial(),
            _ => Ok(()),
        }
    }
}

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::Var(v) => write!(f, "{}", v),
            Self::Num(n) => write!(f, "{}", n),
            Self::Sum(s) => write!(f, "{}", s),
            Self::Product(p) => write!(f, "{}", p),
            Self::Power(p) => write!(f, "{}", p),
            Self::Grouped(e) => write!(f, "({})", e),
            Self::Neg(e) => write!(f, "-{}", e),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Sum {
    pub left: Box<Expr>,
    pub right: Box<Expr>,
}

impl Sum {
    fn validate_polynomial(&self) -> Result<(), String> {
        self.left.validate_polynomial()?;
        self.right.validate_polynomial()
    }
}

impl fmt::Display for Sum {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} + {}", self.left, self.right)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Product {
    pub left: Box<Expr>,
    pub right: Box<Expr>,
}

impl Product {
    fn validate_polynomial(&self) -> Result<(), String> {
        if let Expr::Power(p) = self.right.as_ref() {
            if let Expr::Num(n) = p.exponent.as_ref() {
                if n.value == -1 {
                    return Err("Division not allowed in polynomials".to_string());
                }
            }
        }
        self.left.validate_polynomial()?;
        self.right.validate_polynomial()
    }
}

impl fmt::Display for Product {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} * {}", self.left, self.right)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Power {
    pub base: Box<Expr>,
    pub exponent: Box<Expr>,
}

impl Power {
    fn validate_polynomial(&self) -> Result<(), String> {
        if let Expr::Num(n) = self.exponent.as_ref() {
            if n.value < 0 {
                return Err("Negative exponents not allowed in polynomials".to_string());
            }
        }
        self.base.validate_polynomial()
    }
}

impl fmt::Display for Power {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}^{}", self.base, self.exponent)
    }
}

impl Expr {
    fn from_string(input: String) -> Self {
        todo!()
    }

    // integer coefficients
    fn to_univariate_integer_polynomial(
        var: char,
        expression: Expr,
    ) -> Result<Polynomial<Integer>, String> {
        todo!()
    }

    // rational coefficients
    fn to_univariate_rational_polynomial(
        var: char,
        expression: Expr,
    ) -> Result<Polynomial<Rational>, String> {
        todo!()
    }
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
        // multivariate before it reaches the bad coeficient test
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
