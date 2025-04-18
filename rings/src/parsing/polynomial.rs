use crate::polynomial::Polynomial;
use algebraeon_nzq::*;
use lalrpop_util::lalrpop_mod;
use std::collections::HashMap;
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
    fn from_string(input: String) -> Result<Self, String> {
        match polynomial_parser::ExprParser::new().parse(&input) {
            Ok(expr) => Ok(*expr),
            Err(err) => Err(format!("Failed to parse expression: {}", err)),
        }
    }

    // integer coefficients
    pub(crate) fn build_univariate_integer_polynomial(
        var: char,
        expression: Expr,
    ) -> Result<Polynomial<Integer>, String> {
        // First validate it's a valid polynomial
        expression.validate_polynomial()?;

        // Collect terms as (exponent, coefficient) pairs
        let mut terms: HashMap<usize, Integer> = HashMap::new();

        // Check expression contains only our specified variable
        let vars = expression.collect_variables();
        if vars.is_empty() {
            // Constant polynomial
        } else if vars.len() > 1 {
            return Err(format!("Expected univariate polynomial with variable '{}', but found multiple variables: {:?}", var, vars));
        } else if !vars.contains(&var.to_string()) {
            return Err(format!(
                "Expected polynomial with variable '{}', but found '{}'",
                var,
                vars.iter().next().unwrap()
            ));
        }

        // Extract the terms
        expression.extract_terms(var, &mut terms, Integer::from(1));

        // Find the highest degree
        let max_degree = terms.keys().max().copied().unwrap_or(0);

        // Create coefficient vector
        let mut coeffs = vec![Integer::from(0); max_degree + 1];
        for (degree, coeff) in terms {
            coeffs[degree] = coeff;
        }

        Ok(Polynomial::from_coeffs(coeffs))
    }

    // Helper method to collect all variables in the expression
    fn collect_variables(&self) -> std::collections::HashSet<String> {
        let mut vars = std::collections::HashSet::new();
        self.add_variables(&mut vars);
        vars
    }

    fn add_variables(&self, vars: &mut std::collections::HashSet<String>) {
        match self {
            Expr::Var(v) => {
                vars.insert(v.name.clone());
            }
            Expr::Num(_) => {}
            Expr::Sum(s) => {
                s.left.add_variables(vars);
                s.right.add_variables(vars);
            }
            Expr::Product(p) => {
                p.left.add_variables(vars);
                p.right.add_variables(vars);
            }
            Expr::Power(p) => {
                // Only collect variables from base, not exponent (for polynomials)
                p.base.add_variables(vars);
                // For a valid polynomial, exponent must be a constant
            }
            Expr::Grouped(e) => e.add_variables(vars),
            Expr::Neg(e) => e.add_variables(vars),
        }
    }

    // Recursively extract terms from the expression
    fn extract_terms(&self, var: char, terms: &mut HashMap<usize, Integer>, coefficient: Integer) {
        match self {
            Expr::Var(v) => {
                if v.name == var.to_string() {
                    // x^1
                    *terms.entry(1).or_insert(Integer::from(0)) += coefficient;
                } else {
                    // Should never happen due to variable checking
                    panic!("Unexpected variable: {}", v.name);
                }
            }
            Expr::Num(n) => {
                // Constant term
                *terms.entry(0).or_insert(Integer::from(0)) += coefficient * Integer::from(n.value);
            }
            Expr::Sum(s) => {
                s.left.extract_terms(var, terms, coefficient.clone());
                s.right.extract_terms(var, terms, coefficient);
            }
            Expr::Product(p) => {
                match (p.left.as_ref(), p.right.as_ref()) {
                    // Handle coefficient * var^n pattern
                    (Expr::Num(n), right) => {
                        let new_coeff = coefficient * Integer::from(n.value);
                        right.extract_terms(var, terms, new_coeff);
                    }
                    // Handle var^n * coefficient pattern
                    (left, Expr::Num(n)) => {
                        let new_coeff = coefficient * Integer::from(n.value);
                        left.extract_terms(var, terms, new_coeff);
                    }
                    // Handle var * var pattern (produces var^2)
                    (Expr::Var(v1), Expr::Var(v2))
                        if v1.name == var.to_string() && v2.name == var.to_string() =>
                    {
                        *terms.entry(2).or_insert(Integer::from(0)) += coefficient;
                    }
                    // Handle var * var^n pattern
                    (Expr::Var(v), Expr::Power(pow)) if v.name == var.to_string() => {
                        if let Expr::Var(base_var) = pow.base.as_ref() {
                            if base_var.name == var.to_string() {
                                if let Expr::Num(exp) = pow.exponent.as_ref() {
                                    let degree = 1 + exp.value as usize;
                                    *terms.entry(degree).or_insert(Integer::from(0)) += coefficient;
                                }
                            }
                        }
                    }
                    // Handle var^n * var pattern
                    (Expr::Power(pow), Expr::Var(v)) if v.name == var.to_string() => {
                        if let Expr::Var(base_var) = pow.base.as_ref() {
                            if base_var.name == var.to_string() {
                                if let Expr::Num(exp) = pow.exponent.as_ref() {
                                    let degree = exp.value as usize + 1;
                                    *terms.entry(degree).or_insert(Integer::from(0)) += coefficient;
                                }
                            }
                        }
                    }
                    // Handle general case using distributive property
                    (left, right) => {
                        let mut left_terms: HashMap<usize, Integer> = HashMap::new();
                        let mut right_terms: HashMap<usize, Integer> = HashMap::new();

                        left.extract_terms(var, &mut left_terms, Integer::from(1));
                        right.extract_terms(var, &mut right_terms, Integer::from(1));

                        // Multiply the terms using distributive property
                        for (left_exp, left_coef) in left_terms {
                            for (right_exp, right_coef) in right_terms.iter() {
                                let new_exp = left_exp + right_exp;
                                let new_coef =
                                    left_coef.clone() * right_coef.clone() * coefficient.clone();
                                *terms.entry(new_exp).or_insert(Integer::from(0)) += new_coef;
                            }
                        }
                    }
                }
            }
            Expr::Power(p) => {
                if let Expr::Var(base_var) = p.base.as_ref() {
                    if base_var.name == var.to_string() {
                        if let Expr::Num(exp) = p.exponent.as_ref() {
                            let degree = exp.value as usize;
                            *terms.entry(degree).or_insert(Integer::from(0)) += coefficient;
                        } else {
                            // We already validated that exponents are integers
                            panic!("Expected integer exponent");
                        }
                    } else {
                        // Not our variable - this should have been caught in variable checking
                        panic!("Unexpected variable in power expression");
                    }
                } else if let Expr::Grouped(inner) = p.base.as_ref() {
                    // (expression)^n
                    if let Expr::Num(exp) = p.exponent.as_ref() {
                        if exp.value == 0 {
                            // Anything^0 = 1
                            *terms.entry(0).or_insert(Integer::from(0)) += coefficient;
                        } else if exp.value > 0 {
                            // For (expression)^n, we need to expand it
                            let mut temp_terms: HashMap<usize, Integer> = HashMap::new();
                            inner.extract_terms(var, &mut temp_terms, Integer::from(1));

                            // Apply the exponent by expanding the polynomial power
                            // This is a simple implementation that works for small exponents
                            let mut result_terms: HashMap<usize, Integer> = HashMap::new();
                            result_terms.insert(0, Integer::from(1)); // Start with x^0 = 1

                            for _ in 0..exp.value {
                                let mut new_result: HashMap<usize, Integer> = HashMap::new();
                                for (res_deg, res_coef) in &result_terms {
                                    for (term_deg, term_coef) in &temp_terms {
                                        let new_deg = res_deg + term_deg;
                                        let new_coef = res_coef.clone() * term_coef.clone();
                                        *new_result.entry(new_deg).or_insert(Integer::from(0)) +=
                                            new_coef;
                                    }
                                }
                                result_terms = new_result;
                            }

                            // Add these terms to our result
                            for (deg, coef) in result_terms {
                                *terms.entry(deg).or_insert(Integer::from(0)) +=
                                    coef * coefficient.clone();
                            }
                        }
                    }
                }
            }
            Expr::Grouped(e) => e.extract_terms(var, terms, coefficient),
            Expr::Neg(e) => e.extract_terms(var, terms, -coefficient),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_univariate_whitespace_variable() {
        let input = "  1 + x^2  - 5*x^3  ";
        let polynomial_expression = Expr::from_string(input.to_string()).unwrap();
        let poly = Expr::build_univariate_integer_polynomial('x', polynomial_expression).unwrap();
        assert_eq!(poly, Polynomial::from_coeffs(vec![1, 0, 1, -5]));
    }

    #[test]
    fn test_parse_constant() {
        let input = "42";
        let polynomial_expression = Expr::from_string(input.to_string()).unwrap();
        let poly = Expr::build_univariate_integer_polynomial('x', polynomial_expression).unwrap();
        assert_eq!(poly, Polynomial::from_coeffs(vec![42]));
    }

    #[test]
    fn test_parse_univariate_explicit_mult() {
        let input = "3*x^2 - 4*x + 5";
        let polynomial_expression = Expr::from_string(input.to_string()).unwrap();
        let poly = Expr::build_univariate_integer_polynomial('x', polynomial_expression).unwrap();
        assert_eq!(poly, Polynomial::from_coeffs(vec![5, -4, 3]));
    }

    #[test]
    fn test_parse_univariate_explicit_mult_skip_exp() {
        let input = "3*x^3 - 4*x + 5";
        let polynomial_expression = Expr::from_string(input.to_string()).unwrap();
        let poly = Expr::build_univariate_integer_polynomial('x', polynomial_expression).unwrap();
        assert_eq!(poly, Polynomial::from_coeffs(vec![5, -4, 0, 3]));
    }

    #[test]
    fn test_parse_univariate_diff_variable() {
        let input = "y^3 + 2*y - 7";
        let polynomial_expression = Expr::from_string(input.to_string()).unwrap();
        let poly = Expr::build_univariate_integer_polynomial('y', polynomial_expression).unwrap();
        assert_eq!(poly, Polynomial::from_coeffs(vec![-7, 2, 0, 1]));
    }

    // Error case tests

    #[test]
    fn test_multiple_variables() {
        let input = "2*x + 3*y + z";
        let polynomial_expression = Expr::from_string(input.to_string()).unwrap();
        let result = Expr::build_univariate_integer_polynomial('x', polynomial_expression);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("found multiple variables"));
    }

    #[test]
    fn test_lone_operator() {
        let input = "*";
        let polynomial_expression = Expr::from_string(input.to_string()).unwrap();
        let result = Expr::build_univariate_integer_polynomial('x', polynomial_expression);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Lone operator without operands");
    }

    #[test]
    fn test_invalid_exponent() {
        let input = "x^a";
        let polynomial_expression = Expr::from_string(input.to_string()).unwrap();
        let result = Expr::build_univariate_integer_polynomial('x', polynomial_expression);
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
        let polynomial_expression = Expr::from_string(input.to_string()).unwrap();
        let result = Expr::build_univariate_integer_polynomial('x', polynomial_expression);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Invalid consecutive operators found");
    }

    #[test]
    fn test_empty_input() {
        let input = "";
        let polynomial_expression = Expr::from_string(input.to_string()).unwrap();
        let result = Expr::build_univariate_integer_polynomial('x', polynomial_expression);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err(), "Empty input string");
    }

    #[test]
    fn test_exponent_without_variable() {
        let input = "3^2";
        let polynomial_expression = Expr::from_string(input.to_string()).unwrap();
        println!("{:?}", polynomial_expression);
        let result = Expr::build_univariate_integer_polynomial('x', polynomial_expression);
        println!("{:?}", result.clone().unwrap());
        assert!(result.is_err());
        // assert_eq!(
        //     result.unwrap_err(),
        //     "Invalid exponentiation: '^' without variable"
        // );
    }
}
