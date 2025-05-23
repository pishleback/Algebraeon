use crate::polynomial::{MultiPolynomial, Polynomial, Variable};
use crate::structure::{MetaRing, MetaSemiRing};
use algebraeon_nzq::*;
use lalrpop_util::lalrpop_mod;
use std::collections::HashMap;
use std::fmt;

lalrpop_mod!(polynomial_parser, "/parsing/polynomial_grammar.rs"); // synthesized by LALRPOP

#[derive(Debug, Clone, PartialEq)]
struct ParseVar {
    name: String,
}

impl ParseVar {
    fn formatted(&self) -> String {
        if self.name.len() > 1 {
            format!("{{{}}}", self.name) // Add braces for multi-character variables
        } else {
            self.name.clone()
        }
    }
}

impl fmt::Display for ParseVar {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.formatted())
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Number {
    numerator: Integer,
    denominator: Integer,
}

impl Number {
    fn new(numerator: Integer, denominator: Integer) -> Self {
        if denominator == Integer::from(0) {
            panic!("Denominator cannot be zero");
        }
        Number {
            numerator,
            denominator,
        }
    }
}

impl fmt::Display for Number {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.denominator == Integer::from(1) {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "{}/{}", self.numerator, self.denominator)
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
enum Expr {
    Var(ParseVar),
    Num(Number),
    Sum(Sum),
    Product(Product),
    Power(Power),
    Grouped(Box<Expr>),
    Neg(Box<Expr>),
}

impl Expr {
    fn validate(&self) -> Result<(), String> {
        match self {
            Expr::Power(p) => p.validate(),
            Expr::Product(p) => p.validate(),
            Expr::Sum(s) => s.validate(),
            Expr::Neg(e) => e.validate(),
            Expr::Grouped(e) => e.validate(),
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
struct Sum {
    left: Box<Expr>,
    right: Box<Expr>,
}

impl Sum {
    fn validate(&self) -> Result<(), String> {
        self.left.validate()?;
        self.right.validate()
    }
}

impl fmt::Display for Sum {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.right.as_ref() {
            Expr::Neg(inner) => write!(f, "{} - {}", self.left, inner),
            _ => write!(f, "{} + {}", self.left, self.right),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Product {
    left: Box<Expr>,
    right: Box<Expr>,
}

impl Product {
    fn validate(&self) -> Result<(), String> {
        if let Expr::Power(p) = self.right.as_ref() {
            if let Expr::Num(n) = p.exponent.as_ref() {
                if n.denominator == Integer::from(1) && n.numerator == Integer::from(-1) {
                    return Err("Division not allowed in polynomials".to_string());
                }
            }
        }
        self.left.validate()?;
        self.right.validate()
    }
}

impl fmt::Display for Product {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} * {}", self.left, self.right)
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Power {
    base: Box<Expr>,
    exponent: Box<Expr>,
}

impl Power {
    fn validate(&self) -> Result<(), String> {
        match self.exponent.as_ref() {
            Expr::Num(n) => {
                if n.denominator != Integer::from(1) {
                    return Err("Fractional exponents not allowed in polynomials".to_string());
                }
                if n.numerator < Integer::from(0) {
                    return Err("Negative exponents not allowed in polynomials".to_string());
                }
            }
            _ => return Err("Exponents must be integer constants in polynomials".to_string()),
        }
        self.base.validate()
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
}

impl Expr {
    // Convert expression into univariate integer polynomial
    fn build_univariate_integer_polynomial(
        expression: &Expr,
        var: &str,
    ) -> Result<Polynomial<Integer>, String> {
        // First validate it's a valid polynomial
        expression.validate()?;

        // Collect terms as (exponent, coefficient) pairs
        let mut terms: HashMap<usize, Integer> = HashMap::new();

        // Check expression contains only our specified variable or is constant
        let vars = expression.collect_variables();
        if vars.is_empty() {
            // Constant polynomial
        } else if vars.len() > 1 {
            return Err(format!(
                "Expected univariate polynomial with variable '{}', but found multiple variables: {:?}",
                var, vars
            ));
        } else if !vars.contains(var) {
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

    // New function: Convert expression into univariate rational polynomial
    fn build_univariate_rational_polynomial(
        expression: &Expr,
        var: &str,
    ) -> Result<Polynomial<Rational>, String> {
        // First validate it's a valid polynomial
        expression.validate()?;

        // Collect terms as (exponent, coefficient) pairs
        let mut terms: HashMap<usize, Rational> = HashMap::new();

        // Check expression contains only our specified variable or is constant
        let vars = expression.collect_variables();
        if vars.is_empty() {
            // Constant polynomial
        } else if vars.len() > 1 {
            return Err(format!(
                "Expected univariate polynomial with variable '{}', but found multiple variables: {:?}",
                var, vars
            ));
        } else if !vars.contains(var) {
            return Err(format!(
                "Expected polynomial with variable '{}', but found '{}'",
                var,
                vars.iter().next().unwrap()
            ));
        }

        // Extract the terms
        expression.extract_rational_terms(var, &mut terms, Rational::from(1));

        // Find the highest degree
        let max_degree = terms.keys().max().copied().unwrap_or(0);

        // Create coefficient vector
        let mut coeffs = vec![Rational::from(0); max_degree + 1];
        for (degree, coeff) in terms {
            coeffs[degree] = coeff;
        }

        Ok(Polynomial::from_coeffs(coeffs))
    }

    // Convert expression into multivariate integer polynomial
    fn build_multivariate_integer_polynomial(
        expression: &Expr,
        variable_mapping: HashMap<&str, Variable>,
    ) -> Result<MultiPolynomial<Integer>, String> {
        // First validate it's a valid polynomial
        expression.validate()?;

        // Build the multivariate polynomial recursively
        expression.to_multivariate_integer(&variable_mapping)
    }

    fn to_multivariate_integer(
        &self,
        variable_mapping: &HashMap<&str, Variable>,
    ) -> Result<MultiPolynomial<Integer>, String> {
        match self {
            Expr::Var(v) => match variable_mapping.get(v.name.as_str()) {
                Some(var) => Ok(MultiPolynomial::<Integer>::var(var.clone())),
                None => Err(format!(
                    "Variable '{}' not found in variable mapping",
                    v.name
                )),
            },
            Expr::Num(n) => {
                if n.denominator != Integer::from(1) {
                    return Err("Non-integer coefficient in integer polynomial".to_string());
                }
                Ok(MultiPolynomial::<Integer>::constant(n.numerator.clone()))
            }
            Expr::Sum(s) => {
                let left = s.left.to_multivariate_integer(variable_mapping)?;
                let right = s.right.to_multivariate_integer(variable_mapping)?;
                Ok(MultiPolynomial::add(&left, &right))
            }
            Expr::Product(p) => {
                let left = p.left.to_multivariate_integer(variable_mapping)?;
                let right = p.right.to_multivariate_integer(variable_mapping)?;
                Ok(MultiPolynomial::mul(&left, &right))
            }
            Expr::Power(p) => {
                let base = p.base.to_multivariate_integer(variable_mapping)?;
                match p.exponent.as_ref() {
                    Expr::Num(n) => {
                        if n.denominator != Integer::from(1) {
                            return Err(
                                "Fractional exponents not allowed in polynomials".to_string()
                            );
                        }
                        if n.numerator < Integer::from(0) {
                            return Err("Negative exponents not allowed in polynomials".to_string());
                        }
                        // Convert Integer to Natural for nat_pow
                        let exp_usize: usize = (&n.numerator).try_into().expect("Exponent out of range for usize conversion");
                        let exp_natural = Natural::from(exp_usize);
                        Ok(base.nat_pow(&exp_natural))
                    }
                    _ => Err("Exponents must be integer constants in polynomials".to_string()),
                }
            }
            Expr::Grouped(e) => e.to_multivariate_integer(variable_mapping),
            Expr::Neg(e) => {
                let inner = e.to_multivariate_integer(variable_mapping)?;
                Ok(MultiPolynomial::neg(&inner))
            }
        }
    }

    // Convert expression into multivariate rational polynomial
    fn build_multivariate_rational_polynomial(
        expression: &Expr,
        variable_mapping: HashMap<&str, Variable>,
    ) -> Result<MultiPolynomial<Rational>, String> {
        // First validate it's a valid polynomial
        expression.validate()?;

        // Build the multivariate polynomial recursively
        expression.to_multivariate_rational(&variable_mapping)
    }

    // Convert this expression to a multivariate rational polynomial
    fn to_multivariate_rational(
        &self,
        variable_mapping: &HashMap<&str, Variable>,
    ) -> Result<MultiPolynomial<Rational>, String> {
        match self {
            Expr::Var(v) => match variable_mapping.get(v.name.as_str()) {
                Some(var) => Ok(MultiPolynomial::<Rational>::var(var.clone())),
                None => Err(format!(
                    "Variable '{}' not found in variable mapping",
                    v.name
                )),
            },
            Expr::Num(n) => {
                let rational_coeff =
                    Rational::from_integers(n.numerator.clone(), n.denominator.clone());
                Ok(MultiPolynomial::<Rational>::constant(rational_coeff))
            }
            Expr::Sum(s) => {
                let left = s.left.to_multivariate_rational(variable_mapping)?;
                let right = s.right.to_multivariate_rational(variable_mapping)?;
                Ok(MultiPolynomial::add(&left, &right))
            }
            Expr::Product(p) => {
                let left = p.left.to_multivariate_rational(variable_mapping)?;
                let right = p.right.to_multivariate_rational(variable_mapping)?;
                Ok(MultiPolynomial::mul(&left, &right))
            }
            Expr::Power(p) => {
                let base = p.base.to_multivariate_rational(variable_mapping)?;
                match p.exponent.as_ref() {
                    Expr::Num(n) => {
                        if n.denominator != Integer::from(1) {
                            return Err(
                                "Fractional exponents not allowed in polynomials".to_string()
                            );
                        }
                        if n.numerator < Integer::from(0) {
                            return Err("Negative exponents not allowed in polynomials".to_string());
                        }
                        // Convert Integer to Natural for nat_pow
                        let exp_usize: usize = (&n.numerator).try_into().expect("Exponent out of range for usize conversion");
                        let exp_natural = Natural::from(exp_usize);
                        Ok(base.nat_pow(&exp_natural))
                    }
                    _ => Err("Exponents must be integer constants in polynomials".to_string()),
                }
            }
            Expr::Grouped(e) => e.to_multivariate_rational(variable_mapping),
            Expr::Neg(e) => {
                let inner = e.to_multivariate_rational(variable_mapping)?;
                Ok(MultiPolynomial::neg(&inner))
            }
        }
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

    // Recursively extract terms from the expression with integer coefficients
    fn extract_terms(&self, var: &str, terms: &mut HashMap<usize, Integer>, coefficient: Integer) {
        match self {
            Expr::Var(v) => {
                if v.name == var {
                    // x^1
                    *terms.entry(1).or_insert(Integer::from(0)) += coefficient;
                } else {
                    // Should never happen due to variable checking
                    panic!("Unexpected variable: {}", v.name);
                }
            }
            Expr::Num(n) => {
                // Constant term
                let num_val = n.numerator.clone() * coefficient.clone();
                let denom_val = n.denominator.clone();

                // Handle fractions by ensuring denominator=1 for integer polynomials
                if denom_val != Integer::from(1) {
                    if num_val.clone() % denom_val.clone() != Integer::from(0) {
                        panic!("Non-integer coefficient in integer polynomial");
                    }
                    *terms.entry(0).or_insert(Integer::from(0)) += num_val / denom_val;
                } else {
                    *terms.entry(0).or_insert(Integer::from(0)) += num_val;
                }
            }
            Expr::Sum(s) => {
                s.left.extract_terms(var, terms, coefficient.clone());
                s.right.extract_terms(var, terms, coefficient);
            }
            Expr::Product(p) => {
                match (p.left.as_ref(), p.right.as_ref()) {
                    // Handle coefficient * var^n pattern
                    (Expr::Num(n), right) => {
                        let new_coeff = coefficient * n.numerator.clone() / n.denominator.clone();
                        right.extract_terms(var, terms, new_coeff);
                    }
                    // Handle var^n * coefficient pattern
                    (left, Expr::Num(n)) => {
                        let new_coeff = coefficient * n.numerator.clone() / n.denominator.clone();
                        left.extract_terms(var, terms, new_coeff);
                    }
                    // Handle var * var pattern (produces var^2)
                    (Expr::Var(v1), Expr::Var(v2)) if v1.name == var && v2.name == var => {
                        *terms.entry(2).or_insert(Integer::from(0)) += coefficient;
                    }
                    // Handle var * var^n pattern
                    (Expr::Var(v), Expr::Power(pow)) if v.name == var => {
                        if let Expr::Var(base_var) = pow.base.as_ref() {
                            if base_var.name == var {
                                if let Expr::Num(exp) = pow.exponent.as_ref() {
                                    if exp.denominator == Integer::from(1) {
                                        // Convert the Integer to usize safely
                                        let exp_usize: usize = (&exp.numerator).try_into()
                                            .expect("Exponent out of range for usize conversion");
                                        let degree = 1 + exp_usize;
                                        *terms.entry(degree).or_insert(Integer::from(0)) +=
                                            coefficient;
                                    }
                                }
                            }
                        }
                    }
                    // Handle var^n * var pattern
                    (Expr::Power(pow), Expr::Var(v)) if v.name == var => {
                        if let Expr::Var(base_var) = pow.base.as_ref() {
                            if base_var.name == var {
                                if let Expr::Num(exp) = pow.exponent.as_ref() {
                                    if exp.denominator == Integer::from(1) {
                                        // Convert the Integer to usize safely
                                        let exp_usize: usize = (&exp.numerator).try_into()
                                            .expect("Exponent out of range for usize conversion");
                                        let degree = 1 + exp_usize;
                                        *terms.entry(degree).or_insert(Integer::from(0)) +=
                                            coefficient;
                                    }
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
                    if base_var.name == var {
                        if let Expr::Num(exp) = p.exponent.as_ref() {
                            if exp.denominator == Integer::from(1) {
                                // Convert the Integer to usize safely
                                let degree: usize = (&exp.numerator).try_into()
                                    .expect("Exponent out of range for usize conversion");
                                *terms.entry(degree).or_insert(Integer::from(0)) += coefficient;
                            } else {
                                panic!("Fractional exponent in integer polynomial");
                            }
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
                        if exp.denominator == Integer::from(1) {
                            if exp.numerator == Integer::from(0) {
                                // Anything^0 = 1
                                *terms.entry(0).or_insert(Integer::from(0)) += coefficient;
                            } else if exp.numerator > Integer::from(0) {
                                // For (expression)^n, we need to expand it
                                let mut temp_terms: HashMap<usize, Integer> = HashMap::new();
                                inner.extract_terms(var, &mut temp_terms, Integer::from(1));

                                // Apply the exponent by expanding the polynomial power
                                // This is a simple implementation that works for small exponents
                                let mut result_terms: HashMap<usize, Integer> = HashMap::new();
                                result_terms.insert(0, Integer::from(1)); // Start with x^0 = 1

                                // Convert the Integer to usize safely
                                let exp_usize: usize = (&exp.numerator).try_into()
                                    .expect("Exponent out of range for usize conversion");

                                for _ in 0..exp_usize {
                                    let mut new_result: HashMap<usize, Integer> = HashMap::new();
                                    for (res_deg, res_coef) in &result_terms {
                                        for (term_deg, term_coef) in &temp_terms {
                                            let new_deg = res_deg + term_deg;
                                            let new_coef = res_coef.clone() * term_coef.clone();
                                            *new_result
                                                .entry(new_deg)
                                                .or_insert(Integer::from(0)) += new_coef;
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
                        } else {
                            panic!("Fractional exponent in integer polynomial");
                        }
                    }
                }
            }
            Expr::Grouped(e) => e.extract_terms(var, terms, coefficient),
            Expr::Neg(e) => e.extract_terms(var, terms, -coefficient),
        }
    }

    // Method for rational polynomial extraction
    fn extract_rational_terms(
        &self,
        var: &str,
        terms: &mut HashMap<usize, Rational>,
        coefficient: Rational,
    ) {
        match self {
            Expr::Var(v) => {
                if v.name == var {
                    // x^1
                    *terms.entry(1).or_insert(Rational::from(0)) += coefficient;
                } else {
                    // Should never happen due to variable checking
                    panic!("Unexpected variable: {}", v.name);
                }
            }
            Expr::Num(n) => {
                // Constant term
                let num_rational =
                    Rational::from_integers(n.numerator.clone(), n.denominator.clone());
                *terms.entry(0).or_insert(Rational::from(0)) += coefficient * num_rational;
            }
            Expr::Sum(s) => {
                s.left
                    .extract_rational_terms(var, terms, coefficient.clone());
                s.right.extract_rational_terms(var, terms, coefficient);
            }
            Expr::Product(p) => {
                match (p.left.as_ref(), p.right.as_ref()) {
                    // Handle coefficient * var^n pattern
                    (Expr::Num(n), right) => {
                        let num_rational =
                            Rational::from_integers(n.numerator.clone(), n.denominator.clone());
                        let new_coeff = coefficient * num_rational;
                        right.extract_rational_terms(var, terms, new_coeff);
                    }
                    // Handle var^n * coefficient pattern
                    (left, Expr::Num(n)) => {
                        let num_rational =
                            Rational::from_integers(n.numerator.clone(), n.denominator.clone());
                        let new_coeff = coefficient * num_rational;
                        left.extract_rational_terms(var, terms, new_coeff);
                    }
                    // Handle var * var pattern (produces var^2)
                    (Expr::Var(v1), Expr::Var(v2)) if v1.name == var && v2.name == var => {
                        *terms.entry(2).or_insert(Rational::from(0)) += coefficient;
                    }
                    // Handle var * var^n pattern
                    (Expr::Var(v), Expr::Power(pow)) if v.name == var => {
                        if let Expr::Var(base_var) = pow.base.as_ref() {
                            if base_var.name == var {
                                if let Expr::Num(exp) = pow.exponent.as_ref() {
                                    if exp.denominator == Integer::from(1) {
                                        // Convert the Integer to usize safely
                                        let exp_usize: usize = (&exp.numerator).try_into()
                                            .expect("Exponent out of range for usize conversion");
                                        let degree = 1 + exp_usize;
                                        *terms.entry(degree).or_insert(Rational::from(0)) +=
                                            coefficient;
                                    }
                                }
                            }
                        }
                    }
                    // Handle var^n * var pattern
                    (Expr::Power(pow), Expr::Var(v)) if v.name == var => {
                        if let Expr::Var(base_var) = pow.base.as_ref() {
                            if base_var.name == var {
                                if let Expr::Num(exp) = pow.exponent.as_ref() {
                                    if exp.denominator == Integer::from(1) {
                                        // Convert the Integer to usize safely
                                        let exp_usize: usize = (&exp.numerator).try_into()
                                            .expect("Exponent out of range for usize conversion");
                                        let degree = 1 + exp_usize;
                                        *terms.entry(degree).or_insert(Rational::from(0)) +=
                                            coefficient;
                                    }
                                }
                            }
                        }
                    }
                    // Handle general case using distributive property
                    (left, right) => {
                        let mut left_terms: HashMap<usize, Rational> = HashMap::new();
                        let mut right_terms: HashMap<usize, Rational> = HashMap::new();

                        left.extract_rational_terms(var, &mut left_terms, Rational::from(1));
                        right.extract_rational_terms(var, &mut right_terms, Rational::from(1));

                        // Multiply the terms using distributive property
                        for (left_exp, left_coef) in left_terms {
                            for (right_exp, right_coef) in right_terms.iter() {
                                let new_exp = left_exp + right_exp;
                                let new_coef =
                                    left_coef.clone() * right_coef.clone() * coefficient.clone();
                                *terms.entry(new_exp).or_insert(Rational::from(0)) += new_coef;
                            }
                        }
                    }
                }
            }
            Expr::Power(p) => {
                if let Expr::Var(base_var) = p.base.as_ref() {
                    if base_var.name == var {
                        if let Expr::Num(exp) = p.exponent.as_ref() {
                            if exp.denominator == Integer::from(1) {
                                // Convert the Integer to usize safely
                                let degree: usize = (&exp.numerator).try_into()
                                    .expect("Exponent out of range for usize conversion");
                                *terms.entry(degree).or_insert(Rational::from(0)) += coefficient;
                            } else {
                                panic!("Fractional exponent in polynomial");
                            }
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
                        if exp.denominator == Integer::from(1) {
                            if exp.numerator == Integer::from(0) {
                                // Anything^0 = 1
                                *terms.entry(0).or_insert(Rational::from(0)) += coefficient;
                            } else if exp.numerator > Integer::from(0) {
                                // For (expression)^n, we need to expand it
                                let mut temp_terms: HashMap<usize, Rational> = HashMap::new();
                                inner.extract_rational_terms(
                                    var,
                                    &mut temp_terms,
                                    Rational::from(1),
                                );

                                // Apply the exponent by expanding the polynomial power
                                let mut result_terms: HashMap<usize, Rational> = HashMap::new();
                                result_terms.insert(0, Rational::from(1)); // Start with x^0 = 1

                                let exp_usize: usize = (&exp.numerator).try_into()
                                    .expect("Exponent out of range for usize conversion");

                                for _ in 0..exp_usize {
                                    let mut new_result: HashMap<usize, Rational> = HashMap::new();
                                    for (res_deg, res_coef) in &result_terms {
                                        for (term_deg, term_coef) in &temp_terms {
                                            let new_deg = res_deg + term_deg;
                                            let new_coef = res_coef.clone() * term_coef.clone();
                                            *new_result
                                                .entry(new_deg)
                                                .or_insert(Rational::from(0)) += new_coef;
                                        }
                                    }
                                    result_terms = new_result;
                                }

                                // Add these terms to our result
                                for (deg, coef) in result_terms {
                                    *terms.entry(deg).or_insert(Rational::from(0)) +=
                                        coef * coefficient.clone();
                                }
                            }
                        } else {
                            panic!("Fractional exponent in polynomial");
                        }
                    }
                }
            }
            Expr::Grouped(e) => e.extract_rational_terms(var, terms, coefficient),
            Expr::Neg(e) => e.extract_rational_terms(var, terms, -coefficient),
        }
    }
}

/// Create a polynomial with integer coefficients from a string
pub fn parse_integer_polynomial(polynomial_str: &str, var: &str) -> Result<Polynomial<Integer>, String> {
    match polynomial_parser::ExprParser::new().parse(polynomial_str) {
        Ok(expr) => Expr::build_univariate_integer_polynomial(&expr, var),
        Err(e) => Err(format!("Failed to parse expression: {:?}", e)),
    }
}

/// Create a polynomial with rational coefficients from a string
pub fn parse_rational_polynomial(
    polynomial_str: &str,
    var: &str,
) -> Result<Polynomial<Rational>, String> {
    match polynomial_parser::ExprParser::new().parse(polynomial_str) {
        Ok(expr) => Expr::build_univariate_rational_polynomial(&expr, var),
        Err(e) => Err(format!("Failed to parse expression: {:?}", e)),
    }
}

/// Create a multivariate polynomial with integer coefficients from a string
pub fn parse_multivariate_integer_polynomial(
    polynomial_str: &str,
    variable_mapping: HashMap<&str, Variable>,
) -> Result<MultiPolynomial<Integer>, String> {
    match polynomial_parser::ExprParser::new().parse(polynomial_str) {
        Ok(expr) => Expr::build_multivariate_integer_polynomial(&expr, variable_mapping),
        Err(e) => Err(format!("Failed to parse expression: {:?}", e)),
    }
}

/// Create a multivariate polynomial with integer coefficients from a string
pub fn parse_multivariate_rational_polynomial(
    polynomial_str: &str,
    variable_mapping: HashMap<&str, Variable>,
) -> Result<MultiPolynomial<Rational>, String> {
    match polynomial_parser::ExprParser::new().parse(polynomial_str) {
        Ok(expr) => Expr::build_multivariate_rational_polynomial(&expr, variable_mapping),
        Err(e) => Err(format!("Failed to parse expression: {:?}", e)),
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use super::*;

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

    // TODO: this test case needs a potential grammar update (shelved for now)
    // #[test]
    // fn test_rational_polynomial_divide_polynomial_by_integer() {
    //     // dividing by integers should be ok for polynomials with rational coefficients
    //     // though dividing by polynomial expressions should not be valid
    //
    //     let result = parse_and_build_rational_poly("(x^2 + 3*x + 2)/3", "x").unwrap();
    //
    //     assert_eq!(
    //         result,
    //         Polynomial::from_coeffs(vec![
    //             Rational::from_str("2/3").unwrap(),
    //             Rational::from_str("1").unwrap(),
    //             Rational::from_str("1/3").unwrap()
    //         ])
    //     );
    // }

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

        let result =
            parse_multivariate_integer_polynomial("x^2 + y^3", variable_mapping).unwrap();

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

        let result = parse_multivariate_integer_polynomial(
            "3*x^2*y + 2*x*y^2 - x + 5",
            variable_mapping,
        )
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

        let result =
            parse_multivariate_rational_polynomial("1/2*x + 3/4*y", variable_mapping).unwrap();

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

        let result =
            parse_multivariate_integer_polynomial("(x + y)^2", variable_mapping).unwrap();

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

        let result =
            parse_multivariate_integer_polynomial("(x + y + z)^30", variable_mapping).unwrap();

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

        let result =
            parse_multivariate_integer_polynomial("{x} + y", variable_mapping).unwrap();

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
                Integer::from_str(
                    "123456789123456789123456789123456789123456789123456789123456789"
                )
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
                Rational::from_str(
                    "-123456789123456789123456789123456789123456789123456789123456789"
                )
                .unwrap(),
                Rational::from_str(
                    "123456789123456789123456789123456789123456789123456789123456789"
                )
                .unwrap(),
                Rational::from_str(
                    "123456789123456789123456789123456789123456789123456789123456789"
                )
                .unwrap()
            ])
        );
    }
}
