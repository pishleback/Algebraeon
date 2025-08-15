use crate::parsing::ast::*;
use crate::polynomial::{MultiPolynomial, Polynomial, Variable};
use crate::structure::{MetaRing, MetaSemiRing};
use algebraeon_nzq::*;
use lalrpop_util::lalrpop_mod;
use std::collections::HashMap;

lalrpop_mod!(polynomial_parser, "/parsing/polynomial_grammar.rs");

impl Expr {
    pub fn from_string(input: String) -> Result<Self, String> {
        match polynomial_parser::ExprParser::new().parse(&input) {
            Ok(expr) => Ok(*expr),
            Err(err) => Err(format!("Failed to parse expression: {}", err)),
        }
    }

    // Convert expression into univariate integer polynomial
    pub fn build_univariate_integer_polynomial(
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
    pub fn build_univariate_rational_polynomial(
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
    pub fn build_multivariate_integer_polynomial(
        expression: &Expr,
        variable_mapping: HashMap<&str, Variable>,
    ) -> Result<MultiPolynomial<Integer>, String> {
        // First validate it's a valid polynomial
        expression.validate()?;

        // Build the multivariate polynomial recursively
        expression.to_multivariate_integer(&variable_mapping)
    }

    pub fn to_multivariate_integer(
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
                if s.terms.is_empty() {
                    return Ok(MultiPolynomial::<Integer>::constant(Integer::from(0)));
                }

                // Process the first term
                let first_term = &s.terms[0];
                let mut result = first_term.term.to_multivariate_integer(variable_mapping)?;
                if !first_term.sign {
                    result = MultiPolynomial::neg(&result);
                }

                // Add/subtract the remaining terms
                for term in &s.terms[1..] {
                    let term_poly = term.term.to_multivariate_integer(variable_mapping)?;
                    if term.sign {
                        result = MultiPolynomial::add(&result, &term_poly);
                    } else {
                        result = MultiPolynomial::add(&result, &MultiPolynomial::neg(&term_poly));
                    }
                }

                Ok(result)
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
                        let exp_usize: usize = (&n.numerator)
                            .try_into()
                            .expect("Exponent out of range for usize conversion");
                        let exp_natural = Natural::from(exp_usize);
                        Ok(base.nat_pow(&exp_natural))
                    }
                    _ => Err("Exponents must be integer constants in polynomials".to_string()),
                }
            }
            Expr::Grouped(e) => e.to_multivariate_integer(variable_mapping),
        }
    }

    // Convert expression into multivariate rational polynomial
    pub fn build_multivariate_rational_polynomial(
        expression: &Expr,
        variable_mapping: HashMap<&str, Variable>,
    ) -> Result<MultiPolynomial<Rational>, String> {
        // First validate it's a valid polynomial
        expression.validate()?;

        // Build the multivariate polynomial recursively
        expression.to_multivariate_rational(&variable_mapping)
    }

    // Convert this expression to a multivariate rational polynomial
    pub fn to_multivariate_rational(
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
                if s.terms.is_empty() {
                    return Ok(MultiPolynomial::<Rational>::constant(Rational::from(0)));
                }

                // Process the first term
                let first_term = &s.terms[0];
                let mut result = first_term.term.to_multivariate_rational(variable_mapping)?;
                if !first_term.sign {
                    result = MultiPolynomial::neg(&result);
                }

                // Add/subtract the remaining terms
                for term in &s.terms[1..] {
                    let term_poly = term.term.to_multivariate_rational(variable_mapping)?;
                    if term.sign {
                        result = MultiPolynomial::add(&result, &term_poly);
                    } else {
                        result = MultiPolynomial::add(&result, &MultiPolynomial::neg(&term_poly));
                    }
                }

                Ok(result)
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
                        let exp_usize: usize = (&n.numerator)
                            .try_into()
                            .expect("Exponent out of range for usize conversion");
                        let exp_natural = Natural::from(exp_usize);
                        Ok(base.nat_pow(&exp_natural))
                    }
                    _ => Err("Exponents must be integer constants in polynomials".to_string()),
                }
            }
            Expr::Grouped(e) => e.to_multivariate_rational(variable_mapping),
        }
    }

    // Recursively extract terms from the expression with integer coefficients
    pub fn extract_terms(
        &self,
        var: &str,
        terms: &mut HashMap<usize, Integer>,
        coefficient: Integer,
    ) {
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
                if denom_val == Integer::from(1) {
                    *terms.entry(0).or_insert(Integer::from(0)) += num_val;
                } else {
                    assert!(
                        num_val.clone() % denom_val.clone() == Integer::from(0),
                        "Non-integer coefficient in integer polynomial"
                    );
                    *terms.entry(0).or_insert(Integer::from(0)) += num_val / denom_val;
                }
            }
            Expr::Sum(s) => {
                for (i, term) in s.terms.iter().enumerate() {
                    let term_coeff = if i == 0 {
                        // First term uses the coefficient directly
                        if term.sign {
                            coefficient.clone()
                        } else {
                            -coefficient.clone()
                        }
                    } else {
                        // Subsequent terms are added/subtracted
                        if term.sign {
                            coefficient.clone()
                        } else {
                            -coefficient.clone()
                        }
                    };
                    term.term.extract_terms(var, terms, term_coeff);
                }
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
                        if let Expr::Var(base_var) = pow.base.as_ref()
                            && base_var.name == var
                            && let Expr::Num(exp) = pow.exponent.as_ref()
                            && exp.denominator == Integer::from(1)
                        {
                            // Convert the Integer to usize safely
                            let exp_usize: usize = (&exp.numerator)
                                .try_into()
                                .expect("Exponent out of range for usize conversion");
                            let degree = 1 + exp_usize;
                            *terms.entry(degree).or_insert(Integer::from(0)) += coefficient;
                        }
                    }
                    // Handle var^n * var pattern
                    (Expr::Power(pow), Expr::Var(v)) if v.name == var => {
                        if let Expr::Var(base_var) = pow.base.as_ref()
                            && base_var.name == var
                            && let Expr::Num(exp) = pow.exponent.as_ref()
                            && exp.denominator == Integer::from(1)
                        {
                            // Convert the Integer to usize safely
                            let exp_usize: usize = (&exp.numerator)
                                .try_into()
                                .expect("Exponent out of range for usize conversion");
                            let degree = 1 + exp_usize;
                            *terms.entry(degree).or_insert(Integer::from(0)) += coefficient;
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
                            for (right_exp, right_coef) in &right_terms {
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
                                let degree: usize = (&exp.numerator)
                                    .try_into()
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
                        #[allow(clippy::comparison_chain)]
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
                                let exp_usize: usize = (&exp.numerator)
                                    .try_into()
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
        }
    }

    // Method for rational polynomial extraction
    pub fn extract_rational_terms(
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
                for (i, term) in s.terms.iter().enumerate() {
                    let term_coeff = if i == 0 {
                        // First term uses the coefficient directly
                        if term.sign {
                            coefficient.clone()
                        } else {
                            -coefficient.clone()
                        }
                    } else {
                        // Subsequent terms are added/subtracted
                        if term.sign {
                            coefficient.clone()
                        } else {
                            -coefficient.clone()
                        }
                    };
                    term.term.extract_rational_terms(var, terms, term_coeff);
                }
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
                        if let Expr::Var(base_var) = pow.base.as_ref()
                            && base_var.name == var
                            && let Expr::Num(exp) = pow.exponent.as_ref()
                            && exp.denominator == Integer::from(1)
                        {
                            // Convert the Integer to usize safely
                            let exp_usize: usize = (&exp.numerator)
                                .try_into()
                                .expect("Exponent out of range for usize conversion");
                            let degree = 1 + exp_usize;
                            *terms.entry(degree).or_insert(Rational::from(0)) += coefficient;
                        }
                    }
                    // Handle var^n * var pattern
                    (Expr::Power(pow), Expr::Var(v)) if v.name == var => {
                        if let Expr::Var(base_var) = pow.base.as_ref()
                            && base_var.name == var
                            && let Expr::Num(exp) = pow.exponent.as_ref()
                            && exp.denominator == Integer::from(1)
                        {
                            // Convert the Integer to usize safely
                            let exp_usize: usize = (&exp.numerator)
                                .try_into()
                                .expect("Exponent out of range for usize conversion");
                            let degree = 1 + exp_usize;
                            *terms.entry(degree).or_insert(Rational::from(0)) += coefficient;
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
                            for (right_exp, right_coef) in &right_terms {
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
                                let degree: usize = (&exp.numerator)
                                    .try_into()
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
                            #[allow(clippy::comparison_chain)]
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

                                let exp_usize: usize = (&exp.numerator)
                                    .try_into()
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
        }
    }
}

/// Create a polynomial with integer coefficients from a string
pub fn parse_integer_polynomial(
    polynomial_str: &str,
    var: &str,
) -> Result<Polynomial<Integer>, String> {
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
