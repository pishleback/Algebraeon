use algebraeon_nzq::*;
use std::collections::HashSet;
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct ParseVar {
    pub name: String,
}

impl ParseVar {
    pub fn formatted(&self) -> String {
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
pub struct Number {
    pub numerator: Integer,
    pub denominator: Integer,
}

impl Number {
    pub fn new(numerator: Integer, denominator: Integer) -> Self {
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
pub enum Expr {
    Var(ParseVar),
    Num(Number),
    Sum(Sum),
    Product(Product),
    Power(Power),
    Grouped(Box<Expr>),
    Neg(Box<Expr>),
}

impl Expr {
    pub fn validate(&self) -> Result<(), String> {
        match self {
            Expr::Power(p) => p.validate(),
            Expr::Product(p) => p.validate(),
            Expr::Sum(s) => s.validate(),
            Expr::Neg(e) => e.validate(),
            Expr::Grouped(e) => e.validate(),
            _ => Ok(()),
        }
    }

    // Helper method to collect all variables in the expression
    pub fn collect_variables(&self) -> HashSet<String> {
        let mut vars = HashSet::new();
        self.add_variables(&mut vars);
        vars
    }

    pub fn add_variables(&self, vars: &mut HashSet<String>) {
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
    pub fn validate(&self) -> Result<(), String> {
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
pub struct Product {
    pub left: Box<Expr>,
    pub right: Box<Expr>,
}

impl Product {
    pub fn validate(&self) -> Result<(), String> {
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
pub struct Power {
    pub base: Box<Expr>,
    pub exponent: Box<Expr>,
}

impl Power {
    pub fn validate(&self) -> Result<(), String> {
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
