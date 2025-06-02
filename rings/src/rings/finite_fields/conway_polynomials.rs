use crate::polynomial::Polynomial;
use algebraeon_nzq::Integer;
use std::{
    collections::HashMap,
    str::{Chars, FromStr},
    sync::OnceLock,
};

pub struct ConwayPolynomialDatabase {
    data: HashMap<usize, HashMap<usize, Polynomial<Integer>>>,
}

impl ConwayPolynomialDatabase {
    fn from_file(content: String) -> Result<Self, ()> {
        let content = content.replace("\r\n", "\n");

        let content = &content[24..];
        let mut chars = content.chars();

        let parse_nat_until = |chars: &mut Chars<'_>| -> (char, usize) {
            let mut nat_string = String::new();
            loop {
                let char = chars.next().unwrap();
                if char == ',' || char == '[' || char == ']' {
                    return (char, usize::from_str(&nat_string).unwrap());
                }
                nat_string.push(char);
            }
        };

        let mut data = HashMap::new();

        macro_rules! assert_or_err {
            ($a:expr) => {
                if !$a {
                    return Err(());
                }
            };
        }

        macro_rules! assert_eq_or_err {
            ($a:expr, $b:expr) => {
                if $a != $b {
                    return Err(());
                }
            };
        }

        assert_eq_or_err!(chars.next().unwrap(), '[');
        loop {
            assert_eq_or_err!(chars.next().unwrap(), '\n');
            let d = chars.next().unwrap();
            if d == '0' {
                break;
            }
            assert_eq_or_err!(d, '[');
            let (d, p) = parse_nat_until(&mut chars);
            assert_eq_or_err!(d, ',');
            let (d, n) = parse_nat_until(&mut chars);
            assert_eq_or_err!(d, ',');
            assert_eq_or_err!(chars.next().unwrap(), '[');
            let mut coeffs = vec![];
            loop {
                let (d, c) = parse_nat_until(&mut chars);
                coeffs.push(c);
                if d == ']' {
                    break;
                }
                assert_eq_or_err!(d, ',');
            }
            assert_or_err!(
                data.entry(p)
                    .or_insert(HashMap::new())
                    .insert(
                        n,
                        Polynomial::from_coeffs(coeffs.into_iter().map(Integer::from).collect())
                    )
                    .is_none()
            );
            assert_eq_or_err!(chars.next().unwrap(), ']');
            assert_eq_or_err!(chars.next().unwrap(), ',');
        }

        Ok(Self { data })
    }

    fn empty() -> Self {
        Self {
            data: HashMap::new(),
        }
    }

    fn get_polynomial(&self, p: usize, n: usize) -> Option<&Polynomial<Integer>> {
        self.data.get(&p)?.get(&n)
    }
}

static POLYNOMIAL_LOOKUP: OnceLock<Result<ConwayPolynomialDatabase, ()>> = OnceLock::new();

fn get_polynomial_lookup() -> &'static Result<ConwayPolynomialDatabase, ()> {
    #[allow(unreachable_code)]
    POLYNOMIAL_LOOKUP.get_or_init(|| {
        ConwayPolynomialDatabase::from_file(String::from(include_str!("./conway_polynomials.txt")))
    })
}

pub fn conway_polynomial(p: usize, n: usize) -> Result<&'static Polynomial<Integer>, ()> {
    match get_polynomial_lookup()
        .as_ref()
        .unwrap_or_else(|_err| panic!("Failed to parse Conway polynomial file. Please report this error to `https://github.com/pishleback/Algebraeon/issues`."))
        .get_polynomial(p, n)
    {
        Some(p) => Ok(p),
        None => Err(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{polynomial::Polynomial, structure::IntoErgonomic};

    #[test]
    fn test_conway_polynomial() {
        let x = &Polynomial::<Integer>::var().into_ergonomic();
        assert_eq!(conway_polynomial(2, 1).unwrap(), &(x + 1).into_verbose());
        assert_eq!(
            conway_polynomial(2, 2).unwrap(),
            &(x.pow(2) + x + 1).into_verbose()
        );
        assert_eq!(
            conway_polynomial(2, 3).unwrap(),
            &(x.pow(3) + x + 1).into_verbose()
        );
        assert_eq!(
            conway_polynomial(2, 6).unwrap(),
            &(x.pow(6) + x.pow(4) + x.pow(3) + x + 1).into_verbose()
        );

        assert_eq!(conway_polynomial(3, 1).unwrap(), &(x + 1).into_verbose());
        assert_eq!(
            conway_polynomial(3, 2).unwrap(),
            &(x.pow(2) + 2 * x + 2).into_verbose()
        );
        assert_eq!(
            conway_polynomial(3, 3).unwrap(),
            &(x.pow(3) + 2 * x + 1).into_verbose()
        );
    }
}
