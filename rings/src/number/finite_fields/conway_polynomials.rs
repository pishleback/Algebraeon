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
    fn from_file(content: String) -> Self {
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

        assert_eq!(chars.next().unwrap(), '[');
        loop {
            assert_eq!(chars.next().unwrap(), '\n');
            let d = chars.next().unwrap();
            if d == '0' {
                break;
            }
            assert_eq!(d, '[');
            let (d, p) = parse_nat_until(&mut chars);
            assert_eq!(d, ',');
            let (d, n) = parse_nat_until(&mut chars);
            assert_eq!(d, ',');
            assert_eq!(chars.next().unwrap(), '[');
            let mut coeffs = vec![];
            loop {
                let (d, c) = parse_nat_until(&mut chars);
                coeffs.push(c);
                if d == ']' {
                    break;
                } else {
                    assert_eq!(d, ',');
                }
            }
            assert!(
                data.entry(p)
                    .or_insert(HashMap::new())
                    .insert(
                        n,
                        Polynomial::from_coeffs(
                            coeffs.into_iter().map(|x| Integer::from(x)).collect()
                        )
                    )
                    .is_none()
            );
            assert_eq!(chars.next().unwrap(), ']');
            assert_eq!(chars.next().unwrap(), ',');
        }

        Self { data }
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

static POLYNOMIAL_LOOKUP: OnceLock<ConwayPolynomialDatabase> = OnceLock::new();

fn get_polynomial_lookup() -> &'static ConwayPolynomialDatabase {
    #[allow(unreachable_code)]
    POLYNOMIAL_LOOKUP.get_or_init(|| {
        #[cfg(feature = "conway-polynomials-buildtime-fetch")]
        {
            return ConwayPolynomialDatabase::from_file(String::from(include_str!(concat!(
                env!("OUT_DIR"),
                "/conway_polynomials.txt"
            ))));
        }

        #[cfg(feature = "conway-polynomials-runtime-fetch")]
        {
            let url = include_str!(concat!(env!("OUT_DIR"), "/conway_polynomials_url.txt"));
            return ConwayPolynomialDatabase::from_file(
                reqwest::blocking::get(url).unwrap().text().unwrap(),
            );
        }

        ConwayPolynomialDatabase::empty()
    })
}

pub fn conway_polynomial(p: usize, n: usize) -> Result<&'static Polynomial<Integer>, ()> {
    match get_polynomial_lookup().get_polynomial(p, n) {
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
