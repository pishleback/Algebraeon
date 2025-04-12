use crate::polynomial::Polynomial;
use algebraeon_nzq::{Integer, Natural};
use std::{collections::HashMap, sync::OnceLock};

pub struct ConwayPolynomialDatabase {
    data: HashMap<Natural, HashMap<Natural, Polynomial<Integer>>>,
}

static POLYNOMIAL_LOOKUP: OnceLock<ConwayPolynomialDatabase> = OnceLock::new();

fn get_polynomial_lookup() -> &'static ConwayPolynomialDatabase {
    POLYNOMIAL_LOOKUP.get_or_init(|| {
        let content: String = {
            #[cfg(feature = "runtime-fetch")]
            {
                println!("online");
                let url = include_str!(concat!(env!("OUT_DIR"), "/conway_polynomials.txt"));
                let content = reqwest::blocking::get(url).unwrap().text().unwrap();
                content
            }

            #[cfg(not(feature = "runtime-fetch"))]
            {
                println!("offline");
                String::from(include_str!(concat!(
                    env!("OUT_DIR"),
                    "/conway_polynomials.txt"
                )))
            }
        };

        println!("{:?}", content);

        ConwayPolynomialDatabase {
            data: HashMap::new(),
        }
    })
}

impl ConwayPolynomialDatabase {
    fn get_polynomial(&self, p: &Natural, n: &Natural) -> Option<&Polynomial<Integer>> {
        self.data.get(p)?.get(n)
    }
}

pub fn conway_polynomial(p: &Natural, n: &Natural) -> Result<&'static Polynomial<Integer>, ()> {
    match get_polynomial_lookup().get_polynomial(p, n) {
        Some(p) => Ok(p),
        None => Err(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::Polynomial;

    #[test]
    fn test_conway_polynomial() {
        let p = conway_polynomial(&Natural::from(2u32), &Natural::from(6u32)).unwrap();
        println!("p = {}", p);
        assert_eq!(p, &Polynomial::from_coeffs(vec![1, 1, 0, 1, 1, 0, 1]));
    }
}
