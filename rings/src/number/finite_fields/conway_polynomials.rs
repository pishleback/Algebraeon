use crate::polynomial::Polynomial;
use algebraeon_nzq::{Integer, Natural};
use std::{cell::OnceCell, collections::HashMap, sync::OnceLock};

pub struct ConwayPolynomialDatabase {
    data: HashMap<Natural, HashMap<Natural, Polynomial<Integer>>>,
}

static POLYNOMIAL_LOOKUP: OnceLock<ConwayPolynomialDatabase> = OnceLock::new();

fn get_polynomial_lookup() -> &'static ConwayPolynomialDatabase {
    POLYNOMIAL_LOOKUP.get_or_init(|| {
        #[cfg(feature = "runtime-fetch")]
        {
            println!("online");
        }

        #[cfg(not(feature = "runtime-fetch"))]
        {
            println!("offline");
        }

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

pub fn get_polynomial(p: &Natural, n: &Natural) -> Option<&'static Polynomial<Integer>> {
    get_polynomial_lookup().get_polynomial(p, n)
}
