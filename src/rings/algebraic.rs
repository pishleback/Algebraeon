#![allow(dead_code)]

use malachite_nz::integer::Integer;

use super::nzq::*;
use super::poly::*;
use super::ring::ComRing;
use super::ring::*;

fn root_sum_poly(p: Polynomial<Integer>, q: Polynomial<Integer>) {
    // x, z = sympy.symbols("x, z")
    // p = poly_to_sympy(p, x)
    // q = poly_to_sympy(q, x)
    // r = sympy.Poly(q.compose(sympy.Poly(z - x)), x, domain = "ZZ[z]")
    // sum_poly = sympy.Poly(sympy.resultant(p, r)).sqf_part()
    // return sympy_to_poly(sum_poly, z)
    todo!();
}

fn root_prod_poly(p: Polynomial<Integer>, q: Polynomial<Integer>) {
    // x, t = sympy.symbols("x, t")
    // p = poly_to_sympy(p, t)
    // q = poly_to_sympy(q, x)
    // r = sympy.Poly(q.homogenize(t), t, domain = "ZZ[x]")
    // #x ** q.degree() * q(t * x ** -1)
    // prod_poly = sympy.Poly(sympy.resultant(p, r), x, domain = "ZZ").sqf_part()
    // return sympy_to_poly(prod_poly, x)
    todo!();
}

#[cfg(test)]
mod tests {
    //#[test]
    fn test_root_sum_poly() {
        todo!();
    }

    //#[test]
    fn test_root_prod_poly() {
        todo!();
    }
}
