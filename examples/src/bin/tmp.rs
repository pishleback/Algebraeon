#![allow(dead_code, warnings)]


fn main() {
    use algebraeon_rings::{
        polynomial::polynomial::*,
        structure::{elements::*, structure::*},
    };
    use malachite_nz::integer::Integer;
    
    let x = &Polynomial::<Integer>::var().into_ergonomic();
    let f = (x.pow(2) - 5*x + 6).into_verbose();
    println!("f(λ) = {}", f.factor().unwrap());
    /*
    Output:
        f(λ) = 1 * ((-2)+λ) * ((-3)+λ)
    */
    
    let f = (x.pow(15) - 1).into_verbose();
    println!("f(λ) = {}", f.factor().unwrap());
    /*
    Output:
        f(λ) = 1 * ((-1)+λ) * (1+λ+λ^2) * (1+λ+λ^2+λ^3+λ^4) * (1+(-1)λ+λ^3+(-1)λ^4+λ^5+(-1)λ^7+λ^8)
    */
}
