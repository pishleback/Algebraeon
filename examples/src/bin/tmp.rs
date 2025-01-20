#![allow(dead_code, warnings)]

fn main() {
    use algebraeon_rings::polynomial::{
        multipoly::{MultiPolynomial, Variable},
        polynomial::Polynomial,
    };
    use malachite_nz::integer::Integer;

    let a_var = Variable::new("a");
    let b_var = Variable::new("b");
    let c_var = Variable::new("c");
    let d_var = Variable::new("d");

    let a = MultiPolynomial::<Integer>::var(a_var);
    let b = MultiPolynomial::<Integer>::var(b_var);
    let c = MultiPolynomial::<Integer>::var(c_var);
    let d = MultiPolynomial::<Integer>::var(d_var);

    let p =
        Polynomial::<MultiPolynomial<Integer>>::from_coeffs(vec![c.clone(), b.clone(), a.clone()]);
    println!("p(λ) = {}", p);
    println!("disc(p) = {}", p.discriminant().unwrap());

    println!();

    let p = Polynomial::<MultiPolynomial<Integer>>::from_coeffs(vec![
        d.clone(),
        c.clone(),
        b.clone(),
        a.clone(),
    ]);
    println!("p(λ) = {}", p);
    println!("disc(p) = {}", p.discriminant().unwrap());
}
