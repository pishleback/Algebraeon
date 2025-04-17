#![allow(dead_code, warnings, unused)]

use std::marker::PhantomData;
use std::str::FromStr;
use std::task::Poll;

use algebraeon_groups::examples::symmetric::*;
use algebraeon_groups::structure::*;
use algebraeon_rings::linear::matrix::*;
use algebraeon_rings::polynomial::*;
use algebraeon_rings::rings::finite_fields::modulo::*;
use algebraeon_rings::structure::*;

fn main() {
    let a = algebraeon_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
        let x: usize = x.into();
        if x == 23 {
            23.into()
        } else {
            let x: Modulo<23> = x.into();
            let x: usize = (x.into_ergonomic() + 1).into_verbose().into();
            x.into()
        }
    });

    let a = algebraeon_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
        let x: usize = x.into();
        if x == 23 {
            23.into()
        } else {
            let x: Modulo<23> = x.into();
            let x: usize = (x.into_ergonomic() + 1).into_verbose().into();
            x.into()
        }
    })
    .unwrap();
    let b = algebraeon_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
        let x: usize = x.into();
        if x == 23 {
            23.into()
        } else {
            let x: Modulo<23> = x.into();
            let x: usize = (x.into_ergonomic() * 2).into_verbose().into();
            x.into()
        }
    })
    .unwrap();
    let c = algebraeon_combinatorics::modular_permutations::modular_permutation::<24>(|x| {
        let x: usize = x.into();
        if x == 23 {
            0.into()
        } else if x == 0 {
            23.into()
        } else {
            let x: Modulo<23> = x.into();
            let x: usize = x.inv().unwrap().neg().into();
            x.into()
        }
    })
    .unwrap();
    let d = Permutation::<24>::new_from_cycles(vec![
        vec![14, 17, 11, 19, 22],
        vec![20, 10, 7, 5, 21],
        vec![18, 4, 2, 6, 1],
        vec![8, 16, 13, 9, 12],
    ])
    .unwrap();

    println!("a = {}", a);
    println!("b = {}", b);
    println!("c = {}", c);
    println!("d = {}", d);

    let sg = Permutation::generated_finite_subgroup(vec![a, b, c, d]);
    println!("{:?}", sg.size())
}
