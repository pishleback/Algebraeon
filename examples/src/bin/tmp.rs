#![allow(dead_code, warnings)]

use algebraeon_rings::{number::algebraic::isolated_roots::padic::*, polynomial::polynomial::Polynomial, structure::{elements::*, structure::*}};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

fn main() {
    // let ring = PAdicAlgebraicStructure::new(Natural::from(5u32));

    // let x = Polynomial::<Integer>::var().into_ergonomic();

    // let a = {
    //     let f = (x.pow(3) - 3 * x.pow(2) - x.pow(1) + 1).into_verbose();
    //     let r = f.all_padic_roots(&Natural::from(5u32));
    //     assert_eq!(r.len(), 1);
    //     r.into_iter().next().unwrap().shift_by(-1)
    // };

    // let b = {
    //     let f = (x.pow(4) + x.pow(2) - 2 * x.pow(1) - 1).into_verbose();
    //     let r = f.all_padic_roots(&Natural::from(5u32));
    //     assert_eq!(r.len(), 1);
    //     r.into_iter().next().unwrap().shift_by(-1)
    // };

    // let c = {
    //     let f = (x.pow(5) + x.pow(2) + 2 * x.pow(1) + 1).into_verbose();
    //     let r = f.all_padic_roots(&Natural::from(5u32));
    //     assert_eq!(r.len(), 1);
    //     r.into_iter().next().unwrap().shift_by(-1)
    // };

    // let d = a.clone().shift_by(4);

    // let x = ring
    //     .from_rat(&Rational::from_integers(
    //         Integer::from(2),
    //         Integer::from(3 * 125),
    //     ))
    //     .unwrap();

    // println!("a = {}", a);
    // println!("b = {}", b);
    // println!("c = {}", c);
    // println!("d = {}", d);
    // println!("x = {}", x);

    // println!("-a = {}", ring.neg(&a));
    // debug_assert_eq!(
    //     ring.neg(&a).reduce_modulo_valuation(5).digits(),
    //     (
    //         vec![
    //             Natural::from(3u8),
    //             Natural::from(0u8),
    //             Natural::from(2u8),
    //             Natural::from(3u8),
    //             Natural::from(3u8),
    //             Natural::from(1u8)
    //         ],
    //         -1
    //     )
    // );
    // debug_assert_eq!(a.valuation(), Some(-1));

    // println!("-b = {}", ring.neg(&b));
    // debug_assert_eq!(
    //     ring.neg(&b).reduce_modulo_valuation(5).digits(),
    //     (
    //         vec![
    //             Natural::from(3u8),
    //             Natural::from(1u8),
    //             Natural::from(3u8),
    //             Natural::from(2u8),
    //             Natural::from(4u8),
    //             Natural::from(3u8)
    //         ],
    //         -1
    //     )
    // );

    // println!("-x = {}", ring.neg(&x));
    // debug_assert_eq!(
    //     ring.neg(&x).reduce_modulo_valuation(3).digits(),
    //     (
    //         vec![
    //             Natural::from(1u8),
    //             Natural::from(3u8),
    //             Natural::from(1u8),
    //             Natural::from(3u8),
    //             Natural::from(1u8),
    //             Natural::from(3u8),
    //         ],
    //         -3
    //     )
    // );

    // println!("a+x = {}", ring.add(&a, &x));
    // debug_assert_eq!(
    //     ring.add(&a, &x).reduce_modulo_valuation(6).digits(),
    //     (
    //         vec![
    //             Natural::from(4u8),
    //             Natural::from(1u8),
    //             Natural::from(0u8),
    //             Natural::from(1u8),
    //             Natural::from(1u8),
    //             Natural::from(3u8),
    //             Natural::from(4u8),
    //             Natural::from(4u8),
    //             Natural::from(1u8),
    //         ],
    //         -3
    //     )
    // );

    // println!("b+x = {}", ring.add(&b, &x));
    // debug_assert_eq!(
    //     ring.add(&b, &x).reduce_modulo_valuation(6).digits(),
    //     (
    //         vec![
    //             Natural::from(4u8),
    //             Natural::from(1u8),
    //             Natural::from(0u8),
    //             Natural::from(0u8),
    //             Natural::from(0u8),
    //             Natural::from(4u8),
    //             Natural::from(3u8),
    //             Natural::from(2u8),
    //             Natural::from(2u8),
    //         ],
    //         -3
    //     )
    // );

    // println!("c+x = {}", ring.add(&c, &x));
    // debug_assert_eq!(
    //     ring.add(&c, &x).reduce_modulo_valuation(6).digits(),
    //     (
    //         vec![
    //             Natural::from(4u8),
    //             Natural::from(1u8),
    //             Natural::from(4u8),
    //             Natural::from(2u8),
    //             Natural::from(1u8),
    //             Natural::from(1u8),
    //             Natural::from(0u8),
    //             Natural::from(2u8),
    //             Natural::from(0u8),
    //         ],
    //         -3
    //     )
    // );

    // println!("a+b = {}", ring.add(&a, &b));
    // debug_assert_eq!(
    //     ring.add(&a, &b).reduce_modulo_valuation(6).digits(),
    //     (
    //         vec![
    //             Natural::from(4u8),
    //             Natural::from(2u8),
    //             Natural::from(4u8),
    //             Natural::from(3u8),
    //             Natural::from(1u8),
    //             Natural::from(4u8),
    //             Natural::from(2u8),
    //         ],
    //         -1
    //     )
    // );

    // println!("c+d = {}", ring.add(&c, &d));
    // // debug_assert_eq!(
    // //     ring.add(&a, &b).reduce_modulo_valuation(6).digits(),
    // //     (
    // //         vec![
    // //             Natural::from(4u8),
    // //             Natural::from(2u8),
    // //             Natural::from(4u8),
    // //             Natural::from(3u8),
    // //             Natural::from(1u8),
    // //             Natural::from(4u8),
    // //             Natural::from(2u8),
    // //         ],
    // //         -1
    // //     )
    // // );
}
