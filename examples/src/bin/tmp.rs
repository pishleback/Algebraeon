#![allow(dead_code, warnings)]

fn main() {
    use algebraeon_rings::{polynomial::polynomial::*, structure::elements::*};
    use malachite_nz::{integer::Integer, natural::Natural};
    let x = Polynomial::<Integer>::var().into_ergonomic();
    let f = (x.pow(2) - 17).into_verbose();
    println!("{:?}", f);
    for mut root in f.all_padic_roots(&Natural::from(2u32)) {
        println!("{}", root.string_repr(20)); // Show 20 2-adic digits
    }
    /*
    Output:
        ...00110010011011101001
        ...11001101100100010111
    */

    // Truncating to the last 16 bits, it can be verified that the square of these values is 17 modulo 2^16 using machine arithmetic:
    let a = 0b0010011011101001u16;
    assert_eq!(a.wrapping_mul(a), 17u16);
    let b = 0b1101100100010111u16;
    assert_eq!(b.wrapping_mul(b), 17u16);
}
