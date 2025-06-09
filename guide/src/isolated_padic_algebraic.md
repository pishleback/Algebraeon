# p-Adic Algebraic Numbers

Isolated algebraic numbers in the \\(p\\)-adic fields \\(\mathbb{Q}_p\\).

## Taking Square Roots

It is possible to take square roots of \\(p\\)-adic algebraic numbers by calling `.square_roots`.

```rust
use algebraeon::nzq::{Natural, Rational};
use algebraeon::rings::isolated_algebraic::PAdicAlgebraic;
use algebraeon::rings::{polynomial::*, structure::*};

let two_adic = PAdicAlgebraic::structure(Natural::from(2u32));
for x in two_adic
    .rational_extension()
    .all_roots(&Polynomial::<Rational>::from_str("(x - 3) * (x^2 - 17)", "x").unwrap())
{
    println!("{} is_square={}", x, two_adic.is_square(&x));
    if let Some((r1, r2)) = two_adic.square_roots(&x) {
        println!("    sqrt = {} and {}", r1, r2);
    }
}

/*
Output:
    ...000011 is_square=false
    ...101001 is_square=true
        sqrt = ...101101 and ...010011
    ...010111 is_square=false
*/
```