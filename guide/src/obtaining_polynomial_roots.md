# Obtaining Polynomial Roots

This example finds all roots of a rational polynomial \\(f(x) \in \mathbb{Q}[x]\\) in various extension fields of \\(\mathbb{Q}\\). The key parts are
- Calling `.into_rational_extension()` on a field \\(K\\) to obtain the field extension \\(\mathbb{Q} \to K\\).
- Calling `.all_roots(&f)` on a field extension \\(\mathbb{Q} \to K\\) to obtain all roots of \\(f(x) \in \mathbb{Q}[x]\\) belonging to \\(K\\).

```rust
use algebraeon::nzq::{Natural, Rational};
use algebraeon::rings::isolated_algebraic::ComplexAlgebraic;
use algebraeon::rings::isolated_algebraic::PAdicAlgebraic;
use algebraeon::rings::isolated_algebraic::RealAlgebraic;
use algebraeon::rings::{polynomial::*, structure::*};
use algebraeon::sets::structure::*;

// Find all roots of f in some fields

let f = Polynomial::<Rational>::from_str(
    "(x - 3) * (x^2 - 17) * (x^2 + 1)", "x"
).unwrap();
println!("f = {}", f);

println!();

let two_adic = PAdicAlgebraic::structure(Natural::from(2u32))
    .into_inbound_principal_rational_map();
let three_adic = PAdicAlgebraic::structure(Natural::from(3u32))
    .into_inbound_principal_rational_map();
let complex = ComplexAlgebraic::structure()
    .into_inbound_principal_rational_map();
let real = RealAlgebraic::structure()
    .into_inbound_principal_rational_map();
let anf = Polynomial::from_str("x^2 + 1", "x")
    .unwrap()
    .algebraic_number_field()
    .unwrap()
    .into_inbound_principal_rational_map();

println!("Real roots of f");
for x in real.all_roots(&f) {
    println!("{}", x);
}

println!();

println!("Complex roots of f");
for x in complex.all_roots(&f) {
    println!("{}", x);
}

println!();

println!("2-adic roots of f");
for x in two_adic.all_roots(&f) {
    println!("{}", x);
}

println!();

println!("3-adic roots of f");
for x in three_adic.all_roots(&f) {
    println!("{}", x);
}

println!();

println!("Roots of f in Q[i]");
for x in anf.all_roots(&f) {
    println!(
        "{}",
        x.apply_map_into(MultiPolynomial::constant).evaluate(
            &Rational::structure()
                .into_multivariable_polynomial_ring()
                .var(Variable::new("i"))
        )
    );
}

/*
Output:
    f = (51)+(-17)λ+(48)λ^2+(-16)λ^3+(-3)λ^4+λ^5

    Real roots of f
    3
    -√17
    √17

    Complex roots of f
    3
    -i
    i
    -√17
    √17

    2-adic roots of f
    ...000011
    ...101001
    ...010111

    3-adic roots of f
    ...000010

    Roots of f in Q[i]
    (3)1
    i
    (-1)i
*/
```