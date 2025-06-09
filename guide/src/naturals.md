# Natural Numbers

## Constructing naturals

There are several ways to construct naturals. Some of them are:

- `Natural::ZERO`, `Natural::ONE`, `Natural::TWO`
- from a string: `Natural::from_str()`
- from builtin unsigned types `Natural::from(4u32)`

```rust
use std::str::FromStr;
use algebraeon::nzq::Natural;

let one = Natural::ONE;
let two = Natural::TWO;
let n = Natural::from(42usize);
let big = Natural::from_str("706000565581575429997696139445280900").unwrap();
```

## Basic operations

Natural supports the following operators:

- `+` (addition)

- `*` (multiplication)

- `%` (modulo)

For exponentiation, use the method `.pow(&exp)` instead of `^` (which is xor).

## Available functions

- [`choose`](https://docs.rs/algebraeon-nzq/latest/algebraeon_nzq/fn.choose.html)
- `euler_totient`
- [`factorial`](https://docs.rs/algebraeon-nzq/latest/algebraeon_nzq/struct.Natural.html#method.factorial)
- [`gcd`](https://docs.rs/algebraeon-nzq/latest/algebraeon_nzq/fn.gcd.html)
- `is_prime`
- `is_square`
- [`lcm`](https://docs.rs/algebraeon-nzq/latest/algebraeon_nzq/fn.lcm.html)
- [`nth_root_floor`](https://docs.rs/algebraeon-nzq/latest/algebraeon_nzq/struct.Natural.html#method.nth_root_floor)
- [`nth_root_ceil`](https://docs.rs/algebraeon-nzq/latest/algebraeon_nzq/struct.Natural.html#method.nth_root_ceil)
- [`sqrt_ceil`](https://docs.rs/algebraeon-nzq/latest/algebraeon_nzq/struct.Natural.html#method.sqrt_ceil)
- [`sqrt_floor`](https://docs.rs/algebraeon-nzq/latest/algebraeon_nzq/struct.Natural.html#method.sqrt_floor)
- `sqrt_if_square`

```rust
use algebraeon::nzq::*;
use algebraeon::rings::natural::factorization::NaturalCanonicalFactorizationStructure;
use algebraeon::sets::structure::*;

let a = Natural::from(12u32);
let b = Natural::from(5u32);

// Basic operations
let sum = &a + &b; // 17
let product = &a * &b; // 60
let power = a.pow(&b); // 248832
let modulo = &a % &b; // 2

assert_eq!(sum, Natural::from(17u32));
assert_eq!(product, Natural::from(60u32));
assert_eq!(power, Natural::from(248832u32));
assert_eq!(modulo, Natural::from(2u32));

// Factorial
assert_eq!(b.factorial(), Natural::from(120u32));

// nth_root_floor and nth_root_ceil
assert_eq!(a.nth_root_floor(&b), Natural::from(1u32));
assert_eq!(a.nth_root_ceil(&b), Natural::from(2u32));

// sqrt_floor, sqrt_ceil, sqrt_if_square
assert_eq!(a.sqrt_floor(), Natural::from(3u32));
assert_eq!(a.sqrt_ceil(), Natural::from(4u32));

let square = Natural::from(144u32); // 12²
assert_eq!(square.sqrt_if_square(), Some(a.clone()));
assert_eq!(a.sqrt_if_square(), None);

// is_square
assert!(square.is_square());
assert!(!a.is_square());

// Combinatorics
assert_eq!(choose(&a, &b), Natural::from(792u32));

// GCD and LCM
assert_eq!(gcd(a.clone(), b.clone()), Natural::from(1u32));
assert_eq!(lcm(a.clone(), b.clone()), Natural::from(60u32));

// is_prime
use algebraeon::rings::natural::factorization::primes::is_prime;
assert!(!is_prime(&a)); // 12 is not prime
assert!(is_prime(&b)); // 5 is prime

// Euler's totient function
use algebraeon::rings::natural::factorization::factor;
assert_eq!(
    Natural::structure()
        .factorizations()
        .euler_totient(&factor(a).unwrap()),
    Natural::from(4u32)
); // φ(12) = 4
assert_eq!(
    Natural::structure()
        .factorizations()
        .euler_totient(&factor(b).unwrap()),
    Natural::from(4u32)
); // φ(5) = 4
```

## Factoring

Algebraeon implements [Lenstra elliptic-curve factorization](https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization) for quickly finding prime factors up to around 20 digits.

```rust
# use algebraeon::sets::structure::ToStringSignature;
# use algebraeon::{nzq::Natural, rings::natural::factorization::factor};
# use algebraeon::{
    rings::natural::factorization::NaturalCanonicalFactorizationStructure,
    sets::structure::MetaType,
};
# use std::str::FromStr;
# 
let n = Natural::from_str("706000565581575429997696139445280900").unwrap();
let f = factor(n.clone()).unwrap();
println!(
    "{} = {}",
    n,
    Natural::structure().factorizations().to_string(&f)
);;
/*
Output:
    706000565581575429997696139445280900 = 2^2 × 5^2 × 6988699669998001 × 1010203040506070809
*/
```
