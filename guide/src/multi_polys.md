# Multivariate Polynomials

## Example

```rust 
use algebraeon::nzq::Integer;
use algebraeon::rings::polynomial::*;
use algebraeon::rings::structure::*;
use std::collections::HashMap;

// Define symbols for the variables we wish to use
let x_var = Variable::new("x");
let y_var = Variable::new("y");

// Construct the polynomials x and y from the variables
let x = &MultiPolynomial::<Integer>::var(x_var.clone());
let y = &MultiPolynomial::<Integer>::var(y_var.clone());

// x and y can now be used just like other ring elements in Algebraeon
let f = MultiPolynomial::add(x, y).nat_pow(&3u32.into());

// f = x^3+(3)x^2y+(3)xy^2+y^3
println!("f = {}", f);

// For evaluating f the inputs are specified using the variable symbols
let v = f.evaluate(HashMap::from([(x_var, &1.into()), (y_var, &2.into())]));

// v = f(1, 2) = 27
println!("v = {}", v);
```