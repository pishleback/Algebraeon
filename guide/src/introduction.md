# Introduction
**What is Algebraeon?**

Algebraeon is a computer algebra system written purely in Rust. It implements algorithms for working with matrices, polynomials, algebraic numbers, factorizations, etc. The general focus is on number theory and exact algebraic computations over approximate numerical solutions. A growing subset of the features are available as a [Python module](https://pishleback.github.io/Algebraeon-Python/).

> [!WARNING]
> The API is subject to breaking changes at this time. The hope is to stabilize things more in the future. 

That said, each release should be fully functional and bug-free. Please report any issues [here](https://github.com/pishleback/Algebraeon/issues).

**Source Code and Contributing**

The project is open source and hosted on GitHub [here](https://github.com/pishleback/Algebraeon). See the [contributing](contributing.md) section for more information on how to contribute.

**Links**

 - [GitHub](https://github.com/pishleback/Algebraeon)
 - [Crates.io](https://crates.io/crates/algebraeon)
 - [Algebraeon as a Python module](https://pishleback.github.io/Algebraeon-Python/)

**Crates**

Algebraeon is published to crates.io under an umbrella crate [algebraeon](https://crates.io/crates/algebraeon) and re-exports the sub-crates:

- [algebraeon-sets](https://crates.io/crates/algebraeon-sets)
- [algebraeon-nzq](https://crates.io/crates/algebraeon-nzq)
- [algebraeon-groups](https://crates.io/crates/algebraeon-groups)
- [algebraeon-rings](https://crates.io/crates/algebraeon-rings)
- [algebraeon-geometry](https://crates.io/crates/algebraeon-geometry)

**A Taster**

To give a flavour of what it's about, here is a non-exhaustive list of the things Algebraeon can already do:

- Factoring and primality testing of integers.
- Factoring polynomials defined over finite fields, the integers, the rational numbers, algebraic number fields.
- Polynomial root isolation and arithmetic in the real line, complex plane, and \\(p\\)-adics.
- Linear algebra; matrix image and kernel, solving linear systems, LLL lattice basis rediction, Jordan normal form, Hermite normal form, Smith normal form, boolean operations with linear subspaces.
- Algebraic number fields; computing the ring of integers, working with ideals, factoring ideals.
- Coset enumeration for finitely presented groups.
- Multiplication tables for small finite groups.