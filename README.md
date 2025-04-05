# Algebraeon

Algebraeon is a computer algebra system written purely in Rust. It implements algorithms for working with matrices, polynomials, algebraic numbers, factorizations, etc. The focus is on exact algebraic computations over approximate numerical solutions. Algebraeon is in early stages of development and the API subject to change. Algebraeon uses [Malachite](https://www.malachite.rs/) under the hood for arbitrary sized integer and rational numbers.

See the [user guide](https://pishleback.github.io/Algebraeon/) to get started.

Algebraeon is hosted on crates.io [here](https://crates.io/crates/algebraeon) and the formal documentation can be found [here](https://docs.rs/algebraeon/latest/algebraeon/).

# Contributing

Contributions are welcome. There are two primary ways to contribute:

## Using the issue tracker

Use the issue tracker if you have questions, feature requests, bugs reports, etc.

## Changing the code-base

You should fork this repository, make changes, and submit a pull request. Submitted code should, where applicable, have associated unit tests.

Algebraeon is organized as a [cargo workspace](https://doc.rust-lang.org/book/ch14-03-cargo-workspaces.html). Run `cargo test` in the root directory to build and run all tests.

A suggested workflow for testing new features:

- Create a new binary in `examples/src/bin`, for example `my_main.rs`.
- To run, use `cargo run --bin my_main` in the root directory.
- Test any changes to the codebase with unit tests and/or using `my_main.rs`.

## CLA

Anyone who contributes code is required to sign the CLA. You can sign the CLA when you submit a pull request. The CLA allows for us to relicense future versions of Algebraeon, including your contribution, under any licences the Free Software Foundation classifies as Free Software Licence and which are approved by the Open Source Initiative as Open Source licences. It does not allow us to relicense of your contribution under any other more restrictive licences.
