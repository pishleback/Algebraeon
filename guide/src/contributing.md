# Contributing

Contributions are welcome and should be made via GitHub.

## Using the issue tracker

Use the issue tracker if you have questions, feature requests, bugs reports, etc.

## What to do?

If you're unsure what you can do:
 - Have a look through the issue tracker for unassigned issues and leave a comment.
 - Ask in the discord server.
 - Look at the [project ideas](project_ideas.md).

## Development Environment

Algebraeon is organized as a [cargo workspace](https://doc.rust-lang.org/book/ch14-03-cargo-workspaces.html). Run `cargo test` in the root directory to build and run all tests.

A suggested workflow for getting started:
- Checkout the repository.
- Create a new binary in `examples/src/bin`, for example `my_main.rs`.
- Copy an example into `my_main.rs`.
- Run it with `cargo run --bin my_main` in the root directory.

To submit changes for inclusion in the main branch: Fork this repository, make changes, and submit a pull request.

## Use of AI

Use of AI is fine in contributions subjet to the following requests:
 - Please outline in a PR comment the extent to which AI has been used.
 - Ensure every AI generated line has been understood by you, a human.
 - Don't use AI to help write code _and_ tests for the same code.

## CLA

If you wish to contribute code we ask that you sign a CLA. The CLA allows us to relicense future versions of Algebraeon, including your contribution, under any licences the Free Software Foundation classifies as a Free Software Licence and which are approved by the Open Source Initiative as Open Source licences. It does not allow us to relicense of your contribution under any other more restrictive licences. The CLA can be signed on GitHub once you submit a pull request.
