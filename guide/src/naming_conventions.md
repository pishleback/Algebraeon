# Naming Conventions

This section outlines the naming convention used for methods relating to the structure system. 

These naming convention are not settled and as such may not be used consistently. Nonetheless, this is the set of naming convention we should move towards using. We are open to suggestions on how the naming conventions could be better.

## By Value and By Reference

All methods should be prefixed by `into_` if they consume the structure type by value, and the `into_` prefix should be omitted if they take the structure type by reference. Where both an `into_` and non-`into_` version of a method exist, the resulting structures should behave the same.

## Constructing a Structure Type

When constructing a structure type, it is often the case that not all inputs allowed by the type system are valid. Structure types should implement a constructor of the form `new(..) -> Result<Self, String>` which returns `Ok` on valid inputs and a `String` explaining the invalidity on invalid inputs.

Structure types may also implement a constructor `new_unchecked(..) -> Self` which must panic on all inputs when compiled with `debug_assertions` enabled and which may result in undefined behaviour if called with an invalid input with `debug_assertions` disabled.

A nice pattern for implementing `new` and `new_unchecked` is to provide private methods `check(&self) -> Result<(), String>` and a `new_impl(..) -> Result<Self, String>`. They should be such that calling `new_impl` followed by `check` results in no errors if and only if the inputs were valid. Then `new` calls `new_impl` followed by `check` while propagating any errors, and `new_unchecked` calls `new_impl` with `.unwrap()`, and panics if `check()` fails only when `debug_assertions` are enabled.

## Obtaining a New Structure Type in General

`<>_structure` may be used for methods on a structure type `A` which produce a new structure type `B` whenever `A` and `B` are structures on different things. It is also acceptable for methods of this type to not have the `_structure` suffix.

For example, `.free_module(n)` exists for any ring `R` and returns the free module of finite rank `n`.

## Obtaining a New Structure Type with the Same Set

`<>_restructure` should be used for methods on a structure type `A` which produce a new structure type `B` whenever `A` and `B` both implement `SetSignature` with the same `Set`, meaning both use the same type for `Set` and `is_element` is `Ok` for `A` if and only if it is `Ok` for `B` on every instance of `Set`.

For example, `.abelian_group_restructure()` exists for any order in an algebraic number field (a ring whose elements are represented by `Vec<Integer>`s of a fixed length `n`) and returns the underlying rank `n` free abelian group structure type whose elements are also represented by `Vec<Integer>`s of a fixed length `n`. The resulting structure type is the same as that returned by `.free_module(n)` on the ring of integers.

## Constructing Inbound Morphisms

`inbound_<>` should be used for methods on a structure type `A` which produce a `Morphism` from `B` to `A` for some other structure type `B`.

For example, `.inbound_principal_integer_map()` exists for any ring `R` and returns the unique ring homomorphism from the integers to `R`.

## Constructing Outbound Morphisms

`outbound_<>` should be used for methods on a structure type `A` which produce a `Morphism` from `A` to `B` for some other structure type `B`.
