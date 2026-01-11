# Set Structures

## Elements

 - `Set` for an object representing a set of elements.
 - `Eq : Set` for a set with a binary predicate called equality and denoted `=`. It's required to be an equivalence relation on the set.

## Ordering

 - `PartialOrd : Set + Eq` for sets with a partial ordering.
 - `Ord : Set + PartialOrd` for sets with a total ordering.