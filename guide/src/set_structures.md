# Set Structures

Let \\(X\\) be a set. This section outlines the traits Algebraeon provides for defining elements of \\(X\\) and comparing elements with respect to equality or an ordering on \\(X\\).

## Elements and Equality

 - `Set` for an object representing a set of elements \\(X\\). Set the associated type `Set` to a type `X` to use instances of `X` to represent elements of \\(X\\). It's not necessary that _all_ instances of `X` represent valid elements of \\(X\\). The method `validate_element(a: X) -> Result` should return `Ok` when `a` represents a valid element of \\(X\\) and `Err` when it does not. All other operations with elements of \\(X\\) may exhibit undefined behaviour (ideally panic immediately - at least when building debug mode).
 - `Eq : Set` for a set with a binary predicate \\(=\\) called equality. \\(=\\) is required to be an equivalence relation, meaning
 \\[a = a \quad \forall a \in X\\]
 \\[a = b \implies b = a \quad \forall a, b \in X\\]
 \\[a = b \quad \text{and} \quad b = c \implies a = c \quad \forall a, b \in X\\]
 The method `equal(a: X, b: X) -> bool` indicates whether `a` and `b` are equal.

## Orderings

The ordering structures make use of Rusts `std::cmp::Ordering` enum to relay information about how elements of a set compare with respect to an ordering.

 - `PartialOrd : Set + Eq` for sets with a binary predicate \\(\le\\) satisfying the axioms for a partial order.
\\[a \le a \quad \forall a \in X\\]
\\[a \le b \quad \text{and} \quad b \le a \implies a = b \quad \forall a, b \in X\\]
\\[a \le b \quad \text{and} \quad b \le c \implies a \le c \quad \forall a, b \in X\\]
 The method `partial_ord(a: X, b: X) -> Option<Ordering>` returns information about how `a` and `b` compare wtih respect to \\(=\\) and \\(\le\\).
 - `Ord : Set + PartialOrd` for sets with a total ordering, so in addition to the axioms for a partial ordering, \\(\le\\) also satisfies
\\[a \le b \quad \text{or} \quad b \le a \quad \forall a, b \in X\\]
 The method `ord(a: X, b: X) -> Ordering` returns information about how `a` and `b` compare wtih respect to \\(=\\) and \\(\le\\).
