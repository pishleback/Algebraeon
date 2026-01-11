# Group Structures

Let \\(X\\) be a set. Algebraeon has the following traits for defining group-like structures on \\(X\\).

## Composition

 - `Composition : Set` for a set with a binary operation called composition denoted by \\(\circ\\).
\\[\circ : X \times X \to X : (a, b) \mapsto a \circ b\\]
`.compose(a: X, b: X) -> X` is used to compute \\(a \circ b\\).
 - `AssociativeComposition : Composition` when \\(\circ\\) is associative
\\[a \circ (b \circ c) = (a \circ b) \circ c \quad \forall a, b, c \in X\\]
 - `LeftCancellativeComposition : Composition` when 
\\[x \circ a = x \circ b \implies a = b \quad \forall x, a, b \in X\\]
In that case the solution to \\(x \circ a = b\\) for \\(x\\) given \\(a\\) and \\(b\\) is unique whenever it exists and is computed using `.try_left_difference(a: X, b: X) -> Option<X>`.
 - `RightCancellativeComposition : Composition` when
\\[a \circ x = b \circ x \implies a = b \quad \forall x, a, b \in X\\]
In that case the solution to \\(a \circ x = b\\) for \\(x\\) given \\(a\\) and \\(b\\) is unique whenever it exists and is computed using `.try_right_difference(a: X, b: X) -> Option<X>`.
 
## Identity

 - `Identity : Set` for a set with a special element called the identity element which we'll denote by `1`.
 - `Monoid : Identity + AssociativeComposition` when `a * 1 = 1 * a = a` for all `a`.
 - `TryLeftInverse : Monoid` when the solution to `x * a = 1` for `x` given `a` is unique whenever it exists.
 - `TryRightInverse: Monoid` when the solution to `a * x = 1` for `x` given `a` is unique whenever it exists.
 - `TryTwoSidedInverse: Monoid` when the solution to `x * a = a * x = 1` for `x` given `a` is unique whenever it exists.
 - `Group : TryTwoSidedInverse + TryLeftInverse + TryRightInverse + LeftCancellativeComposition + RightCancellativeComposition + TryTwoSidedInverse` when inverses always exist.

## Commutative Composition

 - `CommutativeComposition : Composition` when `*` is commutative, so `a * b = b * a` for all `a` and `b`.
 - `AbelianGroup := Group + CommutativeComposition`.
