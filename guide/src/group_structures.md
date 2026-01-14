# Group Structures

Let \\(X\\) be a set. This section outlines the traits Algebraeon provides for defining group-like structures on \\(X\\).

## Composition

 - `Composition : Set` for a set with a binary operation called composition denoted by \\(\circ\\).
\\[\circ : X \times X \to X : (a, b) \mapsto a \circ b\\]
`.compose(a: X, b: X) -> X` is used to compute \\(a \circ b\\).
 - `AssociativeComposition : Composition` when \\(\circ\\) is associative
\\[a \circ (b \circ c) = (a \circ b) \circ c \quad \forall a, b, c \in X\\]
 - `CommutativeComposition : Composition` when \\(\circ\\) is commutative
 \\[a \circ b = b \circ a \quad \forall a, b \in X\\]
 - `LeftCancellativeComposition : Composition` when 
\\[a \circ x = a \circ y \implies x = y \quad \forall a, x, y \in X\\]
In that case the solution (or lack thereof) to \\(a = b \circ x\\) for \\(x\\) given \\(a\\) and \\(b\\) is unique whenever it exists and is computed using `.try_left_difference(a: X, b: X) -> Option<X>`.
 - `RightCancellativeComposition : Composition` when 
\\[x \circ a = y \circ a \implies x = y \quad \forall a, x, y \in X\\]
In that case the solution (or lack thereof) to \\(a = x \circ b\\) for \\(x\\) given \\(a\\) and \\(b\\) is unique whenever it exists and is computed using `.try_right_difference(a: X, b: X) -> Option<X>`.
 - `CancellativeComposition := CommutativeComposition + LeftCancellativeComposition + RightCancellativeComposition`. `.try_difference(a: X, b: X) -> Option<X>` can be used as an alias for `try_left_difference` and (equivalently) `try_right_difference`.

## Identity

 - `Identity : Set` for a set with a special element called the identity element which we'll denote by `e`.
 - `TryLeftInverse : Identity + Composition` when the solution to \\[x \circ a = e\\] for \\(x\\) given \\(a\\) is unique whenever it exists. The solution (or lack thereof) is computed using `.try_left_inverse(a: X) -> Option<X>`.
 - `TryRightInverse : Identity + Composition` when the solution to \\[a \circ x = e\\] for \\(x\\) given \\(a\\) is unique whenever it exists. The solution (or lack thereof) is computed using `.try_right_inverse(a: X) -> Option<X>`.
 - `TryInverse: Identity + Composition` when the solution to \\[x \circ a = a \circ x = e\\] for \\(x\\) given \\(a\\) is unique whenever it exists. The solution (or lack thereof) is computed using `.try_inverse(a: X) -> Option<X>`.
 - `Monoid : Identity + AssociativeComposition + TryInverse` when 
\\[a \circ e = e \circ a = a \quad \forall a \in X\\]
 - `Group : TryInverse + TryLeftInverse + TryRightInverse + LeftCancellativeComposition + RightCancellativeComposition` when every element has an inverse. Left-, right-, and two-sided-inverses all coencide in this case and are computed using `.inverse(a: X) -> X`.
 - `AbelianGroup := Group + CommutativeComposition`.

