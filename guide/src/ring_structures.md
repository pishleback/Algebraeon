# Ring Structures

Let \\(X\\) be a set. This section outlines the traits Algebraeon provides for defining group-like structures on \\(X\\).

## Additive Structure

 - `Zero : Set` for a set with a special element \\(0 \in X\\). The method `zero() -> X` returns an instance of `X` representing \\(0\\).
 - `Addition : Set` for a set with an associative and commutative binary operation \\(+\\).
\\[+ : X \times X \to X : (a, b) \mapsto a + b\\]
\\[a + (b + c) = (a + b) + c \quad \forall a, b, c \in X\\]
\\[a + b = b + a \quad \forall a, b \in X\\]
The method `add(a: X, b: X) -> X` computes \\(a + b\\).
 - `AdditiveMonoid : Zero + Addition` when
\\[a + 0 = 0 + a = a \quad \forall a \in X\\]
 - `CancellativeAddition : Composition` when 
\\[a + x = a + y \implies x = y \quad \forall a, x, y \in X\\]
In that case the solution (or lack thereof) to \\(a = b + x\\) for \\(x\\) given \\(a\\) and \\(b\\) is unique whenever it exists and is computed using `.try_sub(a: X, b: X) -> Option<X>`.
 - `TryNegate : AdditiveMonoid` when the solution to \\[x + a = 0\\] for \\(x\\) given \\(a\\) is unique whenever it exists. The solution (or lack thereof) is computed using `.try_neg(a: X) -> Option<X>`.
 - `AdditiveGroup : TryNegate + CancellativeAddition` when every element has an additive inverse which can be computed using `.neg(a: X) -> X`.

## Multiplicative Structure

 - `One : Set` for a set with a special element \\(1 \in X\\). The method `one() -> X` returns an instance of `X` representing \\(1\\).
 - `Multiplication : Set` for a set with an associative binary operation \\(\times\\).
\\[\times : X \times X \to X : (a, b) \mapsto a \times b\\]
\\[a \times (b \times c) = (a \times b) \times c \quad \forall a, b, c \in X\\]
The method `mul(a: X, b: X) -> X` computes \\(a \times b\\).
 - `CommutativeMultiplication : Multiplication` when `*` is commutative, so `a * b = b * a` for all `a` and `b`.
 - `MultiplicativeMonoid : One + Multiplication` when
\\[a \times 1 = 1 \times a = a \quad \forall a \in X\\]

## Combined Additive and Multiplicative Structure

 - `MultiplicativeAbsorptionMonoid : MultiplicativeMonoid + Zero` when
 \\[a \times 0 = 0 \times a = 0 \quad \forall a \in X\\]
 - `LeftDistributiveMultiplicationOverAddition : Addition + Multiplication` when
 \\[a \times (b + c) = (a \times b) + (a \times c) \quad \forall a, b, c \in X\\]
 - `RightDistributiveMultiplicationOverAddition : Addition + Multiplication` when
 \\[(a + b) \times c = (a \times c) + (b \times c) \quad \forall a, b, c \in X\\]
 - `DistributiveMultiplicationOverAddition := LeftDistributiveMultiplicationOverAddition + RightDistributiveMultiplicationOverAddition`.
 - `LeftCancellativeMultiplication : Multiplication + Zero` when 
 \\[a \times x = a \times y \implies x = y \quad \forall x, y \in X \quad a \in X \setminus \\{0\\}\\]
 In that case the solution (or lack thereof) to \\(a = b \times x\\) for \\(x\\) given \\(a\\) and non-zero \\(b\\) is unique whenever it exists and is computed using `.try_left_div(a: X, b: X) -> Option<X>`. `None` is also returned if \\(b = 0\\).
 - `RightCancellativeMultiplication : Multiplication + Zero` when 
 \\[x \times a = y \times a \implies x = y \quad \forall x, y \in X \quad a \in X \setminus \\{0\\}\\]
 In that case the solution (or lack thereof) to \\(a = x \times b\\) for \\(x\\) given \\(a\\) and non-zero \\(b\\) is unique whenever it exists and is computed using `.try_right_div(a: X, b: X) -> Option<X>`. `None` is also returned if \\(b = 0\\).
 - `CancellativeMultiplication := CommutativeMultiplication + LeftCancellativeMultiplication + RightCancellativeMultiplication`. `.try_div(a: X, b: X) -> Option<X>` can be used as an alias for `try_left_div` and (equivalently) `try_right_div`.
 - `TryLeftReciprocal : MultiplicativeMonoid` when the solution to \\[x \times a = 1\\] for \\(x\\) given \\(a\\) is unique whenever it exists. The solution (or lack thereof) is computed using `.try_left_reciprocal(a: X) -> Option<X>`.
 - `TryRightReciprocal : MultiplicativeMonoid` when the solution to \\[a \times x = 1\\] for \\(x\\) given \\(a\\) is unique whenever it exists. The solution (or lack thereof) is computed using `.try_right_reciprocal(a: X) -> Option<X>`.
 - `TryReciprocal : MultiplicativeMonoid` when the solution to \\[x \times a = a \times x = 1\\] for \\(x\\) given \\(a\\) is unique whenever it exists. The solution (or lack thereof) is computed using `.try_reciprocal(a: X) -> Option<X>`.
 - `MultiplicativeIntegralMonoid : MultiplicativeAbsorptionMonoid + TryReciprocal + TryLeftReciprocal + TryRightReciprocal + LeftCancellativeMultiplication + RightCancellativeMultiplication` when 
 \\[a \times b = 0 \implies a = 0 \quad \text{or} \quad b = 0 \quad \forall a, b \in X\\]

## Semirings and Rings

 - `Semiring := AdditiveMonoid + CommutativeMultiplication + MultiplicativeAbsorptionMonoid + DistributiveMultiplicationOverAddition`.
 - `IntegralSemiring := Semiring + MultiplicativeIntegralMonoid`.
 - `Ring := AdditiveGroup + Semiring`.
 - `IntegralDomain := Ring + MultiplicativeIntegralMonoid`.
 - `EuclideanDivision : Semiring` when there exists a norm function
 \\[N : X \to \mathbb{N}\\]
 computed using `norm(a: X) -> Natural` such that for all \\(a \in X\\) and all non-zero \\(b \in X\\) there exists \\(q, r \in X\\) computed using `quorem(a: X, b: X) -> (X, X)` such that
 \\[a = qb + r \quad \text{and either} \quad N(r) < N(b) \quad \text{or} \quad r = 0\\]
 - `EuclideanDomain := IntegralDomain + EuclideanDivision`.