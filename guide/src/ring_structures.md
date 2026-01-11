# Ring Structures

## Additive Structure

 - `Zero : Set` for a set with a special element called zero which we'll denote by `0`.
 - `Addition : Set` for a set with an associative and commutative binary operation called addition which we'll denote by `+`.
 - `AdditiveMonoid : Zero + Addition` when `a + 0 = 0 + a = a` for all `a`.
 - `CancellativeAddition : Addition` when `a + b = a + c` implies `b = c` for all `a`, `b` and `c`. In that case, the solution to `x + a = b` for `x` given `a` and `b` is unique when it exists.
 - `TryNegate : AdditiveMonoid` when solving `x + a = 0` for `x` given `a` has a unique solution whenever it exists.
 - `AdditiveGroup : CancellativeAddition + TryNegate` when finding additive inverses is always possible, so we have an additive group.

## Noncommutative Multiplicative Structure

 - `One : Set` for a set with a special element called one which we'll denote by `1`.
 - `Multiplication : Set` for a set with an associative binary operation called multiplication which we'll denote by `*`.
 - `MultiplicativeMonoid : One + Multiplication` when `a * 1 = 1 * a = a` for all `a`.
 - `LeftCancellativeMultiplication : Multiplication` when `x * a = x * b` implies `a = b` for all `a`, `b` and `x`. In that case, the solution to `x * a = b` for `x` given `a` and `b` is unique when it exists.
 - `RightCancellativeMultiplication : Multiplication` when `a * x = b * x` implies `a = b` for all `a`, `b` and `x`. In that case, the solution to `a * x = b` for `x` given `a` and `b` is uniqe when it exists.
 - `TryLeftReciprocal : MultiplicativeMonoid` when the solution to `x * a = 1` for `x` given `a` is unique whenever it exists.
 - `TryRightReciprocal : MultiplicativeMonoid` when the solution to `a * x = 1` for `x` given `a` is unique whenever it exists.
 - `TryTwoSidedReciprocal : MultiplicativeMonoid` when the solution to `x * a = a * x = 1` for `x` given `a` is unique whenever it exists.

## Combined Additive and Multiplicative Structure

 - `MultiplicativeAbsorptionMonoid : MultiplicativeMonoid + Zero` when we have `a * 0 = 0 * a = 0` for all `a`.
 - `MultiplicativeIntegralMonoid : MultiplicativeAbsorptionMonoid` when `a * b = 0` implies `a = 0` or `b = 0` and the submonoid of all elements not equal to `0` is cancellative.
 - `LeftDistributiveMultiplicationOverAddition : Addition + Multiplication` when we have `a * (b + c) = a * b + a * c` for all `a`, `b` and `c`.
 - `RightDistributiveMultiplicationOverAddition : Addition + Multiplication` when we have `(a + b) * c = a * c + b * c` for all `a`, `b` and `c`.
 - `DistributiveMultiplicationOverAddition : LeftDistributiveMultiplicationOverAddition + RightDistributiveMultiplicationOverAddition`.

## Commutative Semirings and Rings

 - `CommutativeMultiplication : Multiplication` when `*` is commutative, so `a * b = b * a` for all `a` and `b`.
 - `Semiring := AdditiveMonoid + MultiplicativeMonoid + MultiplicativeAbsorptionMonoid + DistributiveMultiplicationOverAddition`.
 - `Ring := AdditiveGroup + MultiplicativeMonoid + MultiplicativeAbsorptionMonoid + DistributiveMultiplicationOverAddition`.
 - `IntegralDomain := Ring + MultiplicativeIntegralMonoid`.
