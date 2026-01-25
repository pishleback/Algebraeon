# Structure Types

## Motivation

In mathematics there are many instances of sets with some additional structure, for example:
 - The set of integers with its ordered ring structure.
 - The set of rational numbers with its ordered field structure.
 - For each natural number \\(n \in \mathbb{N}\\), the finite set of integers modulo \\(n\\) with its ring structure.
 - The set of all algebraic numbers in \\(\mathbb{C}\\) with its field structure.
 - The set of all ideals in a ring with the operations of ideal addition, ideal intersection, and ideal multiplication. Depending on the ring, ideals may be uniquely factorable as a product of prime ideals.

The approach taken by Algebraeon to represent such sets with additional structure is well illustrated by the example of the ring of integers modulo \\(n \in \mathbb{N}\\). This would be done as follows:
 - Define a structure type `IntegersModuloN` whose objects shall represent the ring of integers modulo \\(n\\) for different values of \\(n\\).
 - Implement the desired structure by implementing signature traits on the structure type. In the case of `IntegersModuloN` the required signature traits are:
   - `Signature` the base signature trait.
   - `SetSignature<Set = Integer>` so that instances of `Integer` shall be used to represent elements of the set of integers modulo \\(n\\).
   - `EqSignature` so that pairs of `Integer`s can be tested for equality modulo \\(n\\).
   - `FiniteSetSignature` so that a list of all integers modulo \\(n\\) can be produced.
   - `SemiRingSignature` so that `Integer`s can be added and multiplied modulo \\(n\\).
   - `RingSignature` so that `Integer`s can be subtracted modulo \\(n\\).
  
In practice this looks like
```rust
use algebraeon::nzq::{Integer, Natural};
use algebraeon::rings::structure::*;
use algebraeon::sets::structure::*;

#[derive(Debug, Clone, PartialEq, Eq)]
struct IntegersModuloN {
    n: Natural,
}

impl Signature for IntegersModuloN {}

impl SetSignature for IntegersModuloN {
    type Set = Integer;

    fn validate_element(&self, x: &Integer) -> Result<(), String> {
        if x >= &self.n {
            return Err("too big".to_string());
        }
        Ok(())
    }
}

impl EqSignature for IntegersModuloN {
    fn equal(&self, a: &Integer, b: &Integer) -> bool {
        a == b
    }
}

impl RinglikeSpecializationSignature for IntegersModuloN {}

impl ZeroSignature for IntegersModuloN {
    fn zero(&self) -> Self::Set {
        Integer::ZERO
    }
}

impl AdditionSignature for IntegersModuloN {
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        ((a + b) % &self.n).into()
    }
}

impl CancellativeAdditionSignature for IntegersModuloN {
    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.sub(a, b))
    }
}

impl TryNegateSignature for IntegersModuloN {
    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.neg(a))
    }
}

impl AdditiveMonoidSignature for IntegersModuloN {}

impl AdditiveGroupSignature for IntegersModuloN {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        (-a % &self.n).into()
    }
}

impl OneSignature for IntegersModuloN {
    fn one(&self) -> Self::Set {
        (Integer::ONE % &self.n).into()
    }
}

impl MultiplicationSignature for IntegersModuloN {
    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        ((a * b) % &self.n).into()
    }
}

impl CommutativeMultiplicationSignature for IntegersModuloN {}

impl LeftDistributiveMultiplicationOverAddition for IntegersModuloN {}

impl RightDistributiveMultiplicationOverAddition for IntegersModuloN {}

impl MultiplicativeMonoidSignature for IntegersModuloN {}

impl MultiplicativeAbsorptionMonoidSignature for IntegersModuloN {}

impl SemiRingSignature for IntegersModuloN {}

impl RingSignature for IntegersModuloN {}

let mod_6 = IntegersModuloN { n: 6u32.into() };
// Since we've given `mod_6` the structure of a ring, Algebraeon implements
// the repeated squaring algorithm for taking very large powers modulo `n`.
assert!(mod_6.equal(
    &mod_6.nat_pow(&2.into(), &1000000000000u64.into()),
    &4.into()
));
```