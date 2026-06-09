use crate::num_theory::natural_factorization::factor_nat;
use crate::num_theory::natural_factorization::primes::is_prime_nat;
use crate::structure::{
    Factored, FactoringMonoidSignature, FactoringStructure, MetaFactoringMonoid,
    MultiplicationSignature, MultiplicativeMonoidSignature, SemiRingSignature,
    UniqueFactorizationMonoidSignature,
};
use algebraeon_nzq::traits::ModPow;
use algebraeon_nzq::{Natural, NaturalCanonicalStructure, gcd};
use algebraeon_structures::*;
use itertools::Itertools;

impl UniqueFactorizationMonoidSignature for NaturalCanonicalStructure {
    type FactoredExponent = NaturalCanonicalStructure;

    fn factorization_exponents(&self) -> &Self::FactoredExponent {
        Natural::structure_ref()
    }

    fn into_factorization_exponents(self) -> Self::FactoredExponent {
        Natural::structure()
    }

    fn try_is_irreducible(&self, a: &Self::Elem) -> Option<bool> {
        Some(is_prime_nat(a))
    }

    fn factorization_pow(&self, a: &Self::Elem, k: &Natural) -> Self::Elem {
        self.nat_pow(a, k)
    }
}

impl FactoringMonoidSignature for NaturalCanonicalStructure {
    fn is_irreducible(&self, a: &Self::Elem) -> bool {
        is_prime_nat(a)
    }

    fn factor_unchecked(
        &self,
        a: &Self::Elem,
    ) -> Factored<Self::Elem, <Self::FactoredExponent as SetSignature>::Elem> {
        factor_nat(a.clone())
    }
}

impl<
    ObjectB: BorrowedStructure<NaturalCanonicalStructure>,
    ExponentB: BorrowedStructure<NaturalCanonicalStructure>,
> FactoringStructure<NaturalCanonicalStructure, ObjectB, NaturalCanonicalStructure, ExponentB>
{
    pub fn euler_totient(&self, a: &Factored<Natural, Natural>) -> Natural {
        #[cfg(debug_assertions)]
        self.validate_element(a).unwrap();
        match a {
            Factored::Zero => {
                // The number of units in the quotient ring Z/0Z = Z is 2 i.e. +1 and -1
                self.objects().from_nat(&Natural::TWO)
            }
            Factored::NonZero(a) => {
                let mut prod = a.unit().clone();
                for (p, k) in a.powers() {
                    debug_assert!(p != &Natural::ZERO);
                    debug_assert!(k != &Natural::ZERO);
                    self.objects().mul_mut(
                        &mut prod,
                        &((p - &Natural::ONE) * p.pow(&(k - &Natural::ONE))),
                    );
                }
                prod
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IsPrimitiveRootResult {
    NonUnit,
    No,
    Yes,
}

impl<
    PowersB: BorrowedStructure<NaturalCanonicalStructure>,
    ExponentB: BorrowedStructure<NaturalCanonicalStructure>,
> FactoringStructure<NaturalCanonicalStructure, PowersB, NaturalCanonicalStructure, ExponentB>
{
    /// Return whether x is a primitive root modulo the factorized value
    pub fn is_primitive_root(
        &self,
        x: &Natural,
        n_factored: &Factored<Natural, Natural>,
    ) -> IsPrimitiveRootResult {
        #[cfg(debug_assertions)]
        self.validate_element(n_factored).unwrap();

        let factorizations = Natural::structure_ref().factorizations();
        let n = factorizations.expand(n_factored);
        if gcd(x.clone(), n.clone()) != Natural::ONE {
            IsPrimitiveRootResult::NonUnit
        } else {
            let phi_n = factorizations.euler_totient(n_factored);
            let x_mod_n = x % &n;
            for p in phi_n.clone().factor().distinct_irreducibles().unwrap() {
                if (&x_mod_n).mod_pow(&phi_n / p, &n) == Natural::ONE {
                    return IsPrimitiveRootResult::No;
                }
            }
            IsPrimitiveRootResult::Yes
        }
    }
}

impl<
    PowersB: BorrowedStructure<NaturalCanonicalStructure>,
    ExponentB: BorrowedStructure<NaturalCanonicalStructure>,
> ToStringSignature
    for FactoringStructure<NaturalCanonicalStructure, PowersB, NaturalCanonicalStructure, ExponentB>
{
    fn to_string(&self, elem: &Self::Elem) -> String {
        use std::fmt::Write;
        let mut f = String::new();
        if let Some(powers) = elem.powers() {
            if powers.is_empty() {
                write!(f, "1").unwrap();
            } else {
                for (i, (p, k)) in powers
                    .iter()
                    .sorted_by_cached_key(|(p, _k)| (*p).clone())
                    .enumerate()
                {
                    if i != 0 {
                        write!(f, " × ").unwrap();
                    }
                    write!(f, "{}", p).unwrap();
                    if k != &Natural::ONE {
                        write!(f, "^").unwrap();
                        write!(f, "{}", k).unwrap();
                    }
                }
            }
        } else {
            write!(f, "0").unwrap();
        }
        f
    }
}
