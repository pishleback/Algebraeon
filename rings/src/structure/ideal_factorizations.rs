use super::*;
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::*;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct DedekindDomainPrimeIdeal<Ideal> {
    prime_ideal: Ideal,
}

impl<Ideal> DedekindDomainPrimeIdeal<Ideal> {
    pub fn from_ideal_unchecked(ideal: Ideal) -> Self {
        Self { prime_ideal: ideal }
    }

    pub fn into_ideal(self) -> Ideal {
        self.prime_ideal
    }

    pub fn ideal(&self) -> &Ideal {
        &self.prime_ideal
    }
}

#[derive(Debug, Clone)]
pub struct DedekindDomainIdealFactorization<Ideal> {
    // The prime ideals should be distinct
    // All powers should be non-zero
    factors: Vec<(DedekindDomainPrimeIdeal<Ideal>, Natural)>,
}

impl<Ideal> DedekindDomainIdealFactorization<Ideal> {
    pub fn from_factor_powers(factors: Vec<(DedekindDomainPrimeIdeal<Ideal>, Natural)>) -> Self {
        Self { factors }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DedekindDomainIdealFactorizationStructure<
    RS: DedekindDomainSignature,
    RSB: BorrowedStructure<RS>,
    Ideals: DedekindDomainIdealsSignature<RS, RSB>,
    IdealsB: BorrowedStructure<Ideals>,
> {
    _ring: PhantomData<RS>,
    _ring_borrowed: PhantomData<RSB>,
    _ideals: PhantomData<Ideals>,
    ideals: IdealsB,
}

impl<
    RS: DedekindDomainSignature,
    RSB: BorrowedStructure<RS>,
    Ideals: DedekindDomainIdealsSignature<RS, RSB>,
    IdealsB: BorrowedStructure<Ideals>,
> DedekindDomainIdealFactorizationStructure<RS, RSB, Ideals, IdealsB>
{
    pub fn new(ideals: IdealsB) -> Self {
        Self {
            _ring: PhantomData::default(),
            _ring_borrowed: PhantomData::default(),
            _ideals: PhantomData::default(),
            ideals,
        }
    }

    pub fn ring(&self) -> &RS {
        self.ideals().ring()
    }

    pub fn ideals(&self) -> &Ideals {
        self.ideals.borrow()
    }
}

impl<
    RS: DedekindDomainSignature,
    RSB: BorrowedStructure<RS>,
    Ideals: DedekindDomainIdealsSignature<RS, RSB>,
    IdealsB: BorrowedStructure<Ideals>,
> Signature for DedekindDomainIdealFactorizationStructure<RS, RSB, Ideals, IdealsB>
{
}

impl<
    RS: DedekindDomainSignature,
    RSB: BorrowedStructure<RS>,
    Ideals: DedekindDomainIdealsSignature<RS, RSB>,
    IdealsB: BorrowedStructure<Ideals>,
> SetSignature for DedekindDomainIdealFactorizationStructure<RS, RSB, Ideals, IdealsB>
{
    type Set = DedekindDomainIdealFactorization<Ideals::Set>;

    fn is_element(&self, x: &Self::Set) -> bool {
        for (prime, power) in self.to_powers_unchecked(x) {
            if power == &Natural::ZERO {
                return false;
            }
            if self.try_object_is_prime(prime) == Some(false) {
                return false;
            }
        }
        true
    }
}

impl<
    Ring: DedekindDomainSignature,
    RingB: BorrowedStructure<Ring>,
    Ideals: DedekindDomainIdealsSignature<Ring, RingB>,
    IdealsB: BorrowedStructure<Ideals>,
> FactoredSignature for DedekindDomainIdealFactorizationStructure<Ring, RingB, Ideals, IdealsB>
{
    type PrimeObject = DedekindDomainPrimeIdeal<Ideals::Set>;
    type Object = Ideals::Set;

    fn object_divides(&self, a: &Self::Object, b: &Self::Object) -> bool {
        self.ideals().ideal_contains(a, b)
    }

    fn try_object_is_prime(&self, _object: &Self::PrimeObject) -> Option<bool> {
        None
    }

    fn prime_into_object(&self, prime: Self::PrimeObject) -> Self::Object {
        prime.into_ideal()
    }

    fn object_product(&self, objects: Vec<&Self::Object>) -> Self::Object {
        self.ideals()
            .ideal_product(objects.into_iter().cloned().collect())
    }

    fn new_powers_unchecked(&self, factor_powers: Vec<(Self::PrimeObject, Natural)>) -> Self::Set {
        Self::Set {
            factors: factor_powers,
        }
    }

    fn to_powers_unchecked<'a>(
        &self,
        a: &'a Self::Set,
    ) -> Vec<(&'a Self::PrimeObject, &'a Natural)> {
        a.factors.iter().map(|(p, k)| (p, k)).collect()
    }

    fn into_powers_unchecked(&self, a: Self::Set) -> Vec<(Self::PrimeObject, Natural)> {
        a.factors
    }

    fn expanded(&self, a: &Self::Set) -> Self::Object {
        self.ideals().ideal_product(
            self.to_primes(a)
                .into_iter()
                .map(|p| p.ideal().clone())
                .collect(),
        )
    }

    fn mul(&self, mut a: Self::Set, b: Self::Set) -> Self::Set {
        for (q, l) in self.into_powers(b) {
            'SEARCH_A: {
                for (p, k) in &mut a.factors {
                    if self.object_equivalent(p.ideal(), q.ideal()) {
                        *k += l;
                        break 'SEARCH_A;
                    }
                }
                a.factors.push((q, l));
            }
        }
        a
    }
}

// #[derive(Debug, Clone)]
// pub struct DedekindDomainFractionalIdeal<Ideal> {
//     // The prime ideals should be distinct
//     // All powers should be non-zero
//     factors: Vec<(DedekindDomainPrimeIdeal<Ideal>, Integer)>,
// }

// impl<Ideal> DedekindDomainFractionalIdeal<Ideal> {
//     pub fn from_factor_powers(factors: Vec<(DedekindDomainPrimeIdeal<Ideal>, Integer)>) -> Self {
//         Self { factors }
//     }
// }

// #[derive(Debug, Clone, PartialEq, Eq)]
// pub struct DedekindDomainFractionalIdealStructure<
//     RS: DedekindDomainSignature,
//     RSB: BorrowedStructure<RS>,
//     Ideals: DedekindDomainIdealsSignature<RS, RSB>,
//     IdealsB: BorrowedStructure<Ideals>,
// > {
//     _ring: PhantomData<RS>,
//     _ring_borrowed: PhantomData<RSB>,
//     _ideals: PhantomData<Ideals>,
//     ideals: IdealsB,
// }

// impl<
//     RS: DedekindDomainSignature,
//     RSB: BorrowedStructure<RS>,
//     Ideals: DedekindDomainIdealsSignature<RS, RSB>,
//     IdealsB: BorrowedStructure<Ideals>,
// > DedekindDomainFractionalIdealStructure<RS, RSB, Ideals, IdealsB>
// {
//     pub fn new(ideals: IdealsB) -> Self {
//         Self {
//             _ring: PhantomData::default(),
//             _ring_borrowed: PhantomData::default(),
//             _ideals: PhantomData::default(),
//             ideals,
//         }
//     }

//     pub fn ring(&self) -> &RS {
//         self.ideals().ring()
//     }

//     pub fn ideals(&self) -> &Ideals {
//         self.ideals.borrow()
//     }
// }

// impl<
//     RS: DedekindDomainSignature,
//     RSB: BorrowedStructure<RS>,
//     Ideals: DedekindDomainIdealsSignature<RS, RSB>,
//     IdealsB: BorrowedStructure<Ideals>,
// > Signature for DedekindDomainFractionalIdealStructure<RS, RSB, Ideals, IdealsB>
// {
// }

// impl<
//     RS: DedekindDomainSignature,
//     RSB: BorrowedStructure<RS>,
//     Ideals: DedekindDomainIdealsSignature<RS, RSB>,
//     IdealsB: BorrowedStructure<Ideals>,
// > SetSignature for DedekindDomainFractionalIdealStructure<RS, RSB, Ideals, IdealsB>
// {
//     type Set = DedekindDomainFractionalIdeal<Ideals::Set>;

//     fn is_element(&self, x: &Self::Set) -> bool {
//         for (prime, power) in &x.factors {
//             if power == &Integer::ZERO {
//                 return false;
//             }
//             if self.try_object_is_prime(prime) == Some(false) {
//                 return false;
//             }
//         }
//         true
//     }
// }
