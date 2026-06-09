use crate::{
    num_theory::natural_factorization::primes::is_prime_nat,
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, CancellativeMultiplicationSignature,
        CommutativeMultiplicationSignature, FieldSignature, FiniteFieldSignature,
        IntegralDomainSignature, LeftDistributiveMultiplicationOverAddition,
        MultiplicationSignature, MultiplicativeAbsorptionMonoidSignature,
        MultiplicativeIntegralMonoidSignature, MultiplicativeMonoidSignature,
        MultiplicativeMonoidUnitsStructure, OneSignature, QuotientRingGetPrincipalIdealSignature,
        QuotientRingSignature, RightDistributiveMultiplicationOverAddition, RingSignature,
        RinglikeSpecializationSignature, SemiRingSignature, TryNegateSignature,
        TryReciprocalSignature, ZeroSignature,
    },
};
use algebraeon_structures::*;
use rand::{RngExt, SeedableRng, rngs::StdRng};
use std::{borrow::Cow, ops::Rem};

fn xgcd_u64(mut x: u64, mut y: u64) -> (u64, i64, i64) {
    let mut pa = 1;
    let mut a = 0;
    let mut pb = 0;
    let mut b = 1;
    while y != 0 {
        let (q, r) = (x / y, x % y);
        let new_a = pa - q as i64 * a;
        (a, pa) = (new_a, a);
        let new_b = pb - q as i64 * b;
        (b, pb) = (new_b, b);
        (x, y) = (y, r);
    }
    (x, pa, pb)
}

// Return the inverse of x modulo n if it exists
fn inv_mod_n_u64(x: u64, n: u64) -> Option<u64> {
    let (g, a, _) = xgcd_u64(x, n);
    if g == 1 {
        Some(a.rem_euclid(n as i64).try_into().unwrap())
    } else {
        None
    }
}

/// The ring of integers modulo odd n
/// where elements x are represented in Montgomery form as a 0 <= (x*r % n) < n
/// where r is the smallest power of 2 such that r > n.
#[derive(Debug, Clone)]
pub struct MontgomeryModuloOddStructure {
    n: u64,
    r: u64,
    r_mod_n: u64,
    mod_r_mask: u64, // r-1
    r_pow: u32,
    n_inv_neg_mod_r: u64,
    r_squared_mod_n: u64,
    r_cubed_mod_n: u64,
}

impl MontgomeryModuloOddStructure {
    pub const fn max_modulus() -> u64 {
        1u64 << 31
    }

    pub fn new_unchecked(n: u64) -> Self {
        debug_assert_eq!(n % 2, 1);
        // Need a bound such that all operations remain in u64
        debug_assert!(n < Self::max_modulus());
        let r = n.next_power_of_two();
        let r_mod_n = r % n;
        let mod_r_mask = r - 1;
        debug_assert_eq!(mod_r_mask + 1, r);
        let r_pow = r.ilog2();
        debug_assert_eq!(1 << r_pow, r);
        debug_assert_eq!(r.count_ones(), 1);
        debug_assert!(n < r);
        let n_inv_neg_mod_r = r - inv_mod_n_u64(n, r).unwrap();
        debug_assert_eq!((n * n_inv_neg_mod_r) & (r - 1), r - 1); // n * n_inv_neg_mod_r == -1 mod r
        let r_squared_mod_n = (r_mod_n * r_mod_n) % n;
        debug_assert_eq!((r * r) % n, r_squared_mod_n);
        let r_cubed_mod_n = (r_squared_mod_n * r_mod_n) % n;
        debug_assert_eq!(((r * r) % n * r) % n, r_cubed_mod_n);
        Self {
            n,
            r,
            r_mod_n,
            mod_r_mask,
            r_pow,
            n_inv_neg_mod_r,
            r_squared_mod_n,
            r_cubed_mod_n,
        }
    }

    pub fn n(&self) -> u64 {
        self.n
    }

    pub fn r(&self) -> u64 {
        self.r
    }

    /// Input:
    ///     t in the range [0, n*r)
    /// Output:
    ///     s in the range [0, n) such that s == t/r mod n
    pub fn montgomery_reduction(&self, t: u64) -> u64 {
        // https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
        // The efficiency of this algorithm is the reason Montgomery form is useful.
        // `self.n_inv_neg_mod_r` = n' is such that n * n' == -1 mod r
        // ((t % r) * n') % r
        let m = ((t & self.mod_r_mask) * self.n_inv_neg_mod_r) & self.mod_r_mask;
        // (t + m * n) / r
        let s = (t + m * self.n) >> self.r_pow;
        debug_assert!(s < 2 * self.n);
        // s mod n
        let s = if s >= self.n { s - self.n } else { s };
        debug_assert!(s < self.n);
        debug_assert_eq!((s * self.r) % self.n, t % self.n);
        s
    }

    pub fn to_montgomery_form(&self, x: u64) -> u64 {
        debug_assert!(x < self.n);
        self.montgomery_reduction(x * self.r_squared_mod_n)
    }

    pub fn from_montgomery_form(&self, x: u64) -> u64 {
        debug_assert!(x < self.n);
        self.montgomery_reduction(x)
    }
}

impl PartialEq for MontgomeryModuloOddStructure {
    fn eq(&self, other: &Self) -> bool {
        self.n == other.n
    }
}

impl Eq for MontgomeryModuloOddStructure {}

impl Signature for MontgomeryModuloOddStructure {}

impl SetSignature for MontgomeryModuloOddStructure {
    type Elem = u64; // 0 <= x < n

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if x >= &self.n {
            return Err(format!(
                "{x} is out of the valid range 0 <= {x} < {}",
                self.n
            ));
        }
        Ok(())
    }
}

impl ToStringSignature for MontgomeryModuloOddStructure {
    fn to_string(&self, elem: &Self::Elem) -> String {
        format!("{}", self.from_montgomery_form(*elem))
    }
}

impl EqSignature for MontgomeryModuloOddStructure {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        #[cfg(debug_assertions)]
        {
            self.validate_element(a).unwrap();
            self.validate_element(b).unwrap();
        }
        a == b
    }
}

impl RinglikeSpecializationSignature for MontgomeryModuloOddStructure {}

impl ZeroSignature for MontgomeryModuloOddStructure {
    fn zero(&self) -> Self::Elem {
        0
    }
}

impl AdditionSignature for MontgomeryModuloOddStructure {
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        #[cfg(debug_assertions)]
        {
            self.validate_element(a).unwrap();
            self.validate_element(b).unwrap();
        }
        let s = a + b;
        let s = if s >= self.n { s - self.n } else { s };
        #[cfg(debug_assertions)]
        self.validate_element(&s).unwrap();
        s
    }
}

impl TryNegateSignature for MontgomeryModuloOddStructure {
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl CancellativeAdditionSignature for MontgomeryModuloOddStructure {
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.add(a, &self.neg(b)))
    }
}

impl AdditiveMonoidSignature for MontgomeryModuloOddStructure {}

impl AdditiveGroupSignature for MontgomeryModuloOddStructure {
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        #[cfg(debug_assertions)]
        self.validate_element(a).unwrap();
        let b = if *a == 0 { *a } else { self.n - a };
        #[cfg(debug_assertions)]
        self.validate_element(&b).unwrap();
        b
    }
}

impl OneSignature for MontgomeryModuloOddStructure {
    fn one(&self) -> Self::Elem {
        self.r_mod_n
    }
}

impl MultiplicationSignature for MontgomeryModuloOddStructure {
    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.montgomery_reduction(a * b)
    }
}

impl MultiplicativeMonoidSignature for MontgomeryModuloOddStructure {}

impl CommutativeMultiplicationSignature for MontgomeryModuloOddStructure {}

impl RightDistributiveMultiplicationOverAddition for MontgomeryModuloOddStructure {}

impl LeftDistributiveMultiplicationOverAddition for MontgomeryModuloOddStructure {}

impl MultiplicativeAbsorptionMonoidSignature for MontgomeryModuloOddStructure {}

impl TryReciprocalSignature for MontgomeryModuloOddStructure {
    fn try_reciprocal(&self, a: &Self::Elem) -> Option<Self::Elem> {
        let b = self.montgomery_reduction(inv_mod_n_u64(*a, self.n)? * self.r_cubed_mod_n);
        debug_assert!(self.equal(&self.mul(a, &b), &self.one()));
        Some(b)
    }
}

impl SemiRingSignature for MontgomeryModuloOddStructure {}

impl RingSignature for MontgomeryModuloOddStructure {}

impl QuotientSetSignature<IntegerCanonicalStructure> for MontgomeryModuloOddStructure {
    fn pre_quotient_set(&self) -> &IntegerCanonicalStructure {
        Integer::structure_ref()
    }

    fn project(&self, x: Integer) -> Self::Elem {
        self.to_montgomery_form(x.rem(Natural::from(self.n)).try_into().unwrap())
    }

    fn project_ref(&self, x: &Integer) -> Self::Elem {
        self.to_montgomery_form(x.rem(Natural::from(self.n)).try_into().unwrap())
    }

    fn unproject(&self, x: Self::Elem) -> Integer {
        Integer::from(self.from_montgomery_form(x))
    }

    fn unproject_ref(&self, x: &Self::Elem) -> Integer {
        self.unproject(*x)
    }
}

impl QuotientRingSignature<IntegerCanonicalStructure> for MontgomeryModuloOddStructure {}

impl QuotientRingGetPrincipalIdealSignature<IntegerCanonicalStructure>
    for MontgomeryModuloOddStructure
{
    fn modulus<'a>(&'a self) -> Cow<'a, Integer> {
        Cow::Owned(Integer::from(self.n))
    }
}

impl CountableSetSignature for MontgomeryModuloOddStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> + Clone {
        0..self.n
    }
}

impl FiniteSetSignature for MontgomeryModuloOddStructure {
    fn generate_random_elements(&self, seed: u64) -> impl Iterator<Item = Self::Elem> {
        let mut rng = StdRng::seed_from_u64(seed);
        (0..).map(move |_| rng.random_range(0..self.n))
    }
}

#[derive(Debug, Clone)]
pub struct MontgomeryModuloOddPrimeStructure {
    sup: MontgomeryModuloOddStructure,
    inv_cache: Option<Vec<u64>>,
}

impl MontgomeryModuloOddPrimeStructure {
    pub const fn max_modulus() -> u64 {
        MontgomeryModuloOddStructure::max_modulus()
    }

    pub fn new_unchecked(n: u64) -> Self {
        let sup = MontgomeryModuloOddStructure::new_unchecked(n);
        #[cfg(debug_assertions)]
        debug_assert!(is_prime_nat(&Natural::from(n)));
        Self {
            sup,
            inv_cache: None,
        }
    }

    pub fn populate_inv_cache(&mut self) {
        if self.inv_cache.is_none() {
            self.inv_cache = Some(
                (0..self.sup.n)
                    .map(|a| self.try_reciprocal(&a).unwrap_or(0))
                    .collect::<Vec<_>>(),
            )
        }
    }

    pub fn montgomery_reduction(&self, t: u64) -> u64 {
        self.sup.montgomery_reduction(t)
    }

    pub fn to_montgomery_form(&self, x: u64) -> u64 {
        self.sup.to_montgomery_form(x)
    }

    pub fn from_montgomery_form(&self, x: u64) -> u64 {
        self.sup.from_montgomery_form(x)
    }
}

impl PartialEq for MontgomeryModuloOddPrimeStructure {
    fn eq(&self, other: &Self) -> bool {
        self.sup == other.sup
    }
}

impl Eq for MontgomeryModuloOddPrimeStructure {}

impl Signature for MontgomeryModuloOddPrimeStructure {}

impl SetSignature for MontgomeryModuloOddPrimeStructure {
    type Elem = u64; // 0 <= x < n

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        let n = self.sup.n();
        if !is_prime_nat(&Natural::from(n)) {
            return Err(format!("Modulus {} is not prime", n));
        }
        self.sup.validate_element(x)?;
        Ok(())
    }
}

impl ToStringSignature for MontgomeryModuloOddPrimeStructure {
    fn to_string(&self, elem: &Self::Elem) -> String {
        self.sup.to_string(elem)
    }
}

impl EqSignature for MontgomeryModuloOddPrimeStructure {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        self.sup.equal(a, b)
    }
}

impl RinglikeSpecializationSignature for MontgomeryModuloOddPrimeStructure {}

impl ZeroSignature for MontgomeryModuloOddPrimeStructure {
    fn zero(&self) -> Self::Elem {
        self.sup.zero()
    }
}

impl AdditionSignature for MontgomeryModuloOddPrimeStructure {
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.sup.add(a, b)
    }
}

impl TryNegateSignature for MontgomeryModuloOddPrimeStructure {
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        self.sup.try_neg(a)
    }
}

impl CancellativeAdditionSignature for MontgomeryModuloOddPrimeStructure {
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        self.sup.try_sub(a, b)
    }
}

impl AdditiveMonoidSignature for MontgomeryModuloOddPrimeStructure {}

impl AdditiveGroupSignature for MontgomeryModuloOddPrimeStructure {
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        self.sup.neg(a)
    }
}

impl OneSignature for MontgomeryModuloOddPrimeStructure {
    fn one(&self) -> Self::Elem {
        self.sup.one()
    }
}

impl MultiplicationSignature for MontgomeryModuloOddPrimeStructure {
    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.sup.mul(a, b)
    }
}

impl MultiplicativeMonoidSignature for MontgomeryModuloOddPrimeStructure {}

impl CommutativeMultiplicationSignature for MontgomeryModuloOddPrimeStructure {}

impl RightDistributiveMultiplicationOverAddition for MontgomeryModuloOddPrimeStructure {}

impl LeftDistributiveMultiplicationOverAddition for MontgomeryModuloOddPrimeStructure {}

impl MultiplicativeAbsorptionMonoidSignature for MontgomeryModuloOddPrimeStructure {}

impl SemiRingSignature for MontgomeryModuloOddPrimeStructure {}

impl RingSignature for MontgomeryModuloOddPrimeStructure {}

impl QuotientSetSignature<IntegerCanonicalStructure> for MontgomeryModuloOddPrimeStructure {
    fn pre_quotient_set(&self) -> &IntegerCanonicalStructure {
        self.sup.pre_quotient_set()
    }

    fn project(&self, x: Integer) -> Self::Elem {
        self.sup.project(x)
    }

    fn project_ref(&self, x: &Integer) -> Self::Elem {
        self.sup.project_ref(x)
    }

    fn unproject(&self, x: Self::Elem) -> Integer {
        self.sup.unproject(x)
    }

    fn unproject_ref(&self, x: &Self::Elem) -> Integer {
        self.sup.unproject_ref(x)
    }
}

impl QuotientRingSignature<IntegerCanonicalStructure> for MontgomeryModuloOddPrimeStructure {}

impl QuotientRingGetPrincipalIdealSignature<IntegerCanonicalStructure>
    for MontgomeryModuloOddPrimeStructure
{
    fn modulus<'a>(&'a self) -> Cow<'a, Integer> {
        self.sup.modulus()
    }
}

impl TryReciprocalSignature for MontgomeryModuloOddPrimeStructure {
    fn try_reciprocal(&self, a: &Self::Elem) -> Option<Self::Elem> {
        if let Some(cache) = &self.inv_cache {
            if *a == 0 {
                None
            } else {
                Some(cache[*a as usize])
            }
        } else {
            let b =
                self.montgomery_reduction(inv_mod_n_u64(*a, self.sup.n)? * self.sup.r_cubed_mod_n);
            debug_assert!(self.equal(&self.mul(a, &b), &self.one()));
            Some(b)
        }
    }
}

impl CancellativeMultiplicationSignature for MontgomeryModuloOddPrimeStructure {
    fn try_divide(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.mul(a, &self.try_reciprocal(b)?))
    }
}

impl MultiplicativeIntegralMonoidSignature for MontgomeryModuloOddPrimeStructure {}

impl IntegralDomainSignature for MontgomeryModuloOddPrimeStructure {}

impl FieldSignature for MontgomeryModuloOddPrimeStructure {}

impl CountableSetSignature for MontgomeryModuloOddPrimeStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> + Clone {
        self.sup.generate_all_elements()
    }
}

impl FiniteSetSignature for MontgomeryModuloOddPrimeStructure {
    fn generate_random_elements(&self, seed: u64) -> impl Iterator<Item = Self::Elem> {
        self.sup.generate_random_elements(seed)
    }
}

impl<B: BorrowedStructure<MontgomeryModuloOddPrimeStructure>> CountableSetSignature
    for MultiplicativeMonoidUnitsStructure<MontgomeryModuloOddPrimeStructure, B>
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> + Clone {
        1..self.monoid().sup.n
    }
}

impl<B: BorrowedStructure<MontgomeryModuloOddPrimeStructure>> FiniteSetSignature
    for MultiplicativeMonoidUnitsStructure<MontgomeryModuloOddPrimeStructure, B>
{
    fn generate_random_elements(&self, seed: u64) -> impl Iterator<Item = Self::Elem> {
        let mut rng = StdRng::seed_from_u64(seed);
        (0..).map(move |_| rng.random_range(1..self.monoid().sup.n))
    }
}

impl FiniteFieldSignature for MontgomeryModuloOddPrimeStructure {
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (Natural::from(self.sup.n), Natural::ONE)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn odd_arithmetic() {
        let ring = MontgomeryModuloOddStructure::new_unchecked(13);

        assert!(ring.equal(&ring.zero(), &0));
        assert!(ring.equal(&ring.add(&1, &2), &3));
        assert!(ring.equal(&ring.add(&3, &10), &0));
        assert!(ring.equal(&ring.add(&10, &11), &8));
        assert!(ring.equal(&ring.neg(&0), &0));
        assert!(ring.equal(&ring.neg(&4), &9));
        assert!(ring.equal(
            &ring.add(&ring.to_montgomery_form(2), &ring.to_montgomery_form(3)),
            &ring.to_montgomery_form(5)
        ));
        assert!(ring.equal(
            &ring.add(&ring.to_montgomery_form(10), &ring.to_montgomery_form(11)),
            &ring.to_montgomery_form(8)
        ));
        assert!(ring.equal(
            &ring.mul(&ring.to_montgomery_form(2), &ring.to_montgomery_form(3)),
            &ring.to_montgomery_form(6)
        ));
        assert!(ring.equal(
            &ring.mul(&ring.to_montgomery_form(5), &ring.to_montgomery_form(6)),
            &ring.to_montgomery_form(4)
        ));
    }

    #[test]
    fn odd_prime_arithmetic() {
        let mut ring = MontgomeryModuloOddPrimeStructure::new_unchecked(13);

        assert!(ring.equal(&ring.zero(), &0));
        assert!(ring.equal(&ring.add(&1, &2), &3));
        assert!(ring.equal(&ring.add(&3, &10), &0));
        assert!(ring.equal(&ring.add(&10, &11), &8));
        assert!(ring.equal(&ring.neg(&0), &0));
        assert!(ring.equal(&ring.neg(&4), &9));
        assert!(ring.equal(
            &ring.add(&ring.to_montgomery_form(2), &ring.to_montgomery_form(3)),
            &ring.to_montgomery_form(5)
        ));
        assert!(ring.equal(
            &ring.add(&ring.to_montgomery_form(10), &ring.to_montgomery_form(11)),
            &ring.to_montgomery_form(8)
        ));
        assert!(ring.equal(
            &ring.mul(&ring.to_montgomery_form(2), &ring.to_montgomery_form(3)),
            &ring.to_montgomery_form(6)
        ));
        assert!(ring.equal(
            &ring.mul(&ring.to_montgomery_form(5), &ring.to_montgomery_form(6)),
            &ring.to_montgomery_form(4)
        ));
        assert_eq!(ring.try_reciprocal(&ring.to_montgomery_form(0)), None);
        assert_eq!(
            ring.try_reciprocal(&ring.to_montgomery_form(1)),
            Some(ring.to_montgomery_form(1))
        );
        assert_eq!(
            ring.try_reciprocal(&ring.to_montgomery_form(3)),
            Some(ring.to_montgomery_form(9))
        );
        ring.populate_inv_cache();
        assert_eq!(
            ring.try_reciprocal(&ring.to_montgomery_form(3)),
            Some(ring.to_montgomery_form(9))
        );
        println!("{:?}", ring.inv_cache);
    }
}
