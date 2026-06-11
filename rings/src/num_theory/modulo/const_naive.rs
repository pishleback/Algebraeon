use crate::structure::*;
use algebraeon_macros::repeat_small_primes;
use algebraeon_structures::*;
use std::{borrow::Cow, fmt::Display, hash::Hash};

fn xgcd(mut x: usize, mut y: usize) -> (usize, isize, isize) {
    let mut pa = 1;
    let mut a = 0;
    let mut pb = 0;
    let mut b = 1;
    while y != 0 {
        let (q, r) = (x / y, x % y);
        let new_a = pa - q as isize * a;
        (a, pa) = (new_a, a);
        let new_b = pb - q as isize * b;
        (b, pb) = (new_b, b);
        (x, y) = (y, r);
    }
    (x, pa, pb)
}

fn modulo(a: isize, n: usize) -> usize {
    let mut a = a % n as isize;
    if a < 0 {
        a += n as isize;
    }
    a as usize
}

#[derive(Clone)]
pub struct Modulo<const N: usize> {
    x: usize,
}

impl<const N: usize> Modulo<N> {
    pub const fn new(x: usize) -> Self {
        Self { x }
    }
}

impl<const N: usize> std::fmt::Debug for Modulo<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Modulo")
            .field("x", &self.x)
            .field("N", &N)
            .finish()
    }
}

impl<const N: usize> Modulo<N> {
    pub fn lift_nat(&self) -> Natural {
        Natural::from(self.x)
    }
    pub fn lift_int(&self) -> Integer {
        Integer::from(self.x)
    }
}

impl<const N: usize> From<&Integer> for Modulo<N> {
    fn from(value: &Integer) -> Self {
        let value = value % Integer::from(N);
        let value = value.abs();
        Self {
            x: value.try_into().unwrap(),
        }
    }
}
impl<const N: usize> From<Integer> for Modulo<N> {
    fn from(value: Integer) -> Self {
        let value = value % Integer::from(N);
        let value = value.abs();
        Self {
            x: value.try_into().unwrap(),
        }
    }
}
impl<const N: usize> From<usize> for Modulo<N> {
    fn from(value: usize) -> Self {
        Self { x: value % N }
    }
}
impl<const N: usize> From<isize> for Modulo<N> {
    fn from(value: isize) -> Self {
        Self {
            x: modulo(value, N),
        }
    }
}
impl<const N: usize> From<i32> for Modulo<N> {
    fn from(value: i32) -> Self {
        Self {
            x: modulo(value as isize, N),
        }
    }
}

impl<const N: usize> From<Modulo<N>> for usize {
    fn from(value: Modulo<N>) -> Self {
        value.x
    }
}

impl<const N: usize> PartialEq for Modulo<N> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x
    }
}
impl<const N: usize> Eq for Modulo<N> {}
impl<const N: usize> Hash for Modulo<N> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.x.hash(state);
    }
}

impl<const N: usize> Display for Modulo<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.x)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModuloCanonicalStructure<const N: usize> {}

impl<const N: usize> Signature for ModuloCanonicalStructure<N> {}

impl<const N: usize> SetSignature for ModuloCanonicalStructure<N> {
    type Elem = Modulo<N>;

    fn validate_element(&self, _x: &Self::Elem) -> Result<(), String> {
        Ok(())
    }
}

impl<const N: usize> MetaType for Modulo<N> {
    type Signature = ModuloCanonicalStructure<N>;

    fn structure() -> Self::Signature {
        ModuloCanonicalStructure {}
    }
}

impl<const N: usize> EqSignature for ModuloCanonicalStructure<N> {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        a == b
    }
}

impl<const N: usize> ToStringSignature for ModuloCanonicalStructure<N> {
    fn to_string(&self, elem: &Self::Elem) -> String {
        format!("{}", elem)
    }
}

impl<const N: usize> RinglikeSpecializationSignature for ModuloCanonicalStructure<N> {}

impl<const N: usize> ZeroSignature for ModuloCanonicalStructure<N> {
    fn zero(&self) -> Self::Elem {
        Modulo { x: 0 }
    }
}

impl<const N: usize> AdditionSignature for ModuloCanonicalStructure<N> {
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        Modulo { x: (a.x + b.x) % N }
    }
}

impl<const N: usize> CancellativeAdditionSignature for ModuloCanonicalStructure<N> {
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<const N: usize> TryNegateSignature for ModuloCanonicalStructure<N> {
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<const N: usize> AdditiveMonoidSignature for ModuloCanonicalStructure<N> {}

impl<const N: usize> AdditiveGroupSignature for ModuloCanonicalStructure<N> {
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        if a.x == 0 {
            Modulo { x: 0 }
        } else {
            Modulo { x: N - a.x }
        }
    }
}

impl<const N: usize> OneSignature for ModuloCanonicalStructure<N> {
    fn one(&self) -> Self::Elem {
        if N == 1 {
            Modulo { x: 0 }
        } else {
            Modulo { x: 1 }
        }
    }
}

impl<const N: usize> MultiplicationSignature for ModuloCanonicalStructure<N> {
    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        Modulo { x: (a.x * b.x) % N }
    }
}

impl<const N: usize> CommutativeMultiplicationSignature for ModuloCanonicalStructure<N> {}

impl<const N: usize> MultiplicativeMonoidSignature for ModuloCanonicalStructure<N> {}

impl<const N: usize> MultiplicativeAbsorptionMonoidSignature for ModuloCanonicalStructure<N> {}

impl<const N: usize> LeftDistributiveMultiplicationOverAddition for ModuloCanonicalStructure<N> {}

impl<const N: usize> RightDistributiveMultiplicationOverAddition for ModuloCanonicalStructure<N> {}

impl<const N: usize> SemiRingSignature for ModuloCanonicalStructure<N> {}

impl<const N: usize> RingSignature for ModuloCanonicalStructure<N> {}

impl<const N: usize> CharacteristicSignature for ModuloCanonicalStructure<N> {
    fn characteristic(&self) -> Natural {
        Natural::from(N)
    }
}

impl<const N: usize> TryReciprocalSignature for ModuloCanonicalStructure<N> {
    fn try_reciprocal(&self, x: &Self::Elem) -> Option<Self::Elem> {
        if x == &self.zero() {
            None
        } else {
            let (g, a, b) = xgcd(x.x, N);
            debug_assert_eq!(g as isize, a * x.x as isize + b * N as isize);
            if g == 1 {
                Some(Modulo { x: modulo(a, N) })
            } else {
                None
            }
        }
    }
}

impl<const N: usize> QuotientSetSignature<IntegerCanonicalStructure>
    for ModuloCanonicalStructure<N>
{
    fn pre_quotient_set(&self) -> &IntegerCanonicalStructure {
        Integer::structure_ref()
    }

    fn project(&self, x: Integer) -> Self::Elem {
        x.into()
    }

    fn project_ref(&self, x: &Integer) -> Self::Elem {
        x.into()
    }

    fn unproject(&self, x: Self::Elem) -> Integer {
        x.lift_int()
    }

    fn unproject_ref(&self, x: &Self::Elem) -> Integer {
        x.lift_int()
    }
}

impl<const N: usize> QuotientRingSignature<IntegerCanonicalStructure>
    for ModuloCanonicalStructure<N>
{
}

impl<const N: usize> QuotientRingGetPrincipalIdealSignature<IntegerCanonicalStructure>
    for ModuloCanonicalStructure<N>
{
    fn modulus<'a>(&'a self) -> Cow<'a, Integer> {
        Cow::Owned(Integer::from(N))
    }
}

impl<const N: usize> CountableSetSignature for ModuloCanonicalStructure<N> {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> + Clone {
        (0..N).map(Modulo::new)
    }
}

impl<const N: usize> FiniteSetSignature for ModuloCanonicalStructure<N> {
    fn size(&self) -> usize {
        N
    }
}

macro_rules! impl_field {
    ($N: literal) => {
        impl CancellativeMultiplicationSignature for ModuloCanonicalStructure<$N> {
            fn try_divide(&self, top: &Self::Elem, bot: &Self::Elem) -> Option<Self::Elem> {
                Some(self.mul(top, &self.try_reciprocal(bot)?))
            }
        }

        impl MultiplicativeIntegralMonoidSignature for ModuloCanonicalStructure<$N> {}

        impl IntegralDomainSignature for ModuloCanonicalStructure<$N> {}

        impl FieldSignature for ModuloCanonicalStructure<$N> {}

        impl<B: BorrowedStructure<ModuloCanonicalStructure<$N>>> CountableSetSignature
            for MultiplicativeMonoidUnitsStructure<ModuloCanonicalStructure<$N>, B>
        {
            fn generate_all_elements(
                &self,
            ) -> impl std::iter::Iterator<Item = Modulo<$N>> + std::clone::Clone {
                (1..$N).map(|i| Modulo { x: i })
            }
        }

        impl<B: BorrowedStructure<ModuloCanonicalStructure<$N>>> FiniteSetSignature
            for MultiplicativeMonoidUnitsStructure<ModuloCanonicalStructure<$N>, B>
        {
        }

        impl FiniteFieldSignature for ModuloCanonicalStructure<$N> {
            fn characteristic_and_power(&self) -> (Natural, Natural) {
                (Natural::from($N as usize), Natural::from(1u8))
            }
        }
    };
}

repeat_small_primes!(20, p =>
    impl_field!(p);
);

impl<
    B: BorrowedStructure<IntegerCanonicalStructure>,
    BE: BorrowedStructure<EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, true>>,
> CountableSetSignature
    for MultiplicativeMonoidUnitsStructure<
        EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, true>,
        BE,
    >
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> + Clone {
        self.list_all_elements().into_iter()
    }
}

impl<
    B: BorrowedStructure<IntegerCanonicalStructure>,
    BE: BorrowedStructure<EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, true>>,
> FiniteSetSignature
    for MultiplicativeMonoidUnitsStructure<
        EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, true>,
        BE,
    >
{
    fn list_all_elements(&self) -> Vec<Self::Elem> {
        let mut units = vec![];
        let mut u = Integer::from(1);
        while u < Abs::abs(self.monoid().modulus().as_ref()) {
            units.push(u.clone());
            u += Integer::ONE;
        }
        units
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>> FiniteFieldSignature
    for EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, true>
{
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        (Abs::abs(self.modulus().as_ref()), Natural::ONE)
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>, const IS_FIELD: bool> CountableSetSignature
    for EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, IS_FIELD>
{
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> + Clone {
        (0usize..)
            .map(Integer::from)
            .take_while(|n| n < self.modulus().as_ref())
    }
}

impl<B: BorrowedStructure<IntegerCanonicalStructure>, const IS_FIELD: bool> FiniteSetSignature
    for EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, B, IS_FIELD>
{
    fn size(&self) -> usize {
        Abs::abs(self.modulus().as_ref()).try_into().unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gcd_lcm() {
        for x in 0..100 {
            for y in 0..100 {
                let (g, a, b) = xgcd(x, y);
                assert_eq!(g as isize, a * x as isize + b * y as isize);
            }
        }
    }

    #[test]
    fn count_elements() {
        assert_eq!(Modulo::<26>::structure().list_all_elements().len(), 26);
    }
}
