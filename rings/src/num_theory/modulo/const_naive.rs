use crate::structure::*;
use algebraeon_macros::repeat_small_primes;
use algebraeon_structures::*;
use std::{borrow::Cow, cmp::Ordering, fmt::Display, hash::Hash, rc::Rc};

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

    fn validate_element(self: &Rc<Self>, _x: &Self::Elem) -> Result<(), String> {
        Ok(())
    }
}

impl<const N: usize> MetaType for Modulo<N> {
    type Signature = ModuloCanonicalStructure<N>;

    fn structure() -> Rc<Self::Signature> {
        ModuloCanonicalStructure {}.into()
    }
}

impl<const N: usize> EqSignature for ModuloCanonicalStructure<N> {
    fn equal(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> bool {
        a == b
    }
}

impl<const N: usize> PartialOrdSignature for ModuloCanonicalStructure<N> {
    fn partial_cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        a.x.partial_cmp(&b.x)
    }
}

impl<const N: usize> OrdSignature for ModuloCanonicalStructure<N> {
    fn cmp(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        a.x.cmp(&b.x)
    }
}

impl<const N: usize> ToStringSignature for ModuloCanonicalStructure<N> {
    fn to_string(self: &Rc<Self>, elem: &Self::Elem) -> String {
        format!("{}", elem)
    }
}

impl<const N: usize> RinglikeSpecializationSignature for ModuloCanonicalStructure<N> {}

impl<const N: usize> ZeroSignature for ModuloCanonicalStructure<N> {
    fn zero(self: &Rc<Self>) -> Self::Elem {
        Modulo { x: 0 }
    }
}

impl<const N: usize> AdditionSignature for ModuloCanonicalStructure<N> {
    fn add(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        Modulo { x: (a.x + b.x) % N }
    }
}

impl<const N: usize> CancellativeAdditionSignature for ModuloCanonicalStructure<N> {
    fn try_sub(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.sub(a, b))
    }
}

impl<const N: usize> TryNegateSignature for ModuloCanonicalStructure<N> {
    fn try_neg(self: &Rc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl<const N: usize> AdditiveMonoidSignature for ModuloCanonicalStructure<N> {}

impl<const N: usize> AdditiveGroupSignature for ModuloCanonicalStructure<N> {
    fn neg(self: &Rc<Self>, a: &Self::Elem) -> Self::Elem {
        if a.x == 0 {
            Modulo { x: 0 }
        } else {
            Modulo { x: N - a.x }
        }
    }
}

impl<const N: usize> OneSignature for ModuloCanonicalStructure<N> {
    fn one(self: &Rc<Self>) -> Self::Elem {
        if N == 1 {
            Modulo { x: 0 }
        } else {
            Modulo { x: 1 }
        }
    }
}

impl<const N: usize> MultiplicationSignature for ModuloCanonicalStructure<N> {
    fn mul(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
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
    fn characteristic(self: &Rc<Self>) -> Natural {
        Natural::from(N)
    }
}

impl<const N: usize> TryReciprocalSignature for ModuloCanonicalStructure<N> {
    fn try_reciprocal(self: &Rc<Self>, x: &Self::Elem) -> Option<Self::Elem> {
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
    fn pre_quotient_set(self: &Rc<Self>) -> Rc<IntegerCanonicalStructure> {
        Integer::structure()
    }

    fn project(self: &Rc<Self>, x: Integer) -> Self::Elem {
        x.into()
    }

    fn project_ref(self: &Rc<Self>, x: &Integer) -> Self::Elem {
        x.into()
    }

    fn unproject(self: &Rc<Self>, x: Self::Elem) -> Integer {
        x.lift_int()
    }

    fn unproject_ref(self: &Rc<Self>, x: &Self::Elem) -> Integer {
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
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        (0..N).map(Modulo::new)
    }
}

impl<const N: usize> FiniteSetSignature for ModuloCanonicalStructure<N> {
    fn size(self: &Rc<Self>) -> Natural {
        Natural::from(N)
    }
}

impl<const N: usize> OrderedFiniteSetSignature for ModuloCanonicalStructure<N> {
    fn list_all_elements_ordered(self: &Rc<Self>) -> Vec<Self::Elem> {
        (0..N).map(|x| Modulo { x }).collect()
    }

    fn element_to_enumeration(self: &Rc<Self>, elem: &Self::Elem) -> Natural {
        elem.x.into()
    }

    fn enumeration_to_element(self: &Rc<Self>, num: &Natural) -> Option<Self::Elem> {
        let x: usize = num.try_into().ok()?;
        if x < N { Some(Modulo { x }) } else { None }
    }
}

macro_rules! impl_field {
    ($N: literal) => {
        impl CancellativeMultiplicationSignature for ModuloCanonicalStructure<$N> {
            fn try_divide(
                self: &Rc<Self>,
                top: &Self::Elem,
                bot: &Self::Elem,
            ) -> Option<Self::Elem> {
                Some(self.mul(top, &self.try_reciprocal(bot)?))
            }
        }

        impl MultiplicativeIntegralMonoidSignature for ModuloCanonicalStructure<$N> {}

        impl IntegralDomainSignature for ModuloCanonicalStructure<$N> {}

        impl FieldSignature for ModuloCanonicalStructure<$N> {}

        impl CountableSetSignature
            for MultiplicativeMonoidUnitsStructure<ModuloCanonicalStructure<$N>>
        {
            fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
                (1..$N).map(|i| Modulo { x: i })
            }
        }

        impl FiniteSetSignature
            for MultiplicativeMonoidUnitsStructure<ModuloCanonicalStructure<$N>>
        {
        }

        impl FiniteFieldSignature for ModuloCanonicalStructure<$N> {
            fn characteristic_and_power(self: &Rc<Self>) -> (Natural, Natural) {
                (Natural::from($N as usize), Natural::from(1u8))
            }
        }
    };
}

repeat_small_primes!(20, p =>
    impl_field!(p);
);

impl CountableSetSignature
    for MultiplicativeMonoidUnitsStructure<
        EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, true>,
    >
{
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        self.list_all_elements().into_iter()
    }
}

impl FiniteSetSignature
    for MultiplicativeMonoidUnitsStructure<
        EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, true>,
    >
{
    fn list_all_elements(self: &Rc<Self>) -> Vec<Self::Elem> {
        let mut units = vec![];
        let mut u = Integer::from(1);
        while u < Abs::abs(self.monoid().modulus().as_ref()) {
            units.push(u.clone());
            u += Integer::ONE;
        }
        units
    }
}

impl FiniteFieldSignature for EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, true> {
    fn characteristic_and_power(self: &Rc<Self>) -> (Natural, Natural) {
        (Abs::abs(self.modulus().as_ref()), Natural::ONE)
    }
}

impl<const IS_FIELD: bool> CountableSetSignature
    for EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, IS_FIELD>
{
    fn generate_all_elements(self: Rc<Self>) -> impl Iterator<Item = Self::Elem> {
        let modulus = self.modulus().as_ref().clone();
        (0usize..)
            .map(Integer::from)
            .take_while(move |n| n < &modulus)
    }
}

impl<const IS_FIELD: bool> FiniteSetSignature
    for EuclideanRemainderQuotientStructure<IntegerCanonicalStructure, IS_FIELD>
{
    fn size(self: &Rc<Self>) -> Natural {
        Abs::abs(self.modulus().as_ref())
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

    #[test]
    fn enumeration() {
        algebraeon_structures::assert_enumerated_ord_finite_set!(Modulo::<26>::structure(), 26);
    }
}
