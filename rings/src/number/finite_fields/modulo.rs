use malachite_nz::{integer::Integer, natural::Natural};

use crate::polynomial::polynomial::*;
use crate::structure::factorization::*;
use crate::structure::structure::*;
use algebraeon_sets::structure::*;

use std::rc::Rc;
use std::{fmt::Display, hash::Hash};

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

impl<const N: usize> From<&malachite_nz::integer::Integer> for Modulo<N> {
    fn from(value: &malachite_nz::integer::Integer) -> Self {
        let value = value % Integer::from(N);
        let value = value.unsigned_abs_ref();
        Self {
            x: malachite_base::num::conversion::traits::ExactInto::exact_into(value),
        }
    }
}
impl<const N: usize> From<malachite_nz::integer::Integer> for Modulo<N> {
    fn from(value: malachite_nz::integer::Integer) -> Self {
        let value = value % Integer::from(N);
        let value = value.unsigned_abs_ref();
        Self {
            x: malachite_base::num::conversion::traits::ExactInto::exact_into(value),
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
            x: modulo(value as isize, N),
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

impl<const N: usize> MetaType for Modulo<N> {
    type Structure = CannonicalStructure<Modulo<N>>;

    fn structure() -> Rc<Self::Structure> {
        CannonicalStructure::new().into()
    }
}

impl<const N: usize> SemiRingStructure for CannonicalStructure<Modulo<N>> {
    fn zero(&self) -> Self::Set {
        Modulo { x: 0 }
    }

    fn one(&self) -> Self::Set {
        if N == 1 {
            Modulo { x: 0 }
        } else {
            Modulo { x: 1 }
        }
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        Modulo { x: (a.x + b.x) % N }
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        Modulo {
            x: (&a.x * &b.x) % N,
        }
    }
}

impl<const N: usize> RingStructure for CannonicalStructure<Modulo<N>> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        if a.x == 0 {
            Modulo { x: 0 }
        } else {
            Modulo { x: N - a.x }
        }
    }
}

macro_rules! impl_field {
    ($N: literal) => {
        impl IntegralDomainStructure for CannonicalStructure<Modulo<$N>> {
            fn div(
                &self,
                top: &Self::Set,
                bot: &Self::Set,
            ) -> Result<Self::Set, RingDivisionError> {
                if bot == &self.zero() {
                    Err(RingDivisionError::DivideByZero)
                } else {
                    let (g, a, b) = xgcd(bot.x, $N);
                    debug_assert_eq!(g as isize, a * bot.x as isize + b * $N as isize);
                    //g = a * bot mod N and g divides N
                    //want ?*bot = top
                    if top.x % g == 0 {
                        Ok(Modulo {
                            x: (modulo(a, $N) * (top.x / g)) % $N,
                        })
                    } else {
                        Err(RingDivisionError::NotDivisible)
                    }
                }
            }
        }
        impl FieldStructure for CannonicalStructure<Modulo<$N>> {}
        impl FiniteUnitsStructure for CannonicalStructure<Modulo<$N>> {
            fn all_units(&self) -> Vec<Modulo<$N>> {
                let mut units = vec![];
                for x in 1..$N {
                    units.push(Modulo { x });
                }
                units
            }
        }
        impl FiniteFieldStructure for CannonicalStructure<Modulo<$N>> {
            fn characteristic_and_power(&self) -> (Natural, Natural) {
                (Natural::from($N as usize), Natural::from(1u8))
            }
        }
        impl UniqueFactorizationStructure for PolynomialStructure<CannonicalStructure<Modulo<$N>>> {
            fn factor(&self, p: &Self::Set) -> Option<Factored<Self>> {
                Some(
                    self.factorize_monic(p)?
                        .factorize_squarefree()
                        .factorize_distinct_degree().factorize_cantor_zassenhaus(),
                )
            }
        }
    };
}

impl_field!(2);
impl_field!(3);
impl_field!(5);
impl_field!(7);
impl_field!(11);
impl_field!(13);
impl_field!(17);
impl_field!(19);
impl_field!(23);
impl_field!(29);
impl_field!(31);
impl_field!(37);
impl_field!(41);
impl_field!(43);
impl_field!(47);
impl_field!(53);
impl_field!(59);
impl_field!(61);
impl_field!(67);
impl_field!(71);
impl_field!(73);
impl_field!(79);
impl_field!(83);
impl_field!(89);
impl_field!(97);
impl_field!(101);
impl_field!(103);
impl_field!(107);
impl_field!(109);
impl_field!(113);
impl_field!(127);
impl_field!(131);
impl_field!(137);
impl_field!(139);
impl_field!(149);
impl_field!(151);
impl_field!(157);
impl_field!(163);
impl_field!(167);
impl_field!(173);
impl_field!(179);
impl_field!(181);
impl_field!(191);
impl_field!(193);
impl_field!(197);
impl_field!(199);
impl_field!(211);
impl_field!(223);
impl_field!(227);
impl_field!(229);
impl_field!(233);
impl_field!(239);
impl_field!(241);
impl_field!(251);
impl_field!(257);
impl_field!(263);
impl_field!(269);
impl_field!(271);
impl_field!(277);
impl_field!(281);
impl_field!(283);
impl_field!(293);
impl_field!(307);
impl_field!(311);
impl_field!(313);
impl_field!(317);
impl_field!(331);
impl_field!(337);
impl_field!(347);
impl_field!(349);
impl_field!(353);
impl_field!(359);
impl_field!(367);
impl_field!(373);
impl_field!(379);
impl_field!(383);
impl_field!(389);
impl_field!(397);
impl_field!(401);
impl_field!(409);
impl_field!(419);
impl_field!(421);
impl_field!(431);
impl_field!(433);
impl_field!(439);
impl_field!(443);
impl_field!(449);
impl_field!(457);
impl_field!(461);
impl_field!(463);
impl_field!(467);
impl_field!(479);
impl_field!(487);
impl_field!(491);
impl_field!(499);
impl_field!(503);
impl_field!(509);
impl_field!(521);
impl_field!(523);
impl_field!(541);
impl_field!(547);
impl_field!(557);
impl_field!(563);
impl_field!(569);
impl_field!(571);
impl_field!(577);
impl_field!(587);
impl_field!(593);
impl_field!(599);
impl_field!(601);
impl_field!(607);
impl_field!(613);
impl_field!(617);
impl_field!(619);
impl_field!(631);
impl_field!(641);
impl_field!(643);
impl_field!(647);
impl_field!(653);
impl_field!(659);
impl_field!(661);
impl_field!(673);
impl_field!(677);
impl_field!(683);
impl_field!(691);
impl_field!(701);
impl_field!(709);
impl_field!(719);
impl_field!(727);
impl_field!(733);
impl_field!(739);
impl_field!(743);
impl_field!(751);
impl_field!(757);
impl_field!(761);
impl_field!(769);
impl_field!(773);
impl_field!(787);
impl_field!(797);
impl_field!(809);
impl_field!(811);
impl_field!(821);
impl_field!(823);
impl_field!(827);
impl_field!(829);
impl_field!(839);
impl_field!(853);
impl_field!(857);
impl_field!(859);
impl_field!(863);
impl_field!(877);
impl_field!(881);
impl_field!(883);
impl_field!(887);
impl_field!(907);
impl_field!(911);
impl_field!(919);
impl_field!(929);
impl_field!(937);
impl_field!(941);
impl_field!(947);
impl_field!(953);
impl_field!(967);
impl_field!(971);
impl_field!(977);
impl_field!(983);
impl_field!(991);
impl_field!(997);
//TODO: add more or generate these with a proc macro

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
}
