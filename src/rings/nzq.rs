use super::ring::*;
use malachite_base::num::arithmetic::traits::{DivMod, UnsignedAbs};
use malachite_nz::{integer::Integer, natural::Natural};
use malachite_q::Rational;

impl ComRing for Integer {
    fn zero() -> Self {
        Self::from(0)
    }
    fn one() -> Self {
        Self::from(1)
    }

    fn neg_mut(&mut self) {
        *self *= Integer::from(-1)
    }
    fn neg(self) -> Self {
        -self
    }

    fn add_mut(&mut self, x: &Self) {
        *self += x;
    }
    fn add(a: Self, b: Self) -> Self {
        a + b
    }
    fn add_ref(a: Self, b: &Self) -> Self {
        a + b
    }
    fn add_refs(a: &Self, b: &Self) -> Self {
        a + b
    }

    fn mul_mut(&mut self, x: &Self) {
        *self *= x;
    }
    fn mul(a: Self, b: Self) -> Self {
        a * b
    }
    fn mul_ref(a: Self, b: &Self) -> Self {
        a * b
    }
    fn mul_refs(a: &Self, b: &Self) -> Self {
        a * b
    }

    fn div(a: Self, b: Self) -> Result<Self, RingOppErr> {
        match <Self as EuclideanDomain>::quorem(a, b) {
            Ok((q, r)) => {
                if r == Self::zero() {
                    Ok(q)
                } else {
                    Err(RingOppErr::NotDivisible)
                }
            }
            Err(e) => Err(e),
        }
    }

    fn div_lref(a: &Self, b: Self) -> Result<Self, RingOppErr> {
        match <Self as EuclideanDomain>::quorem_lref(a, b) {
            Ok((q, r)) => {
                if r == Self::zero() {
                    Ok(q)
                } else {
                    Err(RingOppErr::NotDivisible)
                }
            }
            Err(e) => Err(e),
        }
    }

    fn div_rref(a: Self, b: &Self) -> Result<Self, RingOppErr> {
        match <Self as EuclideanDomain>::quorem_rref(a, b) {
            Ok((q, r)) => {
                if r == Self::zero() {
                    Ok(q)
                } else {
                    Err(RingOppErr::NotDivisible)
                }
            }
            Err(e) => Err(e),
        }
    }

    fn div_refs(a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        match <Self as EuclideanDomain>::quorem_refs(a, b) {
            Ok((q, r)) => {
                if r == Self::zero() {
                    Ok(q)
                } else {
                    Err(RingOppErr::NotDivisible)
                }
            }
            Err(e) => Err(e),
        }
    }
}

impl IntegralDomain for Integer {}

impl FavoriteAssociate for Integer {
    fn factor_fav_assoc(self) -> (Self, Self) {
        if self == 0 {
            (Self::one(), Self::zero())
        } else if self < 0 {
            (Integer::from(-1), self.neg())
        } else {
            (Integer::from(1), self)
        }
    }
}

impl EuclideanDomain for Integer {
    fn norm(&self) -> Option<Natural> {
        if self == &Integer::from(0) {
            None
        } else {
            Some(self.unsigned_abs())
        }
    }

    fn quorem(a: Self, b: Self) -> Result<(Self, Self), RingOppErr> {
        if b == Integer::from(0) {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(a.div_mod(b.clone()))
        }
    }
}

impl ComRing for Rational {
    fn zero() -> Self {
        Self::from(0)
    }
    fn one() -> Self {
        Self::from(1)
    }

    fn neg_mut(&mut self) {
        *self *= Rational::from(-1);
    }
    fn neg_ref(&self) -> Self {
        -self
    }
    fn neg(self) -> Self {
        -self
    }

    fn add_mut(&mut self, x: &Self) {
        *self += x;
    }
    fn add(a: Self, b: Self) -> Self {
        a + b
    }
    fn add_ref(a: Self, b: &Self) -> Self {
        a + b
    }
    fn add_refs(a: &Self, b: &Self) -> Self {
        a + b
    }

    fn mul_mut(&mut self, x: &Self) {
        *self *= x;
    }
    fn mul(a: Self, b: Self) -> Self {
        a * b
    }
    fn mul_ref(a: Self, b: &Self) -> Self {
        a * b
    }
    fn mul_refs(a: &Self, b: &Self) -> Self {
        a * b
    }

    fn div(a: Self, b: Self) -> Result<Self, RingOppErr> {
        if b == Rational::from(0) {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(a / b)
        }
    }
}
impl IntegralDomain for Rational {}
impl Field for Rational {
    // fn inv(a: Self) -> Result<Self, OppErr> {
    //     if a.numerator_ref() == &Natural::from(0u8) {
    //         Err(OppErr::DivideByZero)
    //     } else {
    //         Ok(a.reciprocal())
    //     }
    // }
}

impl FieldOfFractions<Integer> for Rational {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_int() {
        //happy div
        {
            let a = Integer::from(18);
            let b = Integer::from(6);
            let c = Integer::div(a, b);
            match c {
                Ok(_) => {}
                Err(_e) => panic!(),
            }
        }

        //sad div
        {
            let a = Integer::from(18);
            let b = Integer::from(7);
            let c = Integer::div(a, b);
            match c {
                Ok(_) => panic!(),
                Err(e) => match e {
                    RingOppErr::DivideByZero => panic!(),
                    RingOppErr::NotDivisible => {}
                },
            }
        }

        //euclidean div
        {
            let a = Integer::from(18);
            let b = Integer::from(7);
            let (q, r) = Integer::quorem_refs(&a, &b).unwrap();
            assert!(r.norm() < b.norm());
            assert_eq!(a, b * q + r);
        }

        //xgcd
        {
            let a = Integer::from(31);
            let b = Integer::from(57);
            let (g, x, y) = Integer::xgcd(a.clone(), b.clone());
            assert_eq!(x * a + y * b, g);
        }
    }
}
