use std::fmt::Debug;

use malachite_base::num::{
    arithmetic::traits::{DivRem, UnsignedAbs},
    logic::traits::BitIterable,
};
use malachite_nz::{integer::Integer, natural::Natural};

#[derive(Debug)]
pub enum RingOppErr {
    DivideByZero,
    NotDivisible,
}

pub trait ComRing: Sized + Clone + PartialEq + Eq + Debug {
    fn zero() -> Self;
    fn one() -> Self;
    fn neg_mut(&mut self);
    fn neg_ref(&self) -> Self {
        self.clone().neg()
    }
    fn neg(mut self) -> Self {
        self.neg_mut();
        self
    }

    fn add_mut(&mut self, x: &Self);
    fn add(mut a: Self, b: Self) -> Self {
        a.add_mut(&b);
        a
    }
    fn add_ref(mut a: Self, b: &Self) -> Self {
        a.add_mut(b);
        a
    }
    fn add_refs(a: &Self, b: &Self) -> Self {
        let mut new_a = a.clone();
        new_a.add_mut(b);
        new_a
    }

    fn mul_mut(&mut self, x: &Self);
    fn mul(mut a: Self, b: Self) -> Self {
        a.mul_mut(&b);
        a
    }
    fn mul_ref(mut a: Self, b: &Self) -> Self {
        a.mul_mut(b);
        a
    }
    fn mul_refs(a: &Self, b: &Self) -> Self {
        let mut new_a = a.clone();
        new_a.mul_mut(b);
        new_a
    }

    fn div(a: Self, b: Self) -> Result<Self, RingOppErr>;
    fn div_lref(a: &Self, b: Self) -> Result<Self, RingOppErr> {
        Self::div(a.clone(), b)
    }
    fn div_rref(a: Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::div(a, b.clone())
    }
    fn div_refs(a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::div(a.clone(), b.clone())
    }

    fn divisible(a : Self, b : Self) -> bool {
        match Self::div(a, b) {
            Ok(_q) => true,
            Err(RingOppErr::NotDivisible) => false,
            Err(RingOppErr::DivideByZero) => false,
        }
    }

    fn sum(elems: Vec<Self>) -> Self {
        let mut ans = Self::zero();
        for elem in elems {
            ans = Self::add(ans, elem);
        }
        ans
    }

    fn product(elems: Vec<Self>) -> Self {
        let mut ans = Self::one();
        for elem in elems {
            ans = Self::mul(ans, elem);
        }
        ans
    }

    fn pow(&self, n: &Natural) -> Self {
        if *n == 0 {
            Self::one()
        } else if *n == 1 {
            self.clone()
        } else {
            debug_assert!(*n >= 2);
            let (q, r) = n.div_rem(Natural::from(2u8));
            Self::mul(self.pow(&q), self.pow(&(&q + r)))
        }
    }

    fn from_int(x: &Integer) -> Self {
        if *x < 0 {
            Self::from_int(&-x).neg()
        } else if *x == 0 {
            Self::zero()
        } else if *x == 1 {
            Self::one()
        } else {
            let two = Self::add(Self::one(), Self::one());
            debug_assert!(*x >= 2);
            let bits: Vec<bool> = x.unsigned_abs().bits().collect();
            let mut ans = Self::zero();
            let mut v = Self::one();
            for i in 0..bits.len() {
                if bits[i] {
                    ans.add_mut(&v);
                }
                v.mul_mut(&two);
            }
            ans
        }
    }

    fn is_unit(self) -> bool {
        match Self::div(Self::one(), self) {
            Ok(_inv) => true,
            Err(RingOppErr::DivideByZero) => false,
            Err(RingOppErr::NotDivisible) => false,
            // Err(_) => panic!(),
        }
    }

    fn inv(self) -> Result<Self, RingOppErr> {
        Self::div(Self::one(), self)
    }

    fn inv_ref(a: &Self) -> Result<Self, RingOppErr> {
        Self::div_rref(Self::one(), a)
    }
}

pub trait IntegralDomain: ComRing {
    //promise that mul(a, b) == 0 implies a == 0 or b == 0
}

pub trait FavoriteAssociate: IntegralDomain {
    //For associate class of elements, choose a unique representative
    //write self=unit*assoc and return (unit, assoc)
    //0 is required to return (1, 0)

    fn factor_fav_assoc(self) -> (Self, Self);
    fn factor_fav_assoc_ref(&self) -> (Self, Self) {
        self.clone().factor_fav_assoc()
    }
}

pub trait GCDDomain: ComRing + FavoriteAssociate {
    //any gcds should be the standard associate representative
    fn gcd(x: Self, y: Self) -> Self;
    fn gcd_list(elems : &Vec<Self>) -> Self {
        let mut ans = Self::zero();
        for x in elems {
            ans = Self::gcd(ans, x.clone());
        }
        ans
    }
    fn xgcd(x: Self, y: Self) -> (Self, Self, Self);
}

// pub trait PrincipalIdealDomain: GCDDomain {
// }

pub trait EuclideanDomain : IntegralDomain {
    //should return None for 0, and Some(norm) for everything else
    fn norm(&self) -> Option<Natural>;
    fn quorem(a: Self, b: Self) -> Result<(Self, Self), RingOppErr>;
    fn quorem_lref(a: &Self, b: Self) -> Result<(Self, Self), RingOppErr> {
        Self::quorem(a.clone(), b)
    }
    fn quorem_rref(a: Self, b: &Self) -> Result<(Self, Self), RingOppErr> {
        Self::quorem(a, b.clone())
    }
    fn quorem_refs(a: &Self, b: &Self) -> Result<(Self, Self), RingOppErr> {
        Self::quorem(a.clone(), b.clone())
    }

    fn quo(a: Self, b: Self) -> Result<Self, RingOppErr> {
        match Self::quorem(a, b) {
            Ok((q, _r)) => Ok(q),
            Err(e) => Err(e),
        }
    }
    fn quo_lref(a: &Self, b: Self) -> Result<Self, RingOppErr> {
        Self::quo(a.clone(), b)
    }
    fn quo_rref(a: Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::quo(a, b.clone())
    }
    fn quo_refs(a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::quo(a.clone(), b.clone())
    }

    fn rem(a: Self, b: Self) -> Result<Self, RingOppErr> {
        match Self::quorem(a, b) {
            Ok((_q, r)) => Ok(r),
            Err(e) => Err(e),
        }
    }
    fn rem_lref(a: &Self, b: Self) -> Result<Self, RingOppErr> {
        Self::rem(a.clone(), b)
    }
    fn rem_rref(a: Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::rem(a, b.clone())
    }
    fn rem_refs(a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        Self::rem(a.clone(), b.clone())
    }
}

impl<R: EuclideanDomain + FavoriteAssociate> GCDDomain for R {
    fn gcd(mut x: Self, mut y: Self) -> Self {
        //Euclidean algorithm
        while y != Self::zero() {
            let r = Self::rem_rref(x, &y).unwrap();
            (x, y) = (y, r)
        };
        let (_unit, assoc) = x.factor_fav_assoc();
        assoc
    }

    fn xgcd(mut x: Self, mut y: Self) -> (Self, Self, Self) {
        let mut pa = Self::one();
        let mut a = Self::zero();
        let mut pb = Self::zero();
        let mut b = Self::one();

        while y != Self::zero() {
            let (q, r) = Self::quorem_rref(x, &y).unwrap();
            let new_a = Self::add(pa, Self::mul_refs(&q, &a).neg());
            (a, pa) = (new_a, a);
            let new_b = Self::add(pb, Self::mul_ref(q, &b).neg());
            (b, pb) = (new_b, b);
            (x, y) = (y, r);
        }
        let (unit, ass_x) = x.factor_fav_assoc();
        // g = u*g_ass
        // g = xa+by
        // xa+by=u*g_ass
        debug_assert!(unit.clone().is_unit());
        (ass_x, R::div_rref(pa, &unit).unwrap(), R::div(pb, unit).unwrap())
    }
}

// impl<R: EuclideanDomain> PrincipalIdealDomain for R {

// }

pub trait Field: IntegralDomain {
    //promise that a/b always works, except unless b=0.
    //in other words, a/b must not return not divisible
}

impl<F: Field> FavoriteAssociate for F {
    fn factor_fav_assoc(self) -> (Self, Self) {
        (self, Self::one())
    }
}

impl<F: Field> EuclideanDomain for F {
    fn norm(&self) -> Option<Natural> {
        if self == &Self::zero() {
            None
        } else {
            Some(Natural::from(1u8))
        }
    }

    fn quorem(a: Self, b: Self) -> Result<(Self, Self), RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok((Self::div(a, b).unwrap(), Self::zero()))
        }
    }

    fn quorem_lref(a: &Self, b: Self) -> Result<(Self, Self), RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok((Self::div_lref(a, b).unwrap(), Self::zero()))
        }
    }

    fn quorem_rref(a: Self, b: &Self) -> Result<(Self, Self), RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok((Self::div_rref(a, b).unwrap(), Self::zero()))
        }
    }

    fn quorem_refs(a: &Self, b: &Self) -> Result<(Self, Self), RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok((Self::div_refs(a, b).unwrap(), Self::zero()))
        }
    }

    fn quo(a: Self, b: Self) -> Result<Self, RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::div(a, b).unwrap())
        }
    }

    fn quo_lref(a: &Self, b: Self) -> Result<Self, RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::div_lref(a, b).unwrap())
        }
    }

    fn quo_rref(a: Self, b: &Self) -> Result<Self, RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::div_rref(a, b).unwrap())
        }
    }

    fn quo_refs(a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::div_refs(a, b).unwrap())
        }
    }

    fn rem(_a: Self, b: Self) -> Result<Self, RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::zero())
        }
    }

    fn rem_lref(_a: &Self, b: Self) -> Result<Self, RingOppErr> {
        if b == Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::zero())
        }
    }

    fn rem_rref(_a: Self, b: &Self) -> Result<Self, RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::zero())
        }
    }

    fn rem_refs(_a: &Self, b: &Self) -> Result<Self, RingOppErr> {
        if b == &Self::zero() {
            Err(RingOppErr::DivideByZero)
        } else {
            Ok(Self::zero())
        }
    }
}

pub trait FieldOfFractions<R: IntegralDomain>: Field + From<R> {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ergonomic<R: ComRing> {
    elem: R,
}

impl<R: ComRing> Ergonomic<R> {
    pub fn new(elem: R) -> Self {
        Self { elem }
    }

    pub fn pow(&self, n: usize) -> Self {
        Self {
            elem: self.elem.pow(&Natural::from(n)),
        }
    }

    pub fn to_elem(self) -> R {
        self.elem
    }

    pub fn elem(&self) -> R {
        self.elem.clone()
    }
}

//val + val
impl<R: ComRing> std::ops::Add for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::add(self.elem, other.elem),
        }
    }
}

//ref + ref
impl<R: ComRing> std::ops::Add for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::add_refs(&self.elem, &other.elem),
        }
    }
}

//val + ref
impl<R: ComRing> std::ops::Add<&Ergonomic<R>> for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::add_ref(self.elem, &other.elem),
        }
    }
}

//ref + val
impl<R: ComRing> std::ops::Add<Ergonomic<R>> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::add_ref(other.elem, &self.elem),
        }
    }
}

//val - val
impl<R: ComRing> std::ops::Sub for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn sub(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::add(self.elem, other.elem.neg()),
        }
    }
}

//ref - ref
impl<R: ComRing> std::ops::Sub for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn sub(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::add_refs(&self.elem, &other.elem.neg_ref()),
        }
    }
}

//val - ref
impl<R: ComRing> std::ops::Sub<&Ergonomic<R>> for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn sub(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::add_ref(self.elem, &other.elem.neg_ref()),
        }
    }
}

//ref - val
impl<R: ComRing> std::ops::Sub<Ergonomic<R>> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn sub(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::add_ref(other.elem.neg(), &self.elem),
        }
    }
}

//-val
impl<R: ComRing> std::ops::Neg for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn neg(self) -> Self::Output {
        Self::Output {
            elem: self.elem.neg(),
        }
    }
}

//-ref
impl<R: ComRing> std::ops::Neg for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn neg(self) -> Self::Output {
        Self::Output {
            elem: self.elem.neg_ref(),
        }
    }
}

//val * val
impl<R: ComRing> std::ops::Mul for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::mul(self.elem, other.elem),
        }
    }
}

//ref * ref
impl<R: ComRing> std::ops::Mul for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::mul_refs(&self.elem, &other.elem),
        }
    }
}

//val * ref
impl<R: ComRing> std::ops::Mul<&Ergonomic<R>> for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::mul_ref(self.elem, &other.elem),
        }
    }
}

//ref * val
impl<R: ComRing> std::ops::Mul<Ergonomic<R>> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            elem: R::mul_ref(other.elem, &self.elem),
        }
    }
}

//val + i32
impl<R: ComRing> std::ops::Add<i32> for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: i32) -> Self::Output {
        self + Ergonomic::new(R::from_int(&Integer::from(other)))
    }
}

//ref + i32
impl<R: ComRing> std::ops::Add<i32> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: i32) -> Self::Output {
        self + Ergonomic::new(R::from_int(&Integer::from(other)))
    }
}

//i32 + val
impl<R: ComRing> std::ops::Add<Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn add(self, other: Ergonomic<R>) -> Self::Output {
        Ergonomic::new(R::from_int(&Integer::from(self))) + other
    }
}

//i32 + ref
impl<R: ComRing> std::ops::Add<&Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn add(self, other: &Ergonomic<R>) -> Self::Output {
        Ergonomic::new(R::from_int(&Integer::from(self))) + other
    }
}

//val - i32
impl<R: ComRing> std::ops::Sub<i32> for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn sub(self, other: i32) -> Self::Output {
        self - Ergonomic::new(R::from_int(&Integer::from(other)))
    }
}

//ref - i32
impl<R: ComRing> std::ops::Sub<i32> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn sub(self, other: i32) -> Self::Output {
        self - Ergonomic::new(R::from_int(&Integer::from(other)))
    }
}

//i32 - val
impl<R: ComRing> std::ops::Sub<Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn sub(self, other: Ergonomic<R>) -> Self::Output {
        Ergonomic::new(R::from_int(&Integer::from(self))) - other
    }
}

//i32 - ref
impl<R: ComRing> std::ops::Sub<&Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn sub(self, other: &Ergonomic<R>) -> Self::Output {
        Ergonomic::new(R::from_int(&Integer::from(self))) - other
    }
}

//val * i32
impl<R: ComRing> std::ops::Mul<i32> for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: i32) -> Self::Output {
        self * Ergonomic::new(R::from_int(&Integer::from(other)))
    }
}

//ref * i32
impl<R: ComRing> std::ops::Mul<i32> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: i32) -> Self::Output {
        self * Ergonomic::new(R::from_int(&Integer::from(other)))
    }
}

//i32 * val
impl<R: ComRing> std::ops::Mul<Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn mul(self, other: Ergonomic<R>) -> Self::Output {
        Ergonomic::new(R::from_int(&Integer::from(self))) * other
    }
}

//i32 * ref
impl<R: ComRing> std::ops::Mul<&Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn mul(self, other: &Ergonomic<R>) -> Self::Output {
        Ergonomic::new(R::from_int(&Integer::from(self))) * other
    }
}
