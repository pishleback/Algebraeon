use malachite_nz::{integer::Integer, natural::Natural};

use super::ring::*;

#[derive(Debug, Clone)]
pub struct Ergonomic<'a, R: ComRing> {
    ring: &'a R,
    elem: R::ElemT,
}

impl<'a, R: ComRing> PartialEq for Ergonomic<'a, R> {
    fn eq(&self, other: &Self) -> bool {
        if self.ring == other.ring {
            let ans = self.ring.equal(&self.elem, &other.elem);
            debug_assert_eq!(ans, other.ring.equal(&self.elem, &other.elem));
            ans
        } else {
            false
        }
    }
}

impl<'a, R: ComRing> Eq for Ergonomic<'a, R> {}

impl<'a, R: ComRing> Ergonomic<'a, R> {
    pub fn new(ring: &'a R, elem: R::ElemT) -> Self {
        Self { ring, elem }
    }

    pub fn pow(&self, n: usize) -> Self {
        Self {
            ring: self.ring,
            elem: self.ring.nat_pow(&self.elem, &Natural::from(n)),
        }
    }

    pub fn to_elem(self) -> R::ElemT {
        self.elem
    }

    pub fn elem(&self) -> R::ElemT {
        self.elem.clone()
    }
}

//val + val
impl<'a, R: ComRing> std::ops::Add for Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn add(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add(self.elem, other.elem),
        }
    }
}

//ref + ref
impl<'a, R: ComRing> std::ops::Add for &Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn add(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add_refs(&self.elem, &other.elem),
        }
    }
}

//val + ref
impl<'a, R: ComRing> std::ops::Add<&Ergonomic<'a, R>> for Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn add(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add_ref(self.elem, &other.elem),
        }
    }
}

//ref + val
impl<'a, R: ComRing> std::ops::Add<Ergonomic<'a, R>> for &Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn add(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add_ref(other.elem, &self.elem),
        }
    }
}

//val - val
impl<'a, R: ComRing> std::ops::Sub for Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn sub(self, other: Ergonomic<'a, R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add(self.elem, self.ring.neg(other.elem)),
        }
    }
}

//ref - ref
impl<'a, R: ComRing> std::ops::Sub for &Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn sub(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self
                .ring
                .add_refs(&self.elem, &self.ring.neg_ref(&other.elem)),
        }
    }
}

//val - ref
impl<'a, R: ComRing> std::ops::Sub<&Ergonomic<'a, R>> for Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn sub(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self
                .ring
                .add_ref(self.elem, &self.ring.neg_ref(&other.elem)),
        }
    }
}

//ref - val
impl<'a, R: ComRing> std::ops::Sub<Ergonomic<'a, R>> for &Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn sub(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add_ref(self.ring.neg(other.elem), &self.elem),
        }
    }
}

//-val
impl<'a, R: ComRing> std::ops::Neg for Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn neg(self) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.neg(self.elem),
        }
    }
}

//-ref
impl<'a, R: ComRing> std::ops::Neg for &Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn neg(self) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.neg_ref(&self.elem),
        }
    }
}

//val * val
impl<'a, R: ComRing> std::ops::Mul for Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn mul(self, other: Ergonomic<'a, R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.mul(self.elem, other.elem),
        }
    }
}

//ref * ref
impl<'a, R: ComRing> std::ops::Mul for &Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn mul(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.mul_refs(&self.elem, &other.elem),
        }
    }
}

//val * ref
impl<'a, R: ComRing> std::ops::Mul<&Ergonomic<'a, R>> for Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn mul(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.mul_ref(self.elem, &other.elem),
        }
    }
}

//ref * val
impl<'a, R: ComRing> std::ops::Mul<Ergonomic<'a, R>> for &Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn mul(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.mul_ref(other.elem, &self.elem),
        }
    }
}

//val + i32
impl<'a, R: ComRing> std::ops::Add<i32> for Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn add(self, other: i32) -> Self::Output {
        &self + Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//ref + i32
impl<'a, R: ComRing> std::ops::Add<i32> for &Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn add(self, other: i32) -> Self::Output {
        self + Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//i32 + val
impl<'a, R: ComRing> std::ops::Add<Ergonomic<'a, R>> for i32 {
    type Output = Ergonomic<'a, R>;

    fn add(self, other: Ergonomic<'a, R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) + other
    }
}

//i32 + ref
impl<'a, R: ComRing> std::ops::Add<&Ergonomic<'a, R>> for i32 {
    type Output = Ergonomic<'a, R>;

    fn add(self, other: &Ergonomic<'a, R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) + other
    }
}

//val - i32
impl<'a, R: ComRing> std::ops::Sub<i32> for Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn sub(self, other: i32) -> Self::Output {
        &self - Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//ref - i32
impl<'a, R: ComRing> std::ops::Sub<i32> for &Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn sub(self, other: i32) -> Self::Output {
        self - Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//i32 - val
impl<'a, R: ComRing> std::ops::Sub<Ergonomic<'a, R>> for i32 {
    type Output = Ergonomic<'a, R>;

    fn sub(self, other: Ergonomic<'a, R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) - other
    }
}

//i32 - ref
impl<'a, R: ComRing> std::ops::Sub<&Ergonomic<'a, R>> for i32 {
    type Output = Ergonomic<'a, R>;

    fn sub(self, other: &Ergonomic<'a, R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) - other
    }
}

//val * i32
impl<'a, R: ComRing> std::ops::Mul<i32> for Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn mul(self, other: i32) -> Self::Output {
        &self * Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//ref * i32
impl<'a, R: ComRing> std::ops::Mul<i32> for &Ergonomic<'a, R> {
    type Output = Ergonomic<'a, R>;

    fn mul(self, other: i32) -> Self::Output {
        self * Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//i32 * val
impl<'a, R: ComRing> std::ops::Mul<Ergonomic<'a, R>> for i32 {
    type Output = Ergonomic<'a, R>;

    fn mul(self, other: Ergonomic<'a, R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) * other
    }
}

//i32 * ref
impl<'a, R: ComRing> std::ops::Mul<&Ergonomic<'a, R>> for i32 {
    type Output = Ergonomic<'a, R>;

    fn mul(self, other: &Ergonomic<'a, R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) * other
    }
}
