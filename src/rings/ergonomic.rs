use malachite_nz::{integer::Integer, natural::Natural};

use super::ring::*;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ergonomic<R: ComRing> {
    ring: R,
    elem: R::ElemT,
}

impl<R: ComRing> Ergonomic<R> {
    pub fn new(ring: R, elem: R::ElemT) -> Self {
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
impl<R: ComRing> std::ops::Add for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add(self.elem, other.elem),
        }
    }
}

//ref + ref
impl<R: ComRing> std::ops::Add for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add_refs(&self.elem, &other.elem),
        }
    }
}

//val + ref
impl<R: ComRing> std::ops::Add<&Ergonomic<R>> for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add_ref(self.elem, &other.elem),
        }
    }
}

//ref + val
impl<R: ComRing> std::ops::Add<Ergonomic<R>> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add_ref(other.elem, &self.elem),
        }
    }
}

//val - val
impl<R: ComRing> std::ops::Sub for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn sub(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add(self.elem, self.ring.neg(other.elem)),
        }
    }
}

//ref - ref
impl<R: ComRing> std::ops::Sub for &Ergonomic<R> {
    type Output = Ergonomic<R>;

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
impl<R: ComRing> std::ops::Sub<&Ergonomic<R>> for Ergonomic<R> {
    type Output = Ergonomic<R>;

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
impl<R: ComRing> std::ops::Sub<Ergonomic<R>> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn sub(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.add_ref(self.ring.neg(other.elem), &self.elem),
        }
    }
}

//-val
impl<R: ComRing> std::ops::Neg for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn neg(self) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.neg(self.elem),
        }
    }
}

//-ref
impl<R: ComRing> std::ops::Neg for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn neg(self) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.neg_ref(&self.elem),
        }
    }
}

//val * val
impl<R: ComRing> std::ops::Mul for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.mul(self.elem, other.elem),
        }
    }
}

//ref * ref
impl<R: ComRing> std::ops::Mul for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.mul_refs(&self.elem, &other.elem),
        }
    }
}

//val * ref
impl<R: ComRing> std::ops::Mul<&Ergonomic<R>> for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: &Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.mul_ref(self.elem, &other.elem),
        }
    }
}

//ref * val
impl<R: ComRing> std::ops::Mul<Ergonomic<R>> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: Ergonomic<R>) -> Self::Output {
        Self::Output {
            ring: self.ring,
            elem: self.ring.mul_ref(other.elem, &self.elem),
        }
    }
}

//val + i32
impl<R: ComRing> std::ops::Add<i32> for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: i32) -> Self::Output {
        self + Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//ref + i32
impl<R: ComRing> std::ops::Add<i32> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn add(self, other: i32) -> Self::Output {
        self + Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//i32 + val
impl<R: ComRing> std::ops::Add<Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn add(self, other: Ergonomic<R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) + other
    }
}

//i32 + ref
impl<R: ComRing> std::ops::Add<&Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn add(self, other: &Ergonomic<R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) + other
    }
}

//val - i32
impl<R: ComRing> std::ops::Sub<i32> for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn sub(self, other: i32) -> Self::Output {
        self - Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//ref - i32
impl<R: ComRing> std::ops::Sub<i32> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn sub(self, other: i32) -> Self::Output {
        self - Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//i32 - val
impl<R: ComRing> std::ops::Sub<Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn sub(self, other: Ergonomic<R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) - other
    }
}

//i32 - ref
impl<R: ComRing> std::ops::Sub<&Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn sub(self, other: &Ergonomic<R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) - other
    }
}

//val * i32
impl<R: ComRing> std::ops::Mul<i32> for Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: i32) -> Self::Output {
        self * Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//ref * i32
impl<R: ComRing> std::ops::Mul<i32> for &Ergonomic<R> {
    type Output = Ergonomic<R>;

    fn mul(self, other: i32) -> Self::Output {
        self * Ergonomic::new(self.ring, self.ring.from_int(&Integer::from(other)))
    }
}

//i32 * val
impl<R: ComRing> std::ops::Mul<Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn mul(self, other: Ergonomic<R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) * other
    }
}

//i32 * ref
impl<R: ComRing> std::ops::Mul<&Ergonomic<R>> for i32 {
    type Output = Ergonomic<R>;

    fn mul(self, other: &Ergonomic<R>) -> Self::Output {
        Ergonomic::new(other.ring, other.ring.from_int(&Integer::from(self))) * other
    }
}
