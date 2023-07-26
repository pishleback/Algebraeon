use malachite_nz::{natural::Natural, integer::Integer};

use super::ring::*;


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
