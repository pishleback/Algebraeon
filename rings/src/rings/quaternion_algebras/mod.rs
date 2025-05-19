use crate::structure::*;
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::*;

#[derive(Debug, Clone)]
pub struct QuaternionAlgebraStructure<Field: FieldSignature> {
    base: Field,
    is_char_2: bool,
    a: Field::Set,
    b: Field::Set,
}

impl<Field: FieldSignature + CharacteristicSignature> QuaternionAlgebraStructure<Field> {
    pub fn new(base: Field, a: Field::Set, b: Field::Set) -> Self {
        let is_char_2 = base.characteristic() == Natural::TWO;
        Self {
            base,
            is_char_2,
            a,
            b,
        }
    }
}

#[derive(Debug, Clone)]
pub struct QuaternionAlgebraElement<Field: FieldSignature> {
    coeffs: [Field::Set; 4],
}

impl<Field: FieldSignature> PartialEq for QuaternionAlgebraStructure<Field> {
    fn eq(&self, other: &Self) -> bool {
        self.base == other.base
            && self.base.equal(&self.a, &other.a)
            && self.base.equal(&self.b, &other.b)
    }
}

impl<Field: FieldSignature> Eq for QuaternionAlgebraStructure<Field> {}

impl<Field: FieldSignature> EqSignature for QuaternionAlgebraStructure<Field> {
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.equal_elements(a, b)
    }
}

impl<Field: FieldSignature> Signature for QuaternionAlgebraStructure<Field> {}

impl<Field: FieldSignature> SetSignature for QuaternionAlgebraStructure<Field> {
    type Set = QuaternionAlgebraElement<Field>;

    fn is_element(&self, _x: &Self::Set) -> bool {
        true
    }
}

impl<Field: FieldSignature> SemiRingSignature for QuaternionAlgebraStructure<Field> {
    fn zero(&self) -> Self::Set {
        QuaternionAlgebraElement {
            coeffs: std::array::from_fn(|_| self.base.zero()),
        }
    }

    fn one(&self) -> Self::Set {
        QuaternionAlgebraElement {
            coeffs: [
                self.base.one(),
                self.base.zero(),
                self.base.zero(),
                self.base.zero(),
            ],
        }
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let mut result = std::array::from_fn(|_| self.base.zero());
        for i in 0..4 {
            result[i] = self.base.add(&a.coeffs[i], &b.coeffs[i]);
        }
        QuaternionAlgebraElement { coeffs: result }
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let (x0, x1, x2, x3) = (&a.coeffs[0], &a.coeffs[1], &a.coeffs[2], &a.coeffs[3]);
        let (y0, y1, y2, y3) = (&b.coeffs[0], &b.coeffs[1], &b.coeffs[2], &b.coeffs[3]);

        let a_param = &self.a;
        let b_param = &self.b;
        let base = &self.base;
        let ab = base.mul(a_param, b_param);

        if self.is_char_2 {
            // Quaternion multiplication in characteristic 2.
            //
            //   i^2 + i = a
            //   j^2 = b
            //   k = ij = j(i + 1)
            unimplemented!("Quaternion multiplication for char=2 is not implemented yet");
        } else {
            // Quaternion multiplication in characteristic ≠ 2.
            //
            //   i^2 = a
            //   j^2 = b
            //   ij = k = -ji

            let z0 = base.sub(
                &base.add(
                    &base.add(&base.mul(x0, y0), &base.mul(&base.mul(x1, y1), a_param)),
                    &base.mul(&base.mul(x2, y2), b_param),
                ),
                &base.mul(&base.mul(x3, y3), &ab),
            );
            let z1 = base.sub(
                &base.add(
                    &base.add(&base.mul(x0, y1), &base.mul(x1, y0)),
                    &base.mul(&base.mul(x2, y3), b_param),
                ),
                &base.mul(&base.mul(x3, y2), b_param),
            );
            let z2 = base.add(
                &base.sub(
                    &base.add(&base.mul(x0, y2), &base.mul(x2, y0)),
                    &base.mul(&base.mul(x1, y3), a_param),
                ),
                &base.mul(&base.mul(x3, y1), a_param),
            );
            let z3 = base.add(
                &base.sub(
                    &base.add(&base.mul(x0, y3), &base.mul(x3, y0)),
                    &base.mul(x2, y1),
                ),
                &base.mul(x1, y2),
            );

            QuaternionAlgebraElement {
                coeffs: [z0, z1, z2, z3],
            }
        }
    }
}

impl<Field: FieldSignature> RingSignature for QuaternionAlgebraStructure<Field> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        QuaternionAlgebraElement {
            coeffs: std::array::from_fn(|i| self.base.neg(&a.coeffs[i])),
        }
    }
}

impl<Field: FieldSignature + CharacteristicSignature> CharacteristicSignature
    for QuaternionAlgebraStructure<Field>
{
    fn characteristic(&self) -> Natural {
        self.base.characteristic()
    }
}

impl<Field: CharZeroFieldSignature> CharZeroRingSignature for QuaternionAlgebraStructure<Field> {
    fn try_to_int(&self, x: &Self::Set) -> Option<algebraeon_nzq::Integer> {
        // The element must be of the form [x0, 0, 0, 0]
        if self.base.is_zero(&x.coeffs[1])
            && self.base.is_zero(&x.coeffs[2])
            && self.base.is_zero(&x.coeffs[3])
        {
            self.base.try_to_int(&x.coeffs[0])
        } else {
            None
        }
    }
}

impl<Field: FieldSignature> QuaternionAlgebraStructure<Field> {
    pub fn i(&self) -> QuaternionAlgebraElement<Field> {
        QuaternionAlgebraElement {
            coeffs: [
                self.base.zero(),
                self.base.one(),
                self.base.zero(),
                self.base.zero(),
            ],
        }
    }

    pub fn j(&self) -> QuaternionAlgebraElement<Field> {
        QuaternionAlgebraElement {
            coeffs: [
                self.base.zero(),
                self.base.zero(),
                self.base.one(),
                self.base.zero(),
            ],
        }
    }

    pub fn k(&self) -> QuaternionAlgebraElement<Field> {
        QuaternionAlgebraElement {
            coeffs: [
                self.base.zero(),
                self.base.zero(),
                self.base.zero(),
                self.base.one(),
            ],
        }
    }

    pub fn conjugate(
        &self,
        a: &QuaternionAlgebraElement<Field>,
    ) -> QuaternionAlgebraElement<Field> {
        let base = &self.base;
        if self.is_char_2 {
            // https://jvoight.github.io/quat-book.pdf paragraph 6.2.6.
            QuaternionAlgebraElement {
                coeffs: [
                    base.add(&a.coeffs[0], &a.coeffs[1]),
                    a.coeffs[1].clone(),
                    a.coeffs[2].clone(),
                    a.coeffs[3].clone(),
                ],
            }
        } else {
            QuaternionAlgebraElement {
                coeffs: [
                    a.coeffs[0].clone(),
                    base.neg(&a.coeffs[1]),
                    base.neg(&a.coeffs[2]),
                    base.neg(&a.coeffs[3]),
                ],
            }
        }
    }

    pub fn reduced_trace(&self, a: &QuaternionAlgebraElement<Field>) -> Field::Set {
        let base = &self.base;
        let a0 = &a.coeffs[0];
        let a1 = &a.coeffs[1];
        if self.is_char_2 {
            // https://jvoight.github.io/quat-book.pdf paragraph 6.2.6.
            a1.clone()
        } else {
            base.add(&a0, &a0) // 2 * a0
        }
    }

    pub fn reduced_norm(&self, a: &QuaternionAlgebraElement<Field>) -> Field::Set {
        let base = &self.base;
        let a0 = &a.coeffs[0];
        let a1 = &a.coeffs[1];
        let a2 = &a.coeffs[2];
        let a3 = &a.coeffs[3];
        let a_param = &self.a;
        let b_param = &self.b;
        let ab = base.mul(a_param, b_param);

        if self.is_char_2 {
            // https://jvoight.github.io/quat-book.pdf equation 6.2.7.
            // t^2 + t·x + a·x^2 + b·y^2 + b·y·z + ab·z^2
            let t2 = base.mul(a0, a0);
            let tx = base.mul(a0, a1);
            let ax2 = base.mul(a_param, &base.mul(a1, a1));
            let by2 = base.mul(b_param, &base.mul(a2, a2));
            let byz = base.mul(b_param, &base.mul(a2, a3));
            let abz2 = base.mul(&ab, &base.mul(a3, a3));

            base.add(
                &base.add(&base.add(&base.add(&base.add(&t2, &tx), &ax2), &by2), &byz),
                &abz2,
            )
        } else {
            let term0 = base.mul(a0, a0);
            let term1 = base.mul(a_param, &base.mul(a1, a1));
            let term2 = base.mul(b_param, &base.mul(a2, a2));
            let term3 = base.mul(&ab, &base.mul(a3, a3));

            base.sub(&base.sub(&base.add(&term0, &term3), &term1), &term2)
        }
    }

    pub fn equal_elements(
        &self,
        a: &QuaternionAlgebraElement<Field>,
        b: &QuaternionAlgebraElement<Field>,
    ) -> bool {
        (0..4).all(|i| self.base.equal(&a.coeffs[i], &b.coeffs[i]))
    }
}

#[cfg(test)]
mod tests {
    use algebraeon_nzq::Rational;
    use algebraeon_nzq::RationalCanonicalStructure;

    use super::*;

    #[test]
    fn test_add_commutativity() {
        // Hamilton quaternion algebra: H = (-1, -1 / QQ)
        let field = RationalCanonicalStructure {};
        let h = QuaternionAlgebraStructure::new(field, -Rational::ONE, -Rational::ONE);

        let i = h.i();
        let j = h.j();
        let i_plus_j = h.add(&i, &j);
        let j_plus_i = h.add(&j, &i);
        let i_times_j = h.mul(&i, &j);
        let j_times_i = h.mul(&j, &i);

        assert!(h.equal_elements(&i_plus_j, &j_plus_i));
        assert!(h.equal_elements(&i_times_j, &h.neg(&j_times_i)));
    }
}
