use crate::structure::*;
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::*;
use std::{borrow::Cow, fmt};

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

impl<Field: FieldSignature + ToStringSignature + fmt::Display> fmt::Display
    for QuaternionAlgebraStructure<Field>
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Quaternion Algebra ({}, {}) over base field {}",
            self.base.to_string(&self.a),
            self.base.to_string(&self.b),
            self.base
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, CanonicalStructure)]
#[canonical_structure(eq)]
pub enum QuaternionAlgebraBasis {
    R,
    I,
    J,
    K,
}

impl CountableSetSignature for QuaternionAlgebraBasisCanonicalStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Set> {
        vec![
            QuaternionAlgebraBasis::R,
            QuaternionAlgebraBasis::I,
            QuaternionAlgebraBasis::J,
            QuaternionAlgebraBasis::K,
        ]
        .into_iter()
    }
}

impl FiniteSetSignature for QuaternionAlgebraBasisCanonicalStructure {
    fn size(&self) -> usize {
        4
    }
}

#[derive(Debug, Clone)]
pub struct QuaternionAlgebraElement<Field: FieldSignature> {
    // represent x + yi + zj + wk
    x: Field::Set,
    y: Field::Set,
    z: Field::Set,
    w: Field::Set,
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
        self.base.equal(&a.x, &b.x)
            && self.base.equal(&a.y, &b.y)
            && self.base.equal(&a.z, &b.z)
            && self.base.equal(&a.w, &b.w)
    }
}

impl<Field: FieldSignature> Signature for QuaternionAlgebraStructure<Field> {}

impl<Field: FieldSignature> SetSignature for QuaternionAlgebraStructure<Field> {
    type Set = QuaternionAlgebraElement<Field>;

    fn is_element(&self, _x: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl<Field: FieldSignature> AdditiveMonoidSignature for QuaternionAlgebraStructure<Field> {
    fn zero(&self) -> Self::Set {
        QuaternionAlgebraElement {
            x: self.base.zero(),
            y: self.base.zero(),
            z: self.base.zero(),
            w: self.base.zero(),
        }
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        QuaternionAlgebraElement {
            x: self.base.add(&a.x, &b.x),
            y: self.base.add(&a.y, &b.y),
            z: self.base.add(&a.z, &b.z),
            w: self.base.add(&a.w, &b.w),
        }
    }
}

impl<Field: FieldSignature> AdditiveGroupSignature for QuaternionAlgebraStructure<Field> {
    fn neg(&self, a: &Self::Set) -> Self::Set {
        QuaternionAlgebraElement {
            x: self.base.neg(&a.x),
            y: self.base.neg(&a.y),
            z: self.base.neg(&a.z),
            w: self.base.neg(&a.w),
        }
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        QuaternionAlgebraElement {
            x: self.base.sub(&a.x, &b.x),
            y: self.base.sub(&a.y, &b.y),
            z: self.base.sub(&a.z, &b.z),
            w: self.base.sub(&a.w, &b.w),
        }
    }
}

impl<Field: FieldSignature> SemiRingSignature for QuaternionAlgebraStructure<Field> {
    fn one(&self) -> Self::Set {
        QuaternionAlgebraElement {
            x: self.base.one(),
            y: self.base.zero(),
            z: self.base.zero(),
            w: self.base.zero(),
        }
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
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
                    &base.add(
                        &base.mul(&a.x, &b.x),
                        &base.mul(&base.mul(&a.y, &b.y), a_param),
                    ),
                    &base.mul(&base.mul(&a.z, &b.z), b_param),
                ),
                &base.mul(&base.mul(&a.w, &b.w), &ab),
            );
            let z1 = base.sub(
                &base.add(
                    &base.add(&base.mul(&a.x, &b.y), &base.mul(&a.y, &b.x)),
                    &base.mul(&base.mul(&a.z, &b.w), b_param),
                ),
                &base.mul(&base.mul(&a.w, &b.z), b_param),
            );
            let z2 = base.add(
                &base.sub(
                    &base.add(&base.mul(&a.x, &b.z), &base.mul(&a.z, &b.x)),
                    &base.mul(&base.mul(&a.y, &b.w), a_param),
                ),
                &base.mul(&base.mul(&a.w, &b.y), a_param),
            );
            let z3 = base.add(
                &base.sub(
                    &base.add(&base.mul(&a.x, &b.w), &base.mul(&a.w, &b.x)),
                    &base.mul(&a.z, &b.y),
                ),
                &base.mul(&a.y, &b.z),
            );

            QuaternionAlgebraElement {
                x: z0,
                y: z1,
                z: z2,
                w: z3,
            }
        }
    }
}

impl<Field: FieldSignature + SemiRingUnitsSignature> SemiRingUnitsSignature
    for QuaternionAlgebraStructure<Field>
{
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, RingDivisionError> {
        let n_inv = self.base.inv(&self.reduced_norm(a))?;
        Ok(self.scalar_mul(&n_inv, &self.conjugate(a)))
    }
}

impl<Field: FieldSignature> RingSignature for QuaternionAlgebraStructure<Field> {}

impl<Field: FieldSignature> SemiModuleSignature<Field> for QuaternionAlgebraStructure<Field> {
    fn ring(&self) -> &Field {
        &self.base
    }

    fn scalar_mul(&self, x: &<Field>::Set, a: &Self::Set) -> Self::Set {
        let base = &self.base;
        QuaternionAlgebraElement {
            x: base.mul(x, &a.x),
            y: base.mul(x, &a.y),
            z: base.mul(x, &a.z),
            w: base.mul(x, &a.w),
        }
    }
}

impl<Field: FieldSignature> AlgebraSignature<Field> for QuaternionAlgebraStructure<Field> {}

impl<Field: FieldSignature> FreeModuleSignature<Field> for QuaternionAlgebraStructure<Field> {
    type Basis = QuaternionAlgebraBasisCanonicalStructure;

    fn basis_set(&self) -> impl std::borrow::Borrow<Self::Basis> {
        QuaternionAlgebraBasisCanonicalStructure {}
    }

    fn to_component<'a>(
        &self,
        b: &QuaternionAlgebraBasis,
        v: &'a Self::Set,
    ) -> Cow<'a, Field::Set> {
        Cow::Borrowed(match b {
            QuaternionAlgebraBasis::R => &v.x,
            QuaternionAlgebraBasis::I => &v.y,
            QuaternionAlgebraBasis::J => &v.z,
            QuaternionAlgebraBasis::K => &v.w,
        })
    }

    fn from_component(&self, b: &QuaternionAlgebraBasis, r: &Field::Set) -> Self::Set {
        let mut v = self.zero();
        match b {
            QuaternionAlgebraBasis::R => v.x = r.clone(),
            QuaternionAlgebraBasis::I => v.y = r.clone(),
            QuaternionAlgebraBasis::J => v.z = r.clone(),
            QuaternionAlgebraBasis::K => v.w = r.clone(),
        }
        v
    }
}

impl<Field: FieldSignature> FinitelyFreeModuleSignature<Field>
    for QuaternionAlgebraStructure<Field>
{
}

impl<Field: FieldSignature + CharacteristicSignature> CharacteristicSignature
    for QuaternionAlgebraStructure<Field>
{
    fn characteristic(&self) -> Natural {
        self.base.characteristic()
    }
}

impl<Field: CharZeroFieldSignature> CharZeroRingSignature for QuaternionAlgebraStructure<Field> {
    fn try_to_int(&self, a: &Self::Set) -> Option<algebraeon_nzq::Integer> {
        // The element must be of the form [a.x, 0, 0, 0]
        if self.base.is_zero(&a.y) && self.base.is_zero(&a.z) && self.base.is_zero(&a.w) {
            self.base.try_to_int(&a.x)
        } else {
            None
        }
    }
}

impl<Field: FieldSignature> QuaternionAlgebraStructure<Field> {
    pub fn i(&self) -> QuaternionAlgebraElement<Field> {
        QuaternionAlgebraElement {
            x: self.base.zero(),
            y: self.base.one(),
            z: self.base.zero(),
            w: self.base.zero(),
        }
    }

    pub fn j(&self) -> QuaternionAlgebraElement<Field> {
        QuaternionAlgebraElement {
            x: self.base.zero(),
            y: self.base.zero(),
            z: self.base.one(),
            w: self.base.zero(),
        }
    }

    pub fn k(&self) -> QuaternionAlgebraElement<Field> {
        QuaternionAlgebraElement {
            x: self.base.zero(),
            y: self.base.zero(),
            z: self.base.zero(),
            w: self.base.one(),
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
                x: base.add(&a.x, &a.y),
                y: a.y.clone(),
                z: a.z.clone(),
                w: a.w.clone(),
            }
        } else {
            QuaternionAlgebraElement {
                x: a.x.clone(),
                y: base.neg(&a.y),
                z: base.neg(&a.z),
                w: base.neg(&a.w),
            }
        }
    }

    pub fn reduced_trace(&self, a: &QuaternionAlgebraElement<Field>) -> Field::Set {
        if self.is_char_2 {
            // https://jvoight.github.io/quat-book.pdf paragraph 6.2.6.
            a.y.clone()
        } else {
            self.base.add(&a.x, &a.x) // 2 * &a.x
        }
    }

    pub fn reduced_norm(&self, a: &QuaternionAlgebraElement<Field>) -> Field::Set {
        let base = &self.base;
        let a_param = &self.a;
        let b_param = &self.b;
        let ab = base.mul(a_param, b_param);

        if self.is_char_2 {
            // https://jvoight.github.io/quat-book.pdf equation 6.2.7.
            // t^2 + t·x + a·x^2 + b·y^2 + b·y·z + ab·z^2
            let t2 = base.mul(&a.x, &a.x);
            let tx = base.mul(&a.x, &a.y);
            let ax2 = base.mul(a_param, &base.mul(&a.y, &a.y));
            let by2 = base.mul(b_param, &base.mul(&a.z, &a.z));
            let byz = base.mul(b_param, &base.mul(&a.z, &a.w));
            let abz2 = base.mul(&ab, &base.mul(&a.w, &a.w));
            base.add(
                &base.add(&base.add(&base.add(&base.add(&t2, &tx), &ax2), &by2), &byz),
                &abz2,
            )
        } else {
            let term0 = base.mul(&a.x, &a.x);
            let term1 = base.mul(a_param, &base.mul(&a.y, &a.y));
            let term2 = base.mul(b_param, &base.mul(&a.z, &a.z));
            let term3 = base.mul(&ab, &base.mul(&a.w, &a.w));

            base.sub(&base.sub(&base.add(&term0, &term3), &term1), &term2)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use algebraeon_nzq::Rational;

    #[test]
    fn test_add_and_mul() {
        // Hamilton quaternion algebra: H = (-1, -1 / QQ)
        let h =
            QuaternionAlgebraStructure::new(Rational::structure(), -Rational::ONE, -Rational::ONE);

        let i = h.i();
        let j = h.j();
        let i_plus_j = h.add(&i, &j);
        let j_plus_i = h.add(&j, &i);
        let i_times_j = h.mul(&i, &j);
        let j_times_i = h.mul(&j, &i);

        assert!(h.equal(&i_plus_j, &j_plus_i));
        assert!(h.equal(&i_times_j, &h.neg(&j_times_i)));
    }

    #[test]
    fn test_reduced_norm_from_sage_example() {
        let h = QuaternionAlgebraStructure::new(
            Rational::structure(),
            Rational::from(-5i32),
            Rational::from(-2i32),
        );

        assert_eq!(h.reduced_norm(&h.i()), Rational::from(5i32));
        assert_eq!(h.reduced_norm(&h.j()), Rational::from(2i32));

        let a = QuaternionAlgebraElement {
            x: Rational::from_integers(1, 3),
            y: Rational::from_integers(1, 5),
            z: Rational::from_integers(1, 7),
            w: Rational::ONE,
        };

        let expected = Rational::from_integers(22826, 2205);
        assert_eq!(h.reduced_norm(&a), expected);
    }
}
