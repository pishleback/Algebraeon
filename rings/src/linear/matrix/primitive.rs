use super::*;

impl<RS: GreatestCommonDivisorStructure> MatrixStructure<RS> {
    pub fn factor_primitive(&self, mut mat: Matrix<RS::Set>) -> Option<(RS::Set, Matrix<RS::Set>)> {
        let entries = mat.entries_list();
        let g = self.ring().gcd_list(entries);
        if self.ring().is_zero(&g) {
            None
        } else {
            for r in 0..mat.rows() {
                for c in 0..mat.cols() {
                    *mat.at_mut(r, c).unwrap() =
                        self.ring().div(mat.at(r, c).unwrap(), &g).unwrap();
                }
            }
            Some((g, mat))
        }
    }

    pub fn primitive_part(&self, mat: Matrix<RS::Set>) -> Option<Matrix<RS::Set>> {
        match self.factor_primitive(mat) {
            Some((_unit, prim)) => Some(prim),
            None => None,
        }
    }
}

pub fn factor_primitive_fof<
    Ring: GreatestCommonDivisorStructure,
    Field: FieldStructure,
    Fof: FieldOfFractionsInclusion<Ring, Field>,
>(
    fof_inclusion: &Fof,
    mat: &Matrix<Field::Set>,
) -> (Field::Set, Matrix<Ring::Set>) {
    let ring = fof_inclusion.domain();
    let field = fof_inclusion.range();
    let mat_ring = MatrixStructure::new(ring.clone());

    let div = ring.lcm_list(
        mat.entries_list()
            .into_iter()
            .map(|c| fof_inclusion.denominator(&c))
            .collect(),
    );

    let (mul, prim) = mat_ring
        .factor_primitive(mat.apply_map(|c| {
            fof_inclusion
                .try_preimage(&field.mul(&fof_inclusion.image(&div), c))
                .unwrap()
        }))
        .unwrap();

    (
        field
            .div(&fof_inclusion.image(&mul), &fof_inclusion.image(&div))
            .unwrap(),
        prim,
    )
}

impl<R: MetaType> Matrix<R>
where
    R::Structure: GreatestCommonDivisorStructure,
{
    pub fn factor_primitive(self) -> Option<(R, Matrix<R>)> {
        Self::structure().factor_primitive(self)
    }

    pub fn primitive_part(self) -> Option<Matrix<R>> {
        Self::structure().primitive_part(self)
    }
}

impl<Field: MetaType> Matrix<Field>
where
    Field::Structure: FieldStructure,
    PrincipalSubringInclusion<Field::Structure>:
        FieldOfFractionsInclusion<IntegerCanonicalStructure, Field::Structure>,
{
    pub fn factor_primitive_fof(&self) -> (Field, Matrix<Integer>) {
        factor_primitive_fof(&PrincipalSubringInclusion::new(Field::structure()), self)
    }
}
