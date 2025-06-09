use crate::{polynomial::*, structure::*};
use algebraeon_nzq::*;
use algebraeon_sets::structure::*;
use itertools::Itertools;

use super::modulo::ModuloCanonicalStructure;

impl<
    FS: FiniteFieldSignature,
    FSB: BorrowedStructure<FS>,
    FSPB: BorrowedStructure<PolynomialStructure<FS, FSB>>,
> FiniteUnitsSignature for FieldExtensionByPolynomialQuotientStructure<FS, FSB, FSPB>
{
    fn all_units(&self) -> Vec<Self::Set> {
        let mut all_base_elements = vec![self.ring().coeff_ring().zero()];
        for unit in self.ring().coeff_ring().all_units() {
            all_base_elements.push(unit);
        }

        let mut all_base_elements_product = (0..self.degree())
            .map(|_| &all_base_elements)
            .multi_cartesian_product();

        // Pop the all-zeros element
        all_base_elements_product.next().unwrap();

        // What remains is the coefficients for all non-zero elements in self
        all_base_elements_product
            .map(|coeffs| {
                self.ring().reduce_poly(Polynomial::from_coeffs(
                    coeffs.into_iter().cloned().collect(),
                ))
            })
            .collect()
    }
}

impl<
    FS: FiniteFieldSignature,
    FSB: BorrowedStructure<FS>,
    FSPB: BorrowedStructure<PolynomialStructure<FS, FSB>>,
> FiniteFieldSignature for FieldExtensionByPolynomialQuotientStructure<FS, FSB, FSPB>
{
    fn characteristic_and_power(&self) -> (Natural, Natural) {
        let (p, t) = self.ring().coeff_ring().characteristic_and_power();
        let d = Natural::from(self.degree());
        (p, d * t)
    }
}

pub fn new_finite_field_extension<FS: FiniteFieldSignature>(
    finite_field: FS,
    poly: Polynomial<FS::Set>,
) -> FieldExtensionByPolynomialQuotientStructure<FS, FS, PolynomialStructure<FS, FS>>
where
    PolynomialStructure<FS, FS>: FactorableSignature<Set = Polynomial<FS::Set>>,
{
    finite_field
        .into_polynomial_ring()
        .into_quotient_field_unchecked(poly)
}

pub(crate) fn f9() -> FieldExtensionByPolynomialQuotientStructure<
    ModuloCanonicalStructure<3>,
    ModuloCanonicalStructure<3>,
    PolynomialStructure<ModuloCanonicalStructure<3>, ModuloCanonicalStructure<3>>,
> {
    use crate::finite_fields::modulo::*;
    new_finite_field_extension::<ModuloCanonicalStructure<3>>(
        Modulo::<3>::structure(),
        Polynomial::from_coeffs(vec![1, 1, 2]),
    )
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_f9_elements() {
        let f9 = f9();

        let (p, t) = f9.characteristic_and_power();
        assert_eq!(p, 3u32.into());
        assert_eq!(t, 2u32.into());

        let mut c = 0;
        for x in f9.all_elements() {
            println!("{:?}", x);
            c += 1;
        }
        assert_eq!(c, 9);
    }
}
