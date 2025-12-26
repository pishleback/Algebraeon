use crate::{
    algebraic_number_field::{
        polynomial_quotient_number_field::AlgebraicNumberFieldPolynomialQuotientStructure,
        quadratic_number_field::{QuadraticNumberFieldElement, QuadraticNumberFieldStructure},
    },
    structure::{MetaFactorableSignature, MetaSemiRing},
};
use algebraeon_nzq::{Integer, Natural, traits::DivMod};
use algebraeon_sets::structure::{BorrowedStructure, Morphism};

#[derive(Debug, Clone)]
pub struct QuadraticNumberFieldIsomorphism<
    AlgebraicNumberFieldPolynomialQuotientStructureBorrowed: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>,
> {
    anf_polyquo: AlgebraicNumberFieldPolynomialQuotientStructureBorrowed,
    anf_quadratic: QuadraticNumberFieldStructure<Integer>,
    // where does the generator of anf_polyquo get mapped to inside anf_quadratic?
    generator_image: QuadraticNumberFieldElement,
}

impl<
    AlgebraicNumberFieldPolynomialQuotientStructureBorrowed: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>,
> QuadraticNumberFieldIsomorphism<AlgebraicNumberFieldPolynomialQuotientStructureBorrowed>
{
    fn new(
        anf_polyquo: AlgebraicNumberFieldPolynomialQuotientStructureBorrowed,
    ) -> Result<Self, ()> {
        let anf_polyquo_borrowed = anf_polyquo.borrow();
        if anf_polyquo_borrowed.degree() == 2 {
            // let g be the generator of this ANF, so we are QQ[g]
            // g is a root of an integer polynomial
            // ax^2 + bx + c
            let poly = anf_polyquo_borrowed.modulus().primitive_part_fof();

            // find s, d such that s^2d = b^2-4ac and d is squarefree
            let (s, d) = {
                let disc = poly.discriminant().unwrap();
                let (mut d, powers) = disc.factor().unwrap().into_unit_and_powers();
                let mut s = Integer::ONE;
                for (p, k) in powers {
                    let (q, r) = k.div_mod(Natural::TWO);
                    s *= p.nat_pow(&q);
                    if r == Natural::ONE {
                        d *= p;
                    }
                }
                (s, d)
            };

            // now g = (-b ± s sqrt(d)) / 2a
            //       = (-b/2a) + (±s/2a) sqrt(d)
            // by the quadratic formula
            // this determines the desired isomorphism

            println!("{:?}, {:?}", s, d);

            compile_error!("todo");
        } else {
            Err(())
        }
    }
}

impl<
    AlgebraicNumberFieldPolynomialQuotientStructureBorrowed: BorrowedStructure<AlgebraicNumberFieldPolynomialQuotientStructure>,
> Morphism<AlgebraicNumberFieldPolynomialQuotientStructure, QuadraticNumberFieldStructure<Integer>>
    for QuadraticNumberFieldIsomorphism<AlgebraicNumberFieldPolynomialQuotientStructureBorrowed>
{
    fn domain(&self) -> &AlgebraicNumberFieldPolynomialQuotientStructure {
        self.anf_polyquo.borrow()
    }

    fn range(&self) -> &QuadraticNumberFieldStructure<Integer> {
        &self.anf_quadratic
    }
}

impl AlgebraicNumberFieldPolynomialQuotientStructure {
    /// Returns an isomorphism from this degree 2 polynomial quotient number field representation to the quadratic number field implementation .
    /// Returns `Err` if this number field is not quadratic.
    fn quadratic_anf_isomorphism<'a>(
        &'a self,
    ) -> Result<QuadraticNumberFieldIsomorphism<&'a Self>, ()> {
        QuadraticNumberFieldIsomorphism::new(self)
    }

    /// Returns an isomorphism from the quadratic number field implementation to this degree 2 polynomial quotient number field representation.
    /// Returns `Err` if this number field is not quadratic.
    fn into_quadratic_anf_isomorphism(self) -> Result<QuadraticNumberFieldIsomorphism<Self>, ()> {
        QuadraticNumberFieldIsomorphism::new(self)
    }
}

impl<D: BorrowedStructure<Integer>> QuadraticNumberFieldStructure<D> {}

#[cfg(test)]
mod tests {
    use crate::parsing::parse_rational_polynomial;

    use super::*;

    // ZZ[i]
    #[test]
    fn qanf_neg_1() {
        let poly = parse_rational_polynomial("x^2-6*x+13", "x").unwrap();
        println!("{}", poly);
        println!("{:?}", poly.clone().discriminant());
        let anf_poly = poly.algebraic_number_field().unwrap();

        let iso = anf_poly.quadratic_anf_isomorphism().unwrap();
    }
}
