use std::marker::PhantomData;

use crate::{polynomial::*, structure::*};
use algebraeon_sets::structure::*;
pub mod conway_finite_fields;
pub mod conway_polynomials;
pub mod extension;
pub mod modulo;
pub mod polynomial;
pub mod quaternary_field;

// #[derive(Debug, Clone, PartialEq, Eq)]
// struct IrreducibleFiniteFieldPolynomialsStructure<
//     FS: FiniteFieldSignature,
//     FSB: BorrowedStructure<FS>,
// > {
//     _field: PhantomData<FS>,
//     field: FSB,
// }

// impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>>
//     IrreducibleFiniteFieldPolynomialsStructure<FS, FSB>
// {
//     fn new(field: FSB) -> Self {
//         Self {
//             _field: PhantomData::default(),
//             field: field,
//         }
//     }

//     fn field(&self) -> &FS {
//         self.field.borrow()
//     }
// }

// impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> Signature
//     for IrreducibleFiniteFieldPolynomialsStructure<FS, FSB>
// {
// }

// impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> SetSignature
//     for IrreducibleFiniteFieldPolynomialsStructure<FS, FSB>
// {
//     type Set = Polynomial<FS::Set>;

//     fn is_element(&self, x: &Self::Set) -> bool {
//         todo!()
//     }
// }

// impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> EqSignature
//     for IrreducibleFiniteFieldPolynomialsStructure<FS, FSB>
// {
//     fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
//         self.field().polynomials().equal(a, b)
//     }
// }

// impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> OrdSignature
//     for IrreducibleFiniteFieldPolynomialsStructure<FS, FSB>
// {
//     fn cmp(&self, a: &Self::Set, b: &Self::Set) -> std::cmp::Ordering {
//         todo!()
//     }
// }

// impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> UniqueFactorizationSignature
//     for PolynomialStructure<FS, FSB>
// {
//     type Irreducibles = IrreducibleFiniteFieldPolynomialsStructure<FS, FS>;

//     type Factorizations<SelfB: BorrowedStructure<Self>> = FactoredRingElementStructure<Self, SelfB>;

//     fn factorizations<'a>(&'a self) -> Self::Factorizations<&'a Self> {
//         FactoredRingElementStructure::new(self)
//     }

//     fn into_factorizations(self) -> Self::Factorizations<Self> {
//         FactoredRingElementStructure::new(self)
//     }

//     fn irreducibles(&self) -> impl std::borrow::Borrow<Self::Irreducibles> {
//         IrreducibleFiniteFieldPolynomialsStructure::new(self.coeff_ring().clone())
//     }
// }

impl<FS: FiniteFieldSignature, FSB: BorrowedStructure<FS>> FactorableSignature
    for PolynomialStructure<FS, FSB>
{
    fn factor(&self, p: &Self::Set) -> Option<crate::structure::FactoredRingElement<Self::Set>> {
        Some(
            self.factorize_monic(p)?
                .factorize_squarefree()
                .factorize_distinct_degree()
                .factorize_cantor_zassenhaus(),
        )
    }
}
