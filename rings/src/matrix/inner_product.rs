use crate::{
    matrix::{Matrix, RingMatricesSignature},
    structure::{
        ComplexConjugateSignature, ComplexSubsetSignature, RealSubsetSignature, RingSignature,
    },
};
use algebraeon_sets::structure::BorrowedStructure;
use std::marker::PhantomData;

pub trait ComplexInnerProduct<Ring: ComplexSubsetSignature> {
    /// # Panics
    /// If the dimensions of `a` and `b` do not match.
    fn inner_product(&self, a: &Matrix<Ring::Set>, b: &Matrix<Ring::Set>) -> Ring::Set;
}

pub trait RealInnerProduct<Ring: RealSubsetSignature>: ComplexInnerProduct<Ring> {}

pub struct StandardInnerProduct<
    Ring: ComplexSubsetSignature + ComplexConjugateSignature,
    RingB: BorrowedStructure<Ring>,
> {
    _ring: PhantomData<Ring>,
    ring: RingB,
}

impl<Ring: ComplexSubsetSignature + ComplexConjugateSignature, RingB: BorrowedStructure<Ring>>
    StandardInnerProduct<Ring, RingB>
{
    pub fn new(ring: RingB) -> Self {
        Self {
            _ring: PhantomData,
            ring,
        }
    }

    pub fn ring(&self) -> &Ring {
        self.ring.borrow()
    }
}

impl<
    Ring: ComplexSubsetSignature + ComplexConjugateSignature + RingSignature,
    RingB: BorrowedStructure<Ring>,
> ComplexInnerProduct<Ring> for StandardInnerProduct<Ring, RingB>
{
    fn inner_product(&self, a: &Matrix<Ring::Set>, b: &Matrix<Ring::Set>) -> Ring::Set {
        self.ring()
            .matrix_structure()
            .dot(a, &self.ring().matrix_structure().conjugate(b))
    }
}

impl<Ring: RealSubsetSignature + RingSignature, RingB: BorrowedStructure<Ring>>
    RealInnerProduct<Ring> for StandardInnerProduct<Ring, RingB>
{
}

// struct RealSymmetricInnerProduct {
//     // symmetric and positive-definite
//     mat: SymmetricMatrix,
// }

// impl InnerProduct for RealSymmetricInnerProduct {}

// struct ComplexHermitianInnerProduct {
//     // Hermitian and positive-definite
//     mat: Matrix,
// }

// impl InnerProduct for ComplexHermitianInnerProduct {}
