use crate::{
    matrix::{Matrix, RingMatricesSignature, SymmetricMatrix},
    structure::{
        ComplexConjugateSignature, ComplexSubsetSignature, OrderedRingSignature,
        RealSubsetSignature, RingSignature,
    },
};
use algebraeon_sets::structure::BorrowedStructure;
use std::{borrow::Borrow, marker::PhantomData};

pub trait ComplexInnerProduct<Ring: ComplexSubsetSignature> {
    /// # Panics
    /// If the dimensions of `a` and `b` do not match.
    fn inner_product(
        &self,
        a: &[impl Borrow<Ring::Set>],
        b: &[impl Borrow<Ring::Set>],
    ) -> Ring::Set;
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
    fn inner_product(
        &self,
        a: &[impl Borrow<Ring::Set>],
        b: &[impl Borrow<Ring::Set>],
    ) -> Ring::Set {
        let n = a.len();
        assert_eq!(n, b.len());
        let mut t = self.ring().zero();
        for i in 0..n {
            self.ring().add_mut(
                &mut t,
                &self
                    .ring()
                    .mul(a[i].borrow(), &self.ring().conjugate(b[i].borrow())),
            );
        }
        t
    }
}

impl<Ring: RealSubsetSignature + RingSignature, RingB: BorrowedStructure<Ring>>
    RealInnerProduct<Ring> for StandardInnerProduct<Ring, RingB>
{
}

pub struct RealSymmetricInnerProduct<Ring: RealSubsetSignature, RingB: BorrowedStructure<Ring>> {
    _ring: PhantomData<Ring>,
    ring: RingB,
    mat: SymmetricMatrix<Ring::Set>, // symmetric and positive-definite
}

fn is_positive_definite<Ring: OrderedRingSignature + RealSubsetSignature>(
    ring: &Ring,
    mat: &SymmetricMatrix<Ring::Set>,
) -> bool {
    let ring_mat = ring.matrix_structure();
    let n = mat.n();
    let mat = Matrix::construct(n, n, |r, c| mat.get(r, c).unwrap().clone());
    for i in 1..n {
        if ring
            .cmp(
                &ring_mat
                    .det_naive(&mat.submatrix((0..i).collect(), (0..i).collect()))
                    .unwrap(),
                &ring.zero(),
            )
            .is_le()
        {
            return false;
        }
    }
    true
}

impl<Ring: RealSubsetSignature, RingB: BorrowedStructure<Ring>>
    RealSymmetricInnerProduct<Ring, RingB>
{
    pub fn new(ring: RingB, mat: SymmetricMatrix<Ring::Set>) -> Self
    where
        Ring: OrderedRingSignature,
    {
        debug_assert!(is_positive_definite(ring.borrow(), &mat));
        Self {
            _ring: PhantomData,
            ring,
            mat,
        }
    }

    pub fn ring(&self) -> &Ring {
        self.ring.borrow()
    }
}

impl<Ring: RingSignature + RealSubsetSignature, RingB: BorrowedStructure<Ring>>
    ComplexInnerProduct<Ring> for RealSymmetricInnerProduct<Ring, RingB>
{
    fn inner_product(
        &self,
        a: &[impl Borrow<Ring::Set>],
        b: &[impl Borrow<Ring::Set>],
    ) -> <Ring>::Set {
        let n = self.mat.n();
        assert_eq!(n, a.len());
        assert_eq!(n, b.len());
        let mut t = self.ring().zero();
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            for j in 0..n {
                self.ring().add_mut(
                    &mut t,
                    &self.ring().mul(
                        self.mat.get(i, j).unwrap(),
                        &self.ring().mul(a[i].borrow(), b[j].borrow()),
                    ),
                );
            }
        }
        t
    }
}

impl<Ring: RealSubsetSignature + RingSignature, RingB: BorrowedStructure<Ring>>
    RealInnerProduct<Ring> for RealSymmetricInnerProduct<Ring, RingB>
{
}

// struct ComplexHermitianInnerProduct {
//     // Hermitian and positive-definite
//     mat: Matrix,
// }

// impl InnerProduct for ComplexHermitianInnerProduct {}
