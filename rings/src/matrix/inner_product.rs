use crate::{
    matrix::{Matrix, RingMatricesSignature, SymmetricMatrix},
    structure::{
        ComplexConjugateSignature, ComplexSubsetSignature, OrderedRingSignature,
        RealSubsetSignature, RingSignature,
    },
};
use algebraeon_structures::*;
use std::{borrow::Borrow, sync::Arc};

pub trait ComplexInnerProduct<Ring: ComplexSubsetSignature> {
    /// # Panics
    /// If the dimensions of `a` and `b` do not match.
    fn inner_product(
        &self,
        a: &[impl Borrow<Ring::Elem>],
        b: &[impl Borrow<Ring::Elem>],
    ) -> Ring::Elem;
}

pub trait RealInnerProduct<Ring: RealSubsetSignature>: ComplexInnerProduct<Ring> {}

pub struct StandardInnerProduct<Ring: ComplexSubsetSignature + ComplexConjugateSignature> {
    ring: Arc<Ring>,
}

impl<Ring: ComplexSubsetSignature + ComplexConjugateSignature> StandardInnerProduct<Ring> {
    pub fn new(ring: Arc<Ring>) -> Arc<Self> {
        Self { ring }.into()
    }

    pub fn ring(&self) -> &Ring {
        self.ring.borrow()
    }
}

impl<Ring: ComplexSubsetSignature + ComplexConjugateSignature + RingSignature>
    ComplexInnerProduct<Ring> for StandardInnerProduct<Ring>
{
    fn inner_product(
        &self,
        a: &[impl Borrow<Ring::Elem>],
        b: &[impl Borrow<Ring::Elem>],
    ) -> Ring::Elem {
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

impl<Ring: RealSubsetSignature + RingSignature> RealInnerProduct<Ring>
    for StandardInnerProduct<Ring>
{
}

pub struct RealSymmetricInnerProduct<Ring: RealSubsetSignature> {
    ring: Arc<Ring>,
    mat: SymmetricMatrix<Ring::Elem>, // symmetric and positive-definite
}

fn is_positive_definite<Ring: OrderedRingSignature + RealSubsetSignature>(
    ring: &Ring,
    mat: &SymmetricMatrix<Ring::Elem>,
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

impl<Ring: RealSubsetSignature> RealSymmetricInnerProduct<Ring> {
    pub fn new(ring: Arc<Ring>, mat: SymmetricMatrix<Ring::Elem>) -> Arc<Self>
    where
        Ring: OrderedRingSignature,
    {
        debug_assert!(is_positive_definite(ring.borrow(), &mat));
        Self { ring, mat }.into()
    }

    pub fn ring(&self) -> &Ring {
        self.ring.borrow()
    }
}

impl<Ring: RingSignature + RealSubsetSignature> ComplexInnerProduct<Ring>
    for RealSymmetricInnerProduct<Ring>
{
    fn inner_product(
        &self,
        a: &[impl Borrow<Ring::Elem>],
        b: &[impl Borrow<Ring::Elem>],
    ) -> <Ring>::Elem {
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

impl<Ring: RealSubsetSignature + RingSignature> RealInnerProduct<Ring>
    for RealSymmetricInnerProduct<Ring>
{
}

// struct ComplexHermitianInnerProduct {
//     // Hermitian and positive-definite
//     mat: Matrix,
// }

// impl InnerProduct for ComplexHermitianInnerProduct {}
