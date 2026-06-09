use crate::{
    matrix::{Matrix, MatrixStructure},
    structure::FieldSignature,
};
use algebraeon_structures::*;
use std::marker::PhantomData;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneralLinearStructure<
    RS: SetSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> {
    _ring: PhantomData<RS>,
    _ringb: PhantomData<RSB>,
    mats: MatB,
    n: usize,
}

impl<
    RS: SetSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> GeneralLinearStructure<RS, RSB, MatB>
{
    pub fn new(mats: MatB, n: usize) -> Self {
        Self {
            _ring: PhantomData,
            _ringb: PhantomData,
            mats,
            n,
        }
    }
}

impl<RS: SetSignature, RSB: BorrowedStructure<RS>> MatrixStructure<RS, RSB> {
    pub fn general_linear_structure(&self, n: usize) -> GeneralLinearStructure<RS, RSB, &Self> {
        GeneralLinearStructure::new(self, n)
    }

    pub fn into_general_linear_structure(self, n: usize) -> GeneralLinearStructure<RS, RSB, Self> {
        GeneralLinearStructure::new(self, n)
    }
}

impl<
    RS: SetSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> GeneralLinearStructure<RS, RSB, MatB>
{
    pub fn ring(&self) -> &RS {
        self.mats().ring()
    }

    pub fn mats(&self) -> &MatrixStructure<RS, RSB> {
        self.mats.borrow()
    }
}

impl<
    RS: SetSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> Signature for GeneralLinearStructure<RS, RSB, MatB>
{
}

impl<
    RS: FieldSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> SetSignature for GeneralLinearStructure<RS, RSB, MatB>
{
    type Elem = Matrix<RS::Elem>;

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        self.mats().validate_element(x)?;
        if x.rows() != self.n || x.cols() != self.n {
            return Err("Wrong dimension".to_string());
        }
        if !self.ring().is_unit(&self.mats().det(x.clone()).unwrap()) {
            return Err("Matrix is not invertible".to_string());
        }
        Ok(())
    }
}

impl<
    RS: FieldSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> IdentitySignature for GeneralLinearStructure<RS, RSB, MatB>
{
    fn identity(&self) -> Self::Elem {
        self.mats().ident(self.n)
    }
}

impl<
    RS: FieldSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> CompositionSignature for GeneralLinearStructure<RS, RSB, MatB>
{
    fn compose(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.mats().mul(a, b).unwrap()
    }
}

impl<
    RS: FieldSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> AssociativeCompositionSignature for GeneralLinearStructure<RS, RSB, MatB>
{
}

impl<
    RS: FieldSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> MonoidSignature for GeneralLinearStructure<RS, RSB, MatB>
{
}

impl<
    RS: FieldSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> TryInverseSignature for GeneralLinearStructure<RS, RSB, MatB>
{
    fn try_inverse(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<
    RS: FieldSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> TryLeftInverseSignature for GeneralLinearStructure<RS, RSB, MatB>
{
    fn try_left_inverse(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<
    RS: FieldSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> TryRightInverseSignature for GeneralLinearStructure<RS, RSB, MatB>
{
    fn try_right_inverse(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<
    RS: FieldSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> LeftCancellativeCompositionSignature for GeneralLinearStructure<RS, RSB, MatB>
{
    fn try_left_difference(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.compose(&self.inverse(b), a))
    }
}

impl<
    RS: FieldSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> RightCancellativeCompositionSignature for GeneralLinearStructure<RS, RSB, MatB>
{
    fn try_right_difference(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.compose(a, &self.inverse(b)))
    }
}

impl<
    RS: FieldSignature,
    RSB: BorrowedStructure<RS>,
    MatB: BorrowedStructure<MatrixStructure<RS, RSB>>,
> GroupSignature for GeneralLinearStructure<RS, RSB, MatB>
{
    fn inverse(&self, a: &Self::Elem) -> Self::Elem {
        self.mats().inv(a.clone()).unwrap()
    }
}
