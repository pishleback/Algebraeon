use crate::{
    matrix::{Matrix, MatrixStructure},
    structure::FieldSignature,
};
use algebraeon_structures::*;
use std::rc::Rc;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneralLinearStructure<RS: SetSignature> {
    mats: Rc<MatrixStructure<RS>>,
    n: usize,
}

impl<RS: SetSignature> GeneralLinearStructure<RS> {
    pub fn new(mats: Rc<MatrixStructure<RS>>, n: usize) -> Rc<Self> {
        Self { mats, n }.into()
    }
}

impl<RS: SetSignature> MatrixStructure<RS> {
    pub fn general_linear_structure(self: &Rc<Self>, n: usize) -> Rc<GeneralLinearStructure<RS>> {
        GeneralLinearStructure::new(self.clone(), n)
    }
}

impl<RS: SetSignature> GeneralLinearStructure<RS> {
    pub fn ring(&self) -> Rc<RS> {
        self.mats().ring()
    }

    pub fn mats(&self) -> &Rc<MatrixStructure<RS>> {
        &self.mats
    }
}

impl<RS: SetSignature> Signature for GeneralLinearStructure<RS> {}

impl<RS: FieldSignature> SetSignature for GeneralLinearStructure<RS> {
    type Elem = Matrix<RS::Elem>;

    fn validate_element(self: &Rc<Self>, x: &Self::Elem) -> Result<(), String> {
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

impl<RS: FieldSignature> IdentitySignature for GeneralLinearStructure<RS> {
    fn identity(self: &Rc<Self>) -> Self::Elem {
        self.mats().ident(self.n)
    }
}

impl<RS: FieldSignature> CompositionSignature for GeneralLinearStructure<RS> {
    fn compose(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        self.mats().mul(a, b).unwrap()
    }
}

impl<RS: FieldSignature> AssociativeCompositionSignature for GeneralLinearStructure<RS> {}

impl<RS: FieldSignature> MonoidSignature for GeneralLinearStructure<RS> {}

impl<RS: FieldSignature> TryInverseSignature for GeneralLinearStructure<RS> {
    fn try_inverse(self: &Rc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<RS: FieldSignature> TryLeftInverseSignature for GeneralLinearStructure<RS> {
    fn try_left_inverse(self: &Rc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<RS: FieldSignature> TryRightInverseSignature for GeneralLinearStructure<RS> {
    fn try_right_inverse(self: &Rc<Self>, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.inverse(a))
    }
}

impl<RS: FieldSignature> LeftCancellativeCompositionSignature for GeneralLinearStructure<RS> {
    fn try_left_difference(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.compose(&self.inverse(b), a))
    }
}

impl<RS: FieldSignature> RightCancellativeCompositionSignature for GeneralLinearStructure<RS> {
    fn try_right_difference(self: &Rc<Self>, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.compose(a, &self.inverse(b)))
    }
}

impl<RS: FieldSignature> GroupSignature for GeneralLinearStructure<RS> {
    fn inverse(self: &Rc<Self>, a: &Self::Elem) -> Self::Elem {
        self.mats().inv(a.clone()).unwrap()
    }
}
