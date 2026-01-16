use std::{fmt::Debug, marker::PhantomData};

use algebraeon_sets::structure::{BorrowedStructure, EqSignature, SetSignature, Signature};

use crate::matrix::MatOppErr;

#[derive(Debug, Clone)]
pub struct SymmetricMatrix<Set: Clone> {
    n: usize,
    elems: Vec<Vec<Set>>, // index by [r][c] with r <= c
}

impl<Set: Debug + Clone + PartialEq> SymmetricMatrix<Set> {
    pub fn construct(n: usize, make_entry: impl Fn(usize, usize) -> Set) -> Self {
        let mut elems = Vec::with_capacity(n);
        for r in 0..n {
            let mut row = Vec::with_capacity(r);
            for c in 0..(r + 1) {
                let a = make_entry(r, c);
                #[cfg(debug_assertions)]
                {
                    let b = make_entry(c, r);
                    assert_eq!(a, b);
                }
                row.push(a);
            }
            elems.push(row);
        }
        Self { n, elems }
    }
}

impl<Set: Clone> SymmetricMatrix<Set> {
    pub fn n(&self) -> usize {
        self.n
    }

    pub fn filled(n: usize, s: Set) -> Self {
        Self::construct_top_right(n, |_, _| s.clone())
    }

    pub fn construct_top_right(n: usize, make_entry: impl Fn(usize, usize) -> Set) -> Self {
        let mut elems = Vec::with_capacity(n);
        for r in 0..n {
            let mut row = Vec::with_capacity(r);
            for c in r..n {
                let a = make_entry(r, c);
                row.push(a);
            }
            elems.push(row);
        }
        Self { n, elems }
    }

    pub fn construct_bottom_left(n: usize, make_entry: impl Fn(usize, usize) -> Set) -> Self {
        let mut elems = Vec::with_capacity(n);
        for r in 0..n {
            let mut row = Vec::with_capacity(r);
            for c in 0..(r + 1) {
                let a = make_entry(r, c);
                row.push(a);
            }
            elems.push(row);
        }
        Self { n, elems }
    }

    pub fn get(&self, mut r: usize, mut c: usize) -> Result<&Set, MatOppErr> {
        if r >= self.n || c >= self.n {
            return Err(MatOppErr::InvalidIndex);
        }
        if r < c {
            (r, c) = (c, r);
        }
        debug_assert!(c <= r);
        Ok(&self.elems[r][c])
    }

    pub fn get_mut(&mut self, mut r: usize, mut c: usize) -> Result<&mut Set, MatOppErr> {
        if r >= self.n || c >= self.n {
            return Err(MatOppErr::InvalidIndex);
        }
        if r < c {
            (r, c) = (c, r);
        }
        debug_assert!(c <= r);
        Ok(&mut self.elems[r][c])
    }

    pub fn set(&mut self, r: usize, c: usize, s: Set) -> Result<(), MatOppErr> {
        *self.get_mut(r, c)? = s;
        Ok(())
    }

    pub fn map<T: Clone>(self, f: impl Fn(Set) -> T) -> SymmetricMatrix<T> {
        SymmetricMatrix {
            n: self.n,
            elems: self
                .elems
                .into_iter()
                .map(|row| row.into_iter().map(&f).collect())
                .collect(),
        }
    }
}

impl<Set: Clone> SymmetricMatrix<Option<Set>> {
    pub fn unwrap_entries(self) -> Option<SymmetricMatrix<Set>> {
        Some(SymmetricMatrix {
            n: self.n,
            elems: self
                .elems
                .into_iter()
                .map(|row| row.into_iter().collect::<Option<Vec<Set>>>())
                .collect::<Option<Vec<Vec<Set>>>>()?,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SymmetricMatrixStructure<RS: SetSignature, RSB: BorrowedStructure<RS>> {
    _set: PhantomData<RS>,
    set: RSB,
}

impl<RS: SetSignature, RSB: BorrowedStructure<RS>> Signature for SymmetricMatrixStructure<RS, RSB> {}

impl<RS: SetSignature, RSB: BorrowedStructure<RS>> SetSignature
    for SymmetricMatrixStructure<RS, RSB>
{
    type Set = SymmetricMatrix<RS::Set>;

    fn is_element(&self, _x: &Self::Set) -> Result<(), String> {
        Ok(())
    }
}

impl<RS: SetSignature, RSB: BorrowedStructure<RS>> SymmetricMatrixStructure<RS, RSB> {
    pub fn new(set: RSB) -> Self {
        Self {
            _set: PhantomData,
            set,
        }
    }

    pub fn set(&self) -> &RS {
        self.set.borrow()
    }
}

pub trait ToSymmetrixMatricesSignature: SetSignature {
    fn symmetric_matrix_structure(&self) -> SymmetricMatrixStructure<Self, &Self> {
        SymmetricMatrixStructure::new(self)
    }

    fn into_symmetric_matrix_structure(self) -> SymmetricMatrixStructure<Self, Self> {
        SymmetricMatrixStructure::new(self)
    }
}
impl<RS: SetSignature> ToSymmetrixMatricesSignature for RS {}

impl<RS: EqSignature, RSB: BorrowedStructure<RS>> SymmetricMatrixStructure<RS, RSB> {
    pub fn equal(&self, a: &SymmetricMatrix<RS::Set>, b: &SymmetricMatrix<RS::Set>) -> bool {
        let n = a.n();
        if n != b.n() {
            false
        } else {
            for c in 0..n {
                for r in c..n {
                    if !self.set().equal(a.get(r, c).unwrap(), b.get(r, c).unwrap()) {
                        return false;
                    }
                }
            }
            true
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn make_and_get() {
        let mat = SymmetricMatrix::construct(3, |r, c| r + c);
        assert_eq!(*mat.get(0, 0).unwrap(), 0);
        assert_eq!(*mat.get(0, 1).unwrap(), 1);
        assert_eq!(*mat.get(0, 2).unwrap(), 2);
        assert_eq!(*mat.get(1, 0).unwrap(), 1);
        assert_eq!(*mat.get(1, 1).unwrap(), 2);
        assert_eq!(*mat.get(1, 2).unwrap(), 3);
        assert_eq!(*mat.get(2, 0).unwrap(), 2);
        assert_eq!(*mat.get(2, 1).unwrap(), 3);
        assert_eq!(*mat.get(2, 2).unwrap(), 4);
    }
}
