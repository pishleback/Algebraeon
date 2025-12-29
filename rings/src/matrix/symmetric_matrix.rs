use std::fmt::Debug;

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

    pub fn at(&self, mut r: usize, mut c: usize) -> Result<&Set, MatOppErr> {
        if r >= self.n || c >= self.n {
            return Err(MatOppErr::InvalidIndex);
        }
        if r < c {
            (r, c) = (c, r);
        }
        debug_assert!(c <= r);
        Ok(&self.elems[r][c])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn make_and_get() {
        let mat = SymmetricMatrix::construct(3, |r, c| r + c);
        assert_eq!(*mat.at(0, 0).unwrap(), 0);
        assert_eq!(*mat.at(0, 1).unwrap(), 1);
        assert_eq!(*mat.at(0, 2).unwrap(), 2);
        assert_eq!(*mat.at(1, 0).unwrap(), 1);
        assert_eq!(*mat.at(1, 1).unwrap(), 2);
        assert_eq!(*mat.at(1, 2).unwrap(), 3);
        assert_eq!(*mat.at(2, 0).unwrap(), 2);
        assert_eq!(*mat.at(2, 1).unwrap(), 3);
        assert_eq!(*mat.at(2, 2).unwrap(), 4);
    }
}
