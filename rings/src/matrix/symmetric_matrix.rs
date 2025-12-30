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
