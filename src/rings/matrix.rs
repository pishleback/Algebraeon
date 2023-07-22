#![allow(dead_code)]

use super::rings::*;

#[derive(Debug)]
pub enum MatOppErr {
    DimMissmatch,
    InvalidIndex,
    NotSquare,
}

#[derive(Debug, Clone)]
pub struct Matrix<R: ComRing> {
    dim1: usize,
    dim2: usize,
    transpose: bool,
    elems: Vec<R>, //length self.rows * self.cols. row r and column c is index c + r * self.cols
}

impl<R: ComRing> PartialEq for Matrix<R> {
    fn eq(&self, other: &Self) -> bool {
        let rows = self.rows();
        let cols = self.cols();
        if rows != other.rows() {
            false
        } else if cols != other.cols() {
            false
        } else {
            for c in 0..cols {
                for r in 0..rows {
                    if self.at(r, c).unwrap() != other.at(r, c).unwrap() {
                        return false;
                    }
                }
            }
            true
        }
    }
}

impl<R: ComRing> Matrix<R> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        if self.elems.len() != self.dim1 * self.dim2 {
            return Err("matrix entries has the wrong length");
        }
        Ok(())
    }

    pub fn from_rows(rows_elems: Vec<Vec<R>>) -> Self {
        let rows = rows_elems.len();
        assert!(rows >= 1);
        let cols = rows_elems[0].len();
        for r in 1..rows {
            assert_eq!(rows_elems[r].len(), cols);
        }
        let mut mat = Self::zero(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                *mat.at_mut(r, c).unwrap() = rows_elems[r][c].clone()
            }
        }
        mat
    }

    pub fn rows(&self) -> usize {
        match self.transpose {
            false => self.dim1,
            true => self.dim2,
        }
    }

    pub fn cols(&self) -> usize {
        match self.transpose {
            false => self.dim2,
            true => self.dim1,
        }
    }

    fn rc_to_idx(&self, r: usize, c: usize) -> usize {
        match self.transpose {
            false => c + r * self.dim2,
            true => r + c * self.dim2,
        }
    }

    pub fn at(&self, r: usize, c: usize) -> Result<&R, MatOppErr> {
        if r >= self.rows() {
            Err(MatOppErr::InvalidIndex)
        } else if c >= self.cols() {
            Err(MatOppErr::InvalidIndex)
        } else {
            let idx = self.rc_to_idx(r, c);
            Ok(&self.elems[idx])
        }
    }

    pub fn at_mut(&mut self, r: usize, c: usize) -> Result<&mut R, MatOppErr> {
        if r >= self.rows() {
            Err(MatOppErr::InvalidIndex)
        } else if c >= self.cols() {
            Err(MatOppErr::InvalidIndex)
        } else {
            let idx = self.rc_to_idx(r, c);
            Ok(&mut self.elems[idx])
        }
    }

    pub fn zero(rows: usize, cols: usize) -> Self {
        let mut elems = Vec::with_capacity(rows * cols);
        for _i in 0..rows * cols {
            elems.push(R::zero());
        }
        Self {
            dim1: rows,
            dim2: cols,
            transpose: false,
            elems,
        }
    }

    pub fn ident(n: usize) -> Self {
        let mut elems = Vec::with_capacity(n * n);
        for r in 0..n {
            for c in 0..n {
                match r == c {
                    true => elems.push(R::one()),
                    false => elems.push(R::zero()),
                }
            }
        }
        Self {
            dim1: n,
            dim2: n,
            transpose: false,
            elems,
        }
    }
}

impl<R: ComRing + std::fmt::Display> Matrix<R> {
    pub fn pprint(&self) {
        let mut str_rows = vec![];
        for r in 0..self.rows() {
            str_rows.push(vec![]);
            for c in 0..self.cols() {
                str_rows[r].push(self.at(r, c).unwrap().to_string());
            }
        }
        let cols_widths: Vec<usize> = (0..self.cols())
            .map(|c| {
                (0..self.rows())
                    .map(|r| str_rows[r][c].len())
                    .fold(0usize, |a, b| a.max(b))
            })
            .collect();

        for r in 0..self.rows() {
            for c in 0..self.cols() {
                while str_rows[r][c].len() < cols_widths[c] {
                    str_rows[r][c].push(' ');
                }
            }
        }
        for r in 0..self.rows() {
            for c in 0..self.cols() {
                if c != 0 {
                    print!("    ");
                }
                print!("{}", str_rows[r][c]);
            }
            print!("\n");
        }
    }
}

impl<R: ComRing> Matrix<R> {
    pub fn transpose(mut self) -> Self {
        self.transpose_mut();
        self
    }

    pub fn transpose_ref(&self) -> Self {
        self.clone().transpose()
    }

    pub fn transpose_mut(&mut self) {
        self.transpose = !self.transpose;
    }

    pub fn add_mut(&mut self, other: &Self) -> Result<(), MatOppErr> {
        if self.rows() != other.rows() || self.cols() != other.cols() {
            Err(MatOppErr::DimMissmatch)
        } else {
            let rows = self.rows();
            let cols = self.cols();
            println!();
            for c in 0..cols {
                for r in 0..rows {
                    println!(
                        "{} {} {:?} {:?}",
                        r,
                        c,
                        self.at(r, c),
                        other.rc_to_idx(r, c)
                    );
                    self.at_mut(r, c).unwrap().add_mut(other.at(r, c).unwrap())
                }
            }
            Ok(())
        }
    }

    pub fn add(mut a: Self, b: Self) -> Result<Self, MatOppErr> {
        match a.add_mut(&b) {
            Ok(()) => Ok(a),
            Err(e) => Err(e),
        }
    }

    pub fn add_ref(mut a: Self, b: &Self) -> Result<Self, MatOppErr> {
        match a.add_mut(b) {
            Ok(()) => Ok(a),
            Err(e) => Err(e),
        }
    }

    pub fn add_refs(a: &Self, b: &Self) -> Result<Self, MatOppErr> {
        let mut new_a = a.clone();
        match new_a.add_mut(&b) {
            Ok(()) => Ok(new_a),
            Err(e) => Err(e),
        }
    }

    // pub fn mul(a: Self, b: Self) -> Result<Self, MatOppErr> {
    //     Self::mul_refs(&a, &b)
    // }

    // pub fn mul_lref(a: &Self, b: Self) -> Result<Self, MatOppErr> {
    //     Self::mul_refs(a, &b)
    // }

    // pub fn mul_rref(a: Self, b: &Self) -> Result<Self, MatOppErr> {
    //     Self::mul_refs(&a, b)
    // }

    pub fn mul_refs(a: &Self, b: &Self) -> Result<Self, MatOppErr> {
        let mids = a.cols();
        if mids != b.rows() {
            return Err(MatOppErr::DimMissmatch);
        }
        let rows = a.rows();
        let cols = b.cols();
        let mut s = Matrix::<R>::zero(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                for m in 0..mids {
                    s.at_mut(r, c)
                        .unwrap()
                        .add_mut(&R::mul_refs(a.at(r, m).unwrap(), b.at(m, c).unwrap()));
                }
            }
        }
        Ok(s)
    }

    pub fn det_naive(&self) -> Result<R, MatOppErr> {
        let n = self.dim1;
        if n != self.dim2 {
            Err(MatOppErr::NotSquare)
        } else {
            let mut det = R::zero();
            for perm in super::super::sets::permutations::all_perms(n) {
                let mut prod = R::one();
                for k in 0..n {
                    prod.mul_mut(self.at(k, perm.call(k).unwrap()).unwrap());
                }
                if !perm.sign() {
                    prod.neg_mut();
                }
                det.add_mut(&prod);
            }
            Ok(det)
        }
    }
}

impl<R: ED> Matrix<R> {
    pub fn row_hermite_algorithm(&self) -> (Self, Self) {
        todo!();
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use super::*;

    #[test]
    fn invariant() {
        let m: Matrix<Integer> = Matrix {
            dim1: 3,
            dim2: 4,
            transpose: false,
            elems: vec![
                Integer::from(1),
                Integer::from(2),
                Integer::from(3),
                Integer::from(4),
                Integer::from(5),
            ],
        };
        match m.check_invariants() {
            Ok(()) => panic!(),
            Err(_) => {}
        }

        let m: Matrix<Integer> = Matrix {
            dim1: 2,
            dim2: 3,
            transpose: true,
            elems: vec![
                Integer::from(1),
                Integer::from(2),
                Integer::from(3),
                Integer::from(4),
                Integer::from(5),
                Integer::from(6),
            ],
        };
        m.check_invariants().unwrap();
    }

    #[test]
    fn transpose_eq() {
        let a: Matrix<Integer> = Matrix {
            dim1: 2,
            dim2: 2,
            transpose: false,
            elems: vec![
                Integer::from(0),
                Integer::from(1),
                Integer::from(2),
                Integer::from(3),
            ],
        };
        a.check_invariants().unwrap();

        let b: Matrix<Integer> = Matrix {
            dim1: 2,
            dim2: 2,
            transpose: true,
            elems: vec![
                Integer::from(0),
                Integer::from(2),
                Integer::from(1),
                Integer::from(3),
            ],
        };
        b.check_invariants().unwrap();

        assert_eq!(a, b);
    }

    #[test]
    fn add() {
        {
            let mut a: Matrix<Integer> = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                elems: vec![
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(4),
                    Integer::from(5),
                    Integer::from(6),
                ],
            };
            a.check_invariants().unwrap();

            let b: Matrix<Integer> = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                elems: vec![
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(1),
                    Integer::from(2),
                ],
            };
            b.check_invariants().unwrap();

            let c: Matrix<Integer> = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                elems: vec![
                    Integer::from(2),
                    Integer::from(4),
                    Integer::from(4),
                    Integer::from(6),
                    Integer::from(6),
                    Integer::from(8),
                ],
            };
            c.check_invariants().unwrap();

            a.add_mut(&b).unwrap();

            assert_eq!(a, c);
        }

        {
            let mut a: Matrix<Integer> = Matrix {
                dim1: 3,
                dim2: 2,
                transpose: false,
                elems: vec![
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(4),
                    Integer::from(5),
                    Integer::from(6),
                ],
            };
            a.check_invariants().unwrap();

            let b: Matrix<Integer> = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: true,
                elems: vec![
                    Integer::from(10),
                    Integer::from(20),
                    Integer::from(30),
                    Integer::from(40),
                    Integer::from(50),
                    Integer::from(60),
                ],
            };
            b.check_invariants().unwrap();

            let c: Matrix<Integer> = Matrix {
                dim1: 3,
                dim2: 2,
                transpose: false,
                elems: vec![
                    Integer::from(11),
                    Integer::from(42),
                    Integer::from(23),
                    Integer::from(54),
                    Integer::from(35),
                    Integer::from(66),
                ],
            };
            c.check_invariants().unwrap();

            a.add_mut(&b).unwrap();

            assert_eq!(a, c);
        }

        {
            let mut a: Matrix<Integer> = Matrix {
                dim1: 3,
                dim2: 2,
                transpose: false,
                elems: vec![
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(4),
                    Integer::from(5),
                    Integer::from(6),
                ],
            };
            a.check_invariants().unwrap();

            let b: Matrix<Integer> = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                elems: vec![
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(1),
                    Integer::from(2),
                ],
            };
            b.check_invariants().unwrap();

            match a.add_mut(&b) {
                Ok(()) => panic!(),
                Err(MatOppErr::DimMissmatch) => {}
                Err(_) => panic!(),
            }
        }

        {
            let a: Matrix<Integer> = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                elems: vec![
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(4),
                    Integer::from(5),
                    Integer::from(6),
                ],
            };
            a.check_invariants().unwrap();

            let b: Matrix<Integer> = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                elems: vec![
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(1),
                    Integer::from(2),
                ],
            };
            b.check_invariants().unwrap();

            let c: Matrix<Integer> = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                elems: vec![
                    Integer::from(2),
                    Integer::from(4),
                    Integer::from(4),
                    Integer::from(6),
                    Integer::from(6),
                    Integer::from(8),
                ],
            };
            c.check_invariants().unwrap();

            assert_eq!(Matrix::<Integer>::add(a, b).unwrap(), c);
        }
    }

    #[test]
    fn mul() {
        {
            let a: Matrix<Integer> = Matrix {
                dim1: 2,
                dim2: 4,
                transpose: false,
                elems: vec![
                    Integer::from(3),
                    Integer::from(2),
                    Integer::from(1),
                    Integer::from(5),
                    Integer::from(9),
                    Integer::from(1),
                    Integer::from(3),
                    Integer::from(0),
                ],
            };
            a.check_invariants().unwrap();

            let b: Matrix<Integer> = Matrix {
                dim1: 4,
                dim2: 3,
                transpose: false,
                elems: vec![
                    Integer::from(2),
                    Integer::from(9),
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(3),
                    Integer::from(5),
                    Integer::from(2),
                    Integer::from(4),
                    Integer::from(7),
                    Integer::from(8),
                    Integer::from(1),
                    Integer::from(5),
                ],
            };
            b.check_invariants().unwrap();

            let c: Matrix<Integer> = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                elems: vec![
                    Integer::from(50),
                    Integer::from(42),
                    Integer::from(42),
                    Integer::from(25),
                    Integer::from(96),
                    Integer::from(26),
                ],
            };
            c.check_invariants().unwrap();

            assert_eq!(Matrix::<Integer>::mul_refs(&a, &b).unwrap(), c);
        }
    }

    #[test]
    fn det_naive() {
        let m: Matrix<Integer> = Matrix::from_rows(vec![
            vec![Integer::from(1), Integer::from(3), Integer::from(2)],
            vec![Integer::from(-3), Integer::from(-1), Integer::from(-3)],
            vec![Integer::from(2), Integer::from(3), Integer::from(1)],
        ]);
        assert_eq!(m.det_naive().unwrap(), Integer::from(-15));
    }
}
