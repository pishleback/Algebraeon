use crate::structure::*;
use algebraeon_nzq::Natural;
use algebraeon_sets::structure::*;
use std::{borrow::Borrow, marker::PhantomData};

#[derive(Debug)]
pub enum MatOppErr {
    DimMissmatch,
    InvalidIndex,
    NotSquare,
    Singular,
}

#[derive(Debug, Clone)]
pub struct Matrix<Set: Clone> {
    dim1: usize,
    dim2: usize,
    transpose: bool,
    flip_rows: bool,
    flip_cols: bool,
    elems: Vec<Set>, //length self.rows * self.cols. row r and column c is index c + r * self.cols
}

impl<Set: Clone> Matrix<Set> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        if self.elems.len() != self.dim1 * self.dim2 {
            return Err("matrix entries has the wrong length");
        }
        Ok(())
    }

    pub fn full(rows: usize, cols: usize, elem: &Set) -> Self {
        let mut elems = Vec::with_capacity(rows * cols);
        for _i in 0..rows * cols {
            elems.push(elem.clone());
        }
        Self {
            dim1: rows,
            dim2: cols,
            transpose: false,
            flip_rows: false,
            flip_cols: false,
            elems,
        }
    }

    pub fn construct(rows: usize, cols: usize, make_entry: impl Fn(usize, usize) -> Set) -> Self {
        let mut elems = Vec::with_capacity(rows * cols);
        for idx in 0..rows * cols {
            let (r, c) = (idx / cols, idx % cols); //idx_to_rc for transpose=false
            elems.push(make_entry(r, c).clone());
        }
        Self {
            dim1: rows,
            dim2: cols,
            transpose: false,
            flip_rows: false,
            flip_cols: false,
            elems,
        }
    }

    pub fn from_rows(rows_elems: Vec<Vec<impl Into<Set> + Clone>>) -> Self {
        let rows = rows_elems.len();
        assert!(rows >= 1);
        let cols = rows_elems[0].len();
        for r in 1..rows {
            assert_eq!(rows_elems[r].len(), cols);
        }
        Self::construct(rows, cols, |r, c| rows_elems[r][c].clone().into())
    }

    pub fn from_cols(cols_elems: Vec<Vec<impl Into<Set> + Clone>>) -> Self {
        Self::from_rows(cols_elems).transpose()
    }

    fn rc_to_idx(&self, mut r: usize, mut c: usize) -> usize {
        if self.flip_rows {
            r = self.rows() - r - 1;
        }

        if self.flip_cols {
            c = self.cols() - c - 1;
        }

        if self.transpose {
            r + c * self.dim2
        } else {
            c + r * self.dim2
        }
    }

    pub fn at(&self, r: usize, c: usize) -> Result<&Set, MatOppErr> {
        if r >= self.rows() || c >= self.cols() {
            Err(MatOppErr::InvalidIndex)
        } else {
            let idx = self.rc_to_idx(r, c);
            Ok(&self.elems[idx])
        }
    }

    pub fn at_mut(&mut self, r: usize, c: usize) -> Result<&mut Set, MatOppErr> {
        if r >= self.rows() || c >= self.cols() {
            Err(MatOppErr::InvalidIndex)
        } else {
            let idx = self.rc_to_idx(r, c);
            Ok(&mut self.elems[idx])
        }
    }

    pub fn rows(&self) -> usize {
        if self.transpose { self.dim2 } else { self.dim1 }
    }

    pub fn cols(&self) -> usize {
        if self.transpose { self.dim1 } else { self.dim2 }
    }

    pub fn submatrix(&self, rows: Vec<usize>, cols: Vec<usize>) -> Self {
        let mut elems = vec![];
        for r in &rows {
            for c in &cols {
                elems.push(self.at(*r, *c).unwrap().clone());
            }
        }
        Matrix {
            dim1: rows.len(),
            dim2: cols.len(),
            transpose: false,
            flip_rows: false,
            flip_cols: false,
            elems,
        }
    }

    pub fn get_row_submatrix(&self, row: usize) -> Self {
        self.submatrix(vec![row], (0..self.cols()).collect())
    }

    pub fn get_col_submatrix(&self, col: usize) -> Self {
        self.submatrix((0..self.rows()).collect(), vec![col])
    }

    pub fn get_row_refs(&self, row: usize) -> Vec<&Set> {
        assert!(row < self.rows());
        (0..self.cols()).map(|c| self.at(row, c).unwrap()).collect()
    }

    pub fn get_col_refs(&self, col: usize) -> Vec<&Set> {
        assert!(col < self.cols());
        (0..self.rows()).map(|r| self.at(r, col).unwrap()).collect()
    }

    pub fn get_row(&self, row: usize) -> Vec<Set> {
        assert!(row < self.rows());
        self.get_row_refs(row).into_iter().cloned().collect()
    }

    pub fn get_col(&self, col: usize) -> Vec<Set> {
        assert!(col < self.cols());
        self.get_col_refs(col).into_iter().cloned().collect()
    }

    pub fn apply_map<NewSet: Clone>(&self, f: impl Fn(&Set) -> NewSet) -> Matrix<NewSet> {
        Matrix {
            dim1: self.dim1,
            dim2: self.dim2,
            transpose: self.transpose,
            flip_rows: self.flip_rows,
            flip_cols: self.flip_cols,
            elems: self.elems.iter().map(f).collect(),
        }
    }

    pub fn transpose(mut self) -> Self {
        self.transpose_mut();
        self
    }
    pub fn transpose_ref(&self) -> Self {
        self.clone().transpose()
    }
    pub fn transpose_mut(&mut self) {
        self.transpose = !self.transpose;
        (self.flip_rows, self.flip_cols) = (self.flip_cols, self.flip_rows);
    }

    pub fn flip_rows(mut self) -> Self {
        self.flip_rows_mut();
        self
    }
    pub fn flip_rows_ref(&self) -> Self {
        self.clone().flip_rows()
    }
    pub fn flip_rows_mut(&mut self) {
        self.flip_rows = !self.flip_rows;
    }

    pub fn flip_cols(mut self) -> Self {
        self.flip_cols_mut();
        self
    }
    pub fn flip_cols_ref(&self) -> Self {
        self.clone().flip_cols()
    }
    pub fn flip_cols_mut(&mut self) {
        self.flip_cols = !self.flip_cols;
    }

    pub fn join_rows<MatT: Borrow<Matrix<Set>>>(cols: usize, mats: Vec<MatT>) -> Matrix<Set> {
        let mut rows = 0;
        for mat in &mats {
            assert_eq!(cols, mat.borrow().cols());
            rows += mat.borrow().rows();
        }
        Matrix::construct(rows, cols, |r, c| {
            //todo use a less cursed method
            let mut row_offset = 0;
            for mat in &mats {
                for mr in 0..mat.borrow().rows() {
                    for mc in 0..cols {
                        if r == row_offset + mr && c == mc {
                            return mat.borrow().at(mr, mc).unwrap().clone();
                        }
                    }
                }
                row_offset += mat.borrow().rows();
            }
            panic!();
        })
    }

    pub fn join_cols<MatT: Borrow<Matrix<Set>>>(rows: usize, mats: Vec<MatT>) -> Matrix<Set> {
        let mut t_mats = vec![];
        for mat in mats {
            t_mats.push(mat.borrow().clone().transpose());
        }
        let joined = Self::join_rows(rows, t_mats.iter().collect());
        joined.transpose()
    }

    pub fn entries_list(&self) -> Vec<&Set> {
        let mut entries = vec![];
        for r in 0..self.rows() {
            for c in 0..self.cols() {
                entries.push(self.at(r, c).unwrap());
            }
        }
        entries
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixStructure<RS: SetSignature, RSB: BorrowedStructure<RS>> {
    _ring: PhantomData<RS>,
    ring: RSB,
}

impl<RS: SetSignature, RSB: BorrowedStructure<RS>> Signature for MatrixStructure<RS, RSB> {}

impl<RS: SetSignature, RSB: BorrowedStructure<RS>> SetSignature for MatrixStructure<RS, RSB> {
    type Set = Matrix<RS::Set>;

    fn is_element(&self, _x: &Self::Set) -> bool {
        true
    }
}

impl<RS: SetSignature, RSB: BorrowedStructure<RS>> MatrixStructure<RS, RSB> {
    pub fn new(ring: RSB) -> Self {
        Self {
            _ring: PhantomData::default(),
            ring,
        }
    }

    pub fn ring(&self) -> &RS {
        self.ring.borrow()
    }
}

impl<RS: EqSignature, RSB: BorrowedStructure<RS>> MatrixStructure<RS, RSB> {
    pub fn equal(&self, a: &Matrix<RS::Set>, b: &Matrix<RS::Set>) -> bool {
        let rows = a.rows();
        let cols = a.cols();
        if rows != b.rows() || cols != b.cols() {
            false
        } else {
            for c in 0..cols {
                for r in 0..rows {
                    if !self.ring().equal(a.at(r, c).unwrap(), b.at(r, c).unwrap()) {
                        return false;
                    }
                }
            }
            true
        }
    }
}

impl<RS: ToStringSignature, RSB: BorrowedStructure<RS>> MatrixStructure<RS, RSB> {
    pub fn pprint(&self, mat: &Matrix<RS::Set>) {
        let mut str_rows = vec![];
        for r in 0..mat.rows() {
            str_rows.push(vec![]);
            for c in 0..mat.cols() {
                str_rows[r].push(self.ring().to_string(mat.at(r, c).unwrap()));
            }
        }
        let cols_widths: Vec<usize> = (0..mat.cols())
            .map(|c| {
                (0..mat.rows())
                    .map(|r| str_rows[r][c].chars().count())
                    .fold(0usize, |a, b| a.max(b))
            })
            .collect();

        for r in 0..mat.rows() {
            for c in 0..mat.cols() {
                while str_rows[r][c].chars().count() < cols_widths[c] {
                    str_rows[r][c].push(' ');
                }
                debug_assert_eq!(str_rows[r][c].chars().count(), cols_widths[c]);
            }
        }

        for r in 0..mat.rows() {
            if mat.rows() == 1 {
                print!("( ");
            } else if r == 0 {
                print!("/ ");
            } else if r == mat.rows() - 1 {
                print!("\\ ");
            } else {
                print!("| ");
            }
            for c in 0..mat.cols() {
                if c != 0 {
                    print!("    ");
                }
                print!("{}", str_rows[r][c]);
            }
            if mat.rows() == 1 {
                print!(" )");
            } else if r == 0 {
                print!(" \\");
            } else if r == mat.rows() - 1 {
                print!(" /");
            } else {
                print!(" |");
            }
            println!();
        }
    }
}

impl<RS: RingSignature, RSB: BorrowedStructure<RS>> MatrixStructure<RS, RSB> {
    pub fn zero(&self, rows: usize, cols: usize) -> Matrix<RS::Set> {
        Matrix::construct(rows, cols, |_r, _c| self.ring().zero())
    }

    pub fn ident(&self, n: usize) -> Matrix<RS::Set> {
        Matrix::construct(n, n, |r, c| {
            if r == c {
                self.ring().one()
            } else {
                self.ring().zero()
            }
        })
    }

    pub fn diag(&self, diag: &Vec<RS::Set>) -> Matrix<RS::Set> {
        Matrix::construct(diag.len(), diag.len(), |r, c| {
            if r == c {
                diag[r].clone()
            } else {
                self.ring().zero()
            }
        })
    }

    pub fn join_diag<MatT: Borrow<Matrix<RS::Set>>>(&self, mats: Vec<MatT>) -> Matrix<RS::Set> {
        if mats.is_empty() {
            Matrix::construct(0, 0, |_r, _c| unreachable!())
        } else if mats.len() == 1 {
            mats[0].borrow().clone()
        } else {
            let i = mats.len() / 2;
            let (first, last) = mats.split_at(i);
            let first = self.join_diag(first.into_iter().map(|m| m.borrow()).collect());
            let last = self.join_diag(last.into_iter().map(|m| m.borrow()).collect());
            Matrix::construct(
                first.rows() + last.rows(),
                first.cols() + last.cols(),
                |r, c| {
                    if r < first.rows() && c < first.cols() {
                        first.at(r, c).unwrap().clone()
                    } else if first.rows() <= r && first.cols() <= c {
                        last.at(r - first.rows(), c - first.cols()).unwrap().clone()
                    } else {
                        self.ring().zero()
                    }
                },
            )
        }
    }

    pub fn dot(&self, a: &Matrix<RS::Set>, b: &Matrix<RS::Set>) -> RS::Set {
        let rows = a.rows();
        let cols = a.cols();
        assert_eq!(rows, b.rows());
        assert_eq!(cols, b.cols());
        let mut tot = self.ring().zero();
        for r in 0..rows {
            for c in 0..cols {
                self.ring().add_mut(
                    &mut tot,
                    &self.ring().mul(a.at(r, c).unwrap(), b.at(r, c).unwrap()),
                );
            }
        }
        tot
    }

    pub fn add_mut(&self, a: &mut Matrix<RS::Set>, b: &Matrix<RS::Set>) -> Result<(), MatOppErr> {
        if a.rows() != b.rows() || a.cols() != b.cols() {
            Err(MatOppErr::DimMissmatch)
        } else {
            let rows = a.rows();
            let cols = a.cols();
            for c in 0..cols {
                for r in 0..rows {
                    self.ring()
                        .add_mut(a.at_mut(r, c).unwrap(), b.at(r, c).unwrap());
                }
            }
            Ok(())
        }
    }

    pub fn add(
        &self,
        a: &Matrix<RS::Set>,
        b: &Matrix<RS::Set>,
    ) -> Result<Matrix<RS::Set>, MatOppErr> {
        let mut new_a = a.clone();
        match self.add_mut(&mut new_a, b) {
            Ok(()) => Ok(new_a),
            Err(e) => Err(e),
        }
    }

    pub fn neg_mut(&self, a: &mut Matrix<RS::Set>) {
        for r in 0..a.rows() {
            for c in 0..a.cols() {
                let neg_elem = self.ring().neg(a.at(r, c).unwrap());
                *a.at_mut(r, c).unwrap() = neg_elem;
            }
        }
    }

    pub fn neg(&self, mut a: Matrix<RS::Set>) -> Matrix<RS::Set> {
        self.neg_mut(&mut a);
        a
    }

    pub fn mul(
        &self,
        a: &Matrix<RS::Set>,
        b: &Matrix<RS::Set>,
    ) -> Result<Matrix<RS::Set>, MatOppErr> {
        let mids = a.cols();
        if mids != b.rows() {
            return Err(MatOppErr::DimMissmatch);
        }
        let rows = a.rows();
        let cols = b.cols();
        let mut s = self.zero(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                for m in 0..mids {
                    self.ring().add_mut(
                        s.at_mut(r, c).unwrap(),
                        &self.ring().mul(a.at(r, m).unwrap(), b.at(m, c).unwrap()),
                    );
                }
            }
        }
        Ok(s)
    }

    pub fn apply_row(&self, mat: &Matrix<RS::Set>, row: &Vec<RS::Set>) -> Vec<RS::Set> {
        assert_eq!(mat.rows(), row.len());
        (0..mat.cols())
            .map(|c| {
                self.ring().sum(
                    (0..mat.rows())
                        .map(|r| self.ring().mul(mat.at(r, c).unwrap(), &row[r]))
                        .collect(),
                )
            })
            .collect()
    }

    pub fn apply_col(&self, mat: &Matrix<RS::Set>, col: &Vec<RS::Set>) -> Vec<RS::Set> {
        assert_eq!(mat.cols(), col.len());
        (0..mat.rows())
            .map(|r| {
                self.ring().sum(
                    (0..mat.cols())
                        .map(|c| self.ring().mul(mat.at(r, c).unwrap(), &col[c]))
                        .collect(),
                )
            })
            .collect()
    }

    pub fn mul_scalar(&self, mut a: Matrix<RS::Set>, scalar: &RS::Set) -> Matrix<RS::Set> {
        for r in 0..a.rows() {
            for c in 0..a.cols() {
                self.ring().mul_mut(a.at_mut(r, c).unwrap(), scalar);
            }
        }
        a
    }

    pub fn mul_scalar_ref(&self, a: &Matrix<RS::Set>, scalar: &RS::Set) -> Matrix<RS::Set> {
        self.mul_scalar(a.clone(), scalar)
    }

    pub fn det_naive(&self, a: &Matrix<RS::Set>) -> Result<RS::Set, MatOppErr> {
        let n = a.rows();
        if n != a.cols() {
            Err(MatOppErr::NotSquare)
        } else {
            let mut det = self.ring().zero();
            for perm in algebraeon_groups::permutation::Permutation::all_permutations(n) {
                let mut prod = self.ring().one();
                for k in 0..n {
                    self.ring()
                        .mul_mut(&mut prod, a.at(k, perm.call(k)).unwrap());
                }
                match perm.sign() {
                    algebraeon_groups::examples::c2::C2::Identity => {}
                    algebraeon_groups::examples::c2::C2::Flip => {
                        prod = self.ring().neg(&prod);
                    }
                }

                self.ring().add_mut(&mut det, &prod);
            }
            Ok(det)
        }
    }

    pub fn trace(&self, a: &Matrix<RS::Set>) -> Result<RS::Set, MatOppErr> {
        let n = a.rows();
        if n != a.cols() {
            Err(MatOppErr::NotSquare)
        } else {
            Ok(self
                .ring()
                .sum((0..n).map(|i| a.at(i, i).unwrap()).collect()))
        }
    }

    pub fn nat_pow(&self, a: &Matrix<RS::Set>, k: &Natural) -> Result<Matrix<RS::Set>, MatOppErr> {
        let n = a.rows();
        if n != a.cols() {
            Err(MatOppErr::NotSquare)
        } else if *k == Natural::ZERO {
            Ok(self.ident(n))
        } else if *k == Natural::ONE {
            Ok(a.clone())
        } else {
            debug_assert!(*k >= Natural::TWO);
            let bits: Vec<_> = k.bits().collect();
            let mut pows = vec![a.clone()];
            while pows.len() < bits.len() {
                pows.push(
                    self.mul(pows.last().unwrap(), pows.last().unwrap())
                        .unwrap(),
                );
            }
            let count = bits.len();
            debug_assert_eq!(count, pows.len());
            let mut ans = self.ident(n);
            for i in 0..count {
                if bits[i] {
                    ans = self.mul(&ans, &pows[i]).unwrap();
                }
            }
            Ok(ans)
        }
    }
}

impl<R: MetaType> MetaType for Matrix<R>
where
    R::Signature: SetSignature,
{
    type Signature = MatrixStructure<R::Signature, R::Signature>;

    fn structure() -> Self::Signature {
        MatrixStructure::new(R::structure())
    }
}

impl<R: MetaType> Matrix<R>
where
    R::Signature: ToStringSignature,
{
    pub fn pprint(&self) {
        Self::structure().pprint(self);
    }
}

impl<R: MetaType> PartialEq for Matrix<R>
where
    R::Signature: RingSignature,
{
    fn eq(&self, other: &Self) -> bool {
        Self::structure().equal(self, other)
    }
}

impl<R: MetaType> Eq for Matrix<R> where R::Signature: RingSignature {}

impl<R: MetaType> Matrix<R>
where
    R::Signature: RingSignature,
{
    pub fn zero(rows: usize, cols: usize) -> Self {
        Self::structure().zero(rows, cols)
    }

    pub fn ident(n: usize) -> Self {
        Self::structure().ident(n)
    }

    pub fn diag(diag: &Vec<R>) -> Self {
        Self::structure().diag(diag)
    }

    pub fn dot(a: &Self, b: &Self) -> R {
        Self::structure().dot(a, b)
    }

    pub fn add_mut(&mut self, b: &Self) -> Result<(), MatOppErr> {
        Self::structure().add_mut(self, b)
    }

    pub fn add(a: &Self, b: &Self) -> Result<Self, MatOppErr> {
        Self::structure().add(a, b)
    }

    pub fn neg_mut(&mut self) {
        Self::structure().neg_mut(self);
    }

    pub fn neg(&self) -> Self {
        Self::structure().neg(self.clone())
    }

    pub fn mul(a: &Self, b: &Self) -> Result<Self, MatOppErr> {
        Self::structure().mul(a, b)
    }

    pub fn apply_row(&self, row: &Vec<R>) -> Vec<R> {
        Self::structure().apply_row(self, row)
    }

    pub fn apply_col(&self, col: &Vec<R>) -> Vec<R> {
        Self::structure().apply_col(self, col)
    }

    pub fn mul_scalar(&self, scalar: &R) -> Matrix<R> {
        Self::structure().mul_scalar(self.clone(), scalar)
    }

    pub fn mul_scalar_ref(&self, scalar: &R) -> Matrix<R> {
        Self::structure().mul_scalar_ref(self, scalar)
    }

    pub fn det_naive(&self) -> Result<R, MatOppErr> {
        Self::structure().det_naive(self)
    }

    pub fn trace(&self) -> Result<R, MatOppErr> {
        Self::structure().trace(self)
    }
}

#[cfg(test)]
mod tests {
    use algebraeon_nzq::Integer;

    use super::*;

    #[test]
    fn test_join_rows() {
        let top = Matrix::<Integer>::from_rows(vec![vec![1, 2, 3], vec![4, 5, 6]]);
        let bot = Matrix::from_rows(vec![vec![7, 8, 9]]);

        let both = Matrix::from_rows(vec![vec![1, 2, 3], vec![4, 5, 6], vec![7, 8, 9]]);

        println!("top");
        top.pprint();
        println!("bot");
        bot.pprint();
        println!("both");
        both.pprint();

        let ans = Matrix::join_rows(3, vec![top, bot]);
        println!("ans");
        ans.pprint();

        assert_eq!(ans, both);
    }

    #[test]
    fn invariants() {
        let m = Matrix {
            dim1: 3,
            dim2: 4,
            transpose: false,
            flip_rows: false,
            flip_cols: false,
            elems: vec![
                Integer::from(1),
                Integer::from(2),
                Integer::from(3),
                Integer::from(4),
                Integer::from(5),
            ],
        };
        if let Ok(()) = m.check_invariants() {
            panic!();
        }

        let m = Matrix {
            dim1: 2,
            dim2: 3,
            transpose: true,
            flip_rows: false,
            flip_cols: false,
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
        let a = Matrix {
            dim1: 2,
            dim2: 2,
            transpose: false,
            flip_rows: false,
            flip_cols: false,
            elems: vec![
                Integer::from(0),
                Integer::from(1),
                Integer::from(2),
                Integer::from(3),
            ],
        };
        a.check_invariants().unwrap();

        let b = Matrix {
            dim1: 2,
            dim2: 2,
            transpose: true,
            flip_rows: false,
            flip_cols: false,
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
    fn flip_axes_eq() {
        let mut a = Matrix::<Integer>::from_rows(vec![vec![1, 2], vec![3, 4]]);
        a.pprint();
        println!("flip rows");
        a.flip_rows_mut();
        a.pprint();
        assert_eq!(
            a,
            Matrix::from_rows(vec![
                vec![Integer::from(3), Integer::from(4)],
                vec![Integer::from(1), Integer::from(2)],
            ])
        );
        println!("transpose");
        a.transpose_mut();
        a.pprint();
        assert_eq!(
            a,
            Matrix::from_rows(vec![
                vec![Integer::from(3), Integer::from(1)],
                vec![Integer::from(4), Integer::from(2)],
            ])
        );
        println!("flip rows");
        a.flip_rows_mut();
        a.pprint();
        assert_eq!(
            a,
            Matrix::from_rows(vec![
                vec![Integer::from(4), Integer::from(2)],
                vec![Integer::from(3), Integer::from(1)],
            ])
        );
        println!("flip cols");
        a.flip_cols_mut();
        a.pprint();
        assert_eq!(
            a,
            Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(4)],
                vec![Integer::from(1), Integer::from(3)],
            ])
        );
        println!("transpose");
        a.transpose_mut();
        a.pprint();
        assert_eq!(
            a,
            Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(1)],
                vec![Integer::from(4), Integer::from(3)],
            ])
        );
        println!("flip cols");
        a.flip_cols_mut();
        a.pprint();
        assert_eq!(
            a,
            Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(2)],
                vec![Integer::from(3), Integer::from(4)],
            ])
        );
    }

    #[test]
    fn add() {
        {
            let mut a = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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

            let b = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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

            let c = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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
            let mut a = Matrix {
                dim1: 3,
                dim2: 2,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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

            let b = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: true,
                flip_rows: false,
                flip_cols: false,
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

            let c = Matrix {
                dim1: 3,
                dim2: 2,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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
            let mut a = Matrix {
                dim1: 3,
                dim2: 2,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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

            let b = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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
            let a = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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

            let b = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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

            let c = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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

            assert_eq!(Matrix::add(&a, &b).unwrap(), c);
        }
    }

    #[test]
    fn mul() {
        {
            let a = Matrix {
                dim1: 2,
                dim2: 4,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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

            let b = Matrix {
                dim1: 4,
                dim2: 3,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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

            let c = Matrix {
                dim1: 2,
                dim2: 3,
                transpose: false,
                flip_rows: false,
                flip_cols: false,
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

            assert_eq!(Matrix::mul(&a, &b).unwrap(), c);
        }
    }

    #[test]
    fn matrix_apply_row_and_col_test() {
        let m = Matrix::<Integer>::from_rows(vec![
            vec![Integer::from(1), Integer::from(2), Integer::from(3)],
            vec![Integer::from(6), Integer::from(5), Integer::from(4)],
        ]);

        assert_eq!(
            m.apply_row(&vec![Integer::from(1), Integer::from(0)]),
            vec![Integer::from(1), Integer::from(2), Integer::from(3)]
        );

        assert_eq!(
            m.apply_row(&vec![Integer::from(0), Integer::from(1)]),
            vec![Integer::from(6), Integer::from(5), Integer::from(4)]
        );

        assert_eq!(
            m.apply_row(&vec![Integer::from(1), Integer::from(1)]),
            vec![Integer::from(7), Integer::from(7), Integer::from(7)]
        );

        assert_eq!(
            m.apply_col(&vec![Integer::from(1), Integer::from(0), Integer::from(0)]),
            vec![Integer::from(1), Integer::from(6)]
        );

        assert_eq!(
            m.apply_col(&vec![Integer::from(0), Integer::from(1), Integer::from(0)]),
            vec![Integer::from(2), Integer::from(5)]
        );

        assert_eq!(
            m.apply_col(&vec![Integer::from(0), Integer::from(0), Integer::from(1)]),
            vec![Integer::from(3), Integer::from(4)]
        );

        assert_eq!(
            m.apply_col(&vec![Integer::from(1), Integer::from(1), Integer::from(1)]),
            vec![Integer::from(6), Integer::from(15)]
        );
    }

    #[test]
    fn det_naive() {
        let m = Matrix::<Integer>::from_rows(vec![
            vec![Integer::from(1), Integer::from(3)],
            vec![Integer::from(4), Integer::from(2)],
        ]);
        println!("{}", m.det_naive().unwrap());
        assert_eq!(m.det_naive().unwrap(), Integer::from(-10));

        let m = Matrix::<Integer>::from_rows(vec![
            vec![Integer::from(1), Integer::from(3), Integer::from(2)],
            vec![Integer::from(-3), Integer::from(-1), Integer::from(-3)],
            vec![Integer::from(2), Integer::from(3), Integer::from(1)],
        ]);
        println!("{}", m.det_naive().unwrap());
        assert_eq!(m.det_naive().unwrap(), Integer::from(-15));
    }
}
