use std::borrow::Borrow;
use std::rc::Rc;

use itertools::Itertools;
use malachite_nz::natural::Natural;

use super::super::ring_structure::structure::*;
use super::super::structure::*;
use super::subspace::*;
use crate::polynomial::polynomial::*;

#[derive(Debug)]
pub enum MatOppErr {
    DimMissmatch,
    InvalidIndex,
    NotSquare,
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

        match self.transpose {
            false => c + r * self.dim2,
            true => r + c * self.dim2,
        }
    }

    pub fn at(&self, r: usize, c: usize) -> Result<&Set, MatOppErr> {
        if r >= self.rows() {
            Err(MatOppErr::InvalidIndex)
        } else if c >= self.cols() {
            Err(MatOppErr::InvalidIndex)
        } else {
            let idx = self.rc_to_idx(r, c);
            Ok(&self.elems[idx])
        }
    }

    pub fn at_mut(&mut self, r: usize, c: usize) -> Result<&mut Set, MatOppErr> {
        if r >= self.rows() {
            Err(MatOppErr::InvalidIndex)
        } else if c >= self.cols() {
            Err(MatOppErr::InvalidIndex)
        } else {
            let idx = self.rc_to_idx(r, c);
            Ok(&mut self.elems[idx])
        }
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

    pub fn get_row(&self, row: usize) -> Self {
        self.submatrix(vec![row], (0..self.cols()).collect())
    }

    pub fn get_col(&self, col: usize) -> Self {
        self.submatrix((0..self.rows()).collect(), vec![col])
    }

    pub fn apply_map<NewSet: Clone>(&self, f: impl Fn(&Set) -> NewSet) -> Matrix<NewSet> {
        Matrix {
            dim1: self.dim1,
            dim2: self.dim2,
            transpose: self.transpose,
            flip_rows: self.flip_rows,
            flip_cols: self.flip_cols,
            elems: self.elems.iter().map(|x| f(x)).collect(),
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
pub struct MatrixStructure<RS: Structure> {
    ring: Rc<RS>,
}

impl<RS: Structure> Structure for MatrixStructure<RS> {
    type Set = Matrix<RS::Set>;
}

impl<RS: Structure> MatrixStructure<RS> {
    pub fn new(ring: Rc<RS>) -> Self {
        Self { ring }
    }

    pub fn ring(&self) -> Rc<RS> {
        self.ring.clone()
    }
}

impl<RS: EqualityStructure> MatrixStructure<RS> {
    fn equal(&self, a: &Matrix<RS::Set>, b: &Matrix<RS::Set>) -> bool {
        let rows = a.rows();
        let cols = a.cols();
        if rows != b.rows() {
            false
        } else if cols != b.cols() {
            false
        } else {
            for c in 0..cols {
                for r in 0..rows {
                    if !self.ring.equal(a.at(r, c).unwrap(), b.at(r, c).unwrap()) {
                        return false;
                    }
                }
            }
            true
        }
    }
}

impl<RS: DisplayableStructure> MatrixStructure<RS> {
    pub fn pprint(&self, mat: &Matrix<RS::Set>) {
        let mut str_rows = vec![];
        for r in 0..mat.rows() {
            str_rows.push(vec![]);
            for c in 0..mat.cols() {
                str_rows[r].push(self.ring.elem_to_string(mat.at(r, c).unwrap()));
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
            print!("\n");
        }
    }
}

impl<RS: RingStructure> MatrixStructure<RS> {
    pub fn zero(&self, rows: usize, cols: usize) -> Matrix<RS::Set> {
        Matrix::construct(rows, cols, |_r, _c| self.ring.zero())
    }

    pub fn ident(&self, n: usize) -> Matrix<RS::Set> {
        Matrix::construct(n, n, |r, c| {
            if r == c {
                self.ring.one()
            } else {
                self.ring.zero()
            }
        })
    }

    pub fn diag(&self, diag: &Vec<RS::Set>) -> Matrix<RS::Set> {
        Matrix::construct(diag.len(), diag.len(), |r, c| {
            if r == c {
                diag[r].clone()
            } else {
                self.ring.zero()
            }
        })
    }

    pub fn join_diag<MatT: Borrow<Matrix<RS::Set>>>(&self, mats: Vec<MatT>) -> Matrix<RS::Set> {
        if mats.len() == 0 {
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
        let mut tot = self.ring.zero();
        for r in 0..rows {
            for c in 0..cols {
                self.ring.add_mut(
                    &mut tot,
                    &self.ring.mul(a.at(r, c).unwrap(), b.at(r, c).unwrap()),
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
                    self.ring
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
        match self.add_mut(&mut new_a, &b) {
            Ok(()) => Ok(new_a),
            Err(e) => Err(e),
        }
    }

    pub fn neg_mut(&self, a: &mut Matrix<RS::Set>) {
        for r in 0..a.rows() {
            for c in 0..a.cols() {
                let neg_elem = self.ring.neg(a.at(r, c).unwrap());
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
                    self.ring.add_mut(
                        s.at_mut(r, c).unwrap(),
                        &self.ring.mul(a.at(r, m).unwrap(), b.at(m, c).unwrap()),
                    );
                }
            }
        }
        Ok(s)
    }

    pub fn mul_scalar(&self, mut a: Matrix<RS::Set>, scalar: &RS::Set) -> Matrix<RS::Set> {
        for r in 0..a.rows() {
            for c in 0..a.cols() {
                self.ring.mul_mut(a.at_mut(r, c).unwrap(), scalar);
            }
        }
        a
    }

    pub fn mul_scalar_ref(&self, a: &Matrix<RS::Set>, scalar: &RS::Set) -> Matrix<RS::Set> {
        self.mul_scalar(a.clone(), scalar)
    }

    pub fn det_naive(&self, a: &Matrix<RS::Set>) -> Result<RS::Set, MatOppErr> {
        let n = a.dim1;
        if n != a.dim2 {
            Err(MatOppErr::NotSquare)
        } else {
            let mut det = self.ring.zero();
            for perm in orthoclase_groups::permutation::Permutation::all_permutations(n) {
                let mut prod = self.ring.one();
                for k in 0..n {
                    self.ring.mul_mut(&mut prod, a.at(k, perm.call(k)).unwrap());
                }
                match perm.sign() {
                    orthoclase_groups::examples::c2::C2::Ident => {}
                    orthoclase_groups::examples::c2::C2::Flip => {
                        prod = self.ring.neg(&prod);
                    }
                }

                self.ring.add_mut(&mut det, &prod);
            }
            Ok(det)
        }
    }

    pub fn trace(&self, a: &Matrix<RS::Set>) -> Result<RS::Set, MatOppErr> {
        let n = a.dim1;
        if n != a.dim2 {
            Err(MatOppErr::NotSquare)
        } else {
            Ok(self.ring.sum((0..n).map(|i| a.at(i, i).unwrap()).collect()))
        }
    }

    pub fn nat_pow(&self, a: &Matrix<RS::Set>, k: &Natural) -> Matrix<RS::Set> {
        use malachite_base::num::logic::traits::BitIterable;
        let n = a.rows();
        assert_eq!(n, a.cols());
        if *k == 0 {
            self.ident(n)
        } else if *k == 1 {
            a.clone()
        } else {
            debug_assert!(*k >= 2);
            let bits: Vec<_> = k.bits().collect();
            let mut pows = vec![a.clone()];
            while pows.len() < bits.len() {
                pows.push(
                    self.mul(&pows.last().unwrap(), &pows.last().unwrap())
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
            ans
        }
    }
}

#[derive(Debug)]
enum ElementaryOppType<RS: RingStructure> {
    //swap distinct rows
    Swap(usize, usize),
    //multiply a row by a unit
    UnitMul {
        row: usize,
        unit: RS::Set,
    },
    //row(i) -> row(i) + x*row(j)
    AddRowMul {
        i: usize,
        j: usize,
        x: RS::Set,
    },
    //apply invertible row operations to two rows
    // /a b\
    // \c d/
    //such that ad-bc is a unit
    TwoInv {
        i: usize,
        j: usize,
        a: RS::Set,
        b: RS::Set,
        c: RS::Set,
        d: RS::Set,
    },
}

struct ElementaryOpp<RS: RingStructure> {
    ring: Rc<RS>,
    transpose: bool, //false = row opp, true = column opp
    opp: ElementaryOppType<RS>,
}

impl<RS: BezoutDomainStructure> ElementaryOpp<RS> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        match &self.opp {
            ElementaryOppType::Swap(i, j) => {
                if i == j {
                    return Err("can only swap distinct rows");
                }
            }
            ElementaryOppType::AddRowMul { i, j, x: _x } => {
                if i == j {
                    return Err("can only add a multiple of a row to a distinct row");
                }
            }
            ElementaryOppType::UnitMul { row: _row, unit } => {
                if !self.ring.is_unit(unit) {
                    return Err("can only multiply a row by a unit");
                }
            }
            ElementaryOppType::TwoInv { i, j, a, b, c, d } => {
                if i == j {
                    return Err("rows must be distinct");
                }
                let m = Matrix::<RS::Set> {
                    dim1: 2,
                    dim2: 2,
                    transpose: false,
                    flip_rows: false,
                    flip_cols: false,
                    elems: vec![a.clone(), b.clone(), c.clone(), d.clone()],
                };
                if !self.ring.is_unit(
                    &MatrixStructure::new(self.ring.clone())
                        .det_naive(&m)
                        .unwrap(),
                ) {
                    return Err("can only apply an invertible row opperation to two rows");
                }
            }
        }
        Ok(())
    }

    fn new_row_opp(ring: Rc<RS>, opp: ElementaryOppType<RS>) -> Self {
        Self {
            ring,
            transpose: false,
            opp,
        }
    }

    fn new_col_opp(ring: Rc<RS>, opp: ElementaryOppType<RS>) -> Self {
        Self {
            ring,
            transpose: true,
            opp,
        }
    }

    fn det(&self) -> RS::Set {
        match &self.opp {
            ElementaryOppType::Swap(_i, _j) => self.ring.neg(&self.ring.one()),
            ElementaryOppType::UnitMul { row: _row, unit } => unit.clone(),
            ElementaryOppType::AddRowMul {
                i: _i,
                j: _j,
                x: _x,
            } => self.ring.one(),
            ElementaryOppType::TwoInv {
                i: _i,
                j: _j,
                a,
                b,
                c,
                d,
            } => self
                .ring
                .add(&self.ring.mul(a, d), &self.ring.neg(&self.ring.mul(b, c))),
        }
    }

    fn apply(&self, m: &mut Matrix<RS::Set>) {
        debug_assert!(self.check_invariants().is_ok());
        if self.transpose {
            m.transpose_mut();
        }
        match &self.opp {
            // /0 1\
            // \1 0/
            ElementaryOppType::Swap(i, j) => {
                for col in 0..m.cols() {
                    let tmp = m.at(*i, col).unwrap().clone();
                    *m.at_mut(*i, col).unwrap() = m.at(*j, col).unwrap().clone();
                    *m.at_mut(*j, col).unwrap() = tmp;
                }
            }
            // /1 x\
            // \0 1/
            ElementaryOppType::AddRowMul { i, j, x } => {
                for col in 0..m.cols() {
                    let offset = self.ring.mul(m.at(*j, col).unwrap(), x);
                    self.ring.add_mut(m.at_mut(*i, col).unwrap(), &offset)
                }
            }
            // /u 0\
            // \0 1/
            ElementaryOppType::UnitMul { row, unit } => {
                for col in 0..m.cols() {
                    self.ring.mul_mut(m.at_mut(*row, col).unwrap(), unit)
                }
            }
            // /a b\
            // \c d/
            ElementaryOppType::TwoInv { i, j, a, b, c, d } => {
                for col in 0..m.cols() {
                    // tmp = c*row(i) + d*row(j)
                    let tmp = self.ring.add(
                        &self.ring.mul(c, m.at(*i, col).unwrap()),
                        &self.ring.mul(d, m.at(*j, col).unwrap()),
                    );
                    // row(i) = a*row(i) + b*row(j)
                    *m.at_mut(*i, col).unwrap() = self.ring.add(
                        &self.ring.mul(a, m.at(*i, col).unwrap()),
                        &self.ring.mul(b, m.at(*j, col).unwrap()),
                    );
                    // row(j) = tmp
                    *m.at_mut(*j, col).unwrap() = tmp;
                }
            }
        };
        if self.transpose {
            m.transpose_mut();
        }
    }
}

impl<RS: BezoutDomainStructure> MatrixStructure<RS> {
    pub fn row_span(&self, a: Matrix<RS::Set>) -> LinearLattice<RS::Set> {
        LinearLatticeStructure::new(self.ring.clone()).from_span(
            1,
            a.cols(),
            (0..a.rows())
                .map(|r| a.submatrix(vec![r], (0..a.cols()).collect()))
                .collect(),
        )
    }

    pub fn col_span(&self, a: Matrix<RS::Set>) -> LinearLattice<RS::Set> {
        LinearLatticeStructure::new(self.ring.clone()).from_span(
            a.rows(),
            1,
            (0..a.cols())
                .map(|c| a.submatrix((0..a.rows()).collect(), vec![c]))
                .collect(),
        )
    }

    pub fn row_affine_span(&self, a: Matrix<RS::Set>) -> AffineLattice<RS::Set> {
        let affine_lattice_structure = AffineLatticeStructure::new(self.ring.clone());
        if a.rows() == 0 {
            affine_lattice_structure.empty(1, a.cols())
        } else {
            let offset = a.get_row(0);

            let b = Matrix::construct(a.rows() - 1, a.cols(), |r, c| {
                self.ring.add(
                    &self.ring.neg(offset.at(0, c).unwrap()),
                    a.at(r + 1, c).unwrap(),
                )
            });

            let linlat = self.row_span(b);

            affine_lattice_structure.from_offset_and_linear_lattice(1, a.cols(), offset, linlat)
        }
    }

    pub fn col_affine_span(&self, a: Matrix<RS::Set>) -> AffineLattice<RS::Set> {
        let affine_lattice_structure = AffineLatticeStructure::new(self.ring.clone());
        if a.cols() == 0 {
            affine_lattice_structure.empty(a.rows(), 1)
        } else {
            let offset = a.get_col(0);

            let b = Matrix::construct(a.rows(), a.cols() - 1, |r, c| {
                self.ring.add(
                    &self.ring.neg(offset.at(r, 0).unwrap()),
                    a.at(r, c + 1).unwrap(),
                )
            });

            let linlat = self.col_span(b);

            affine_lattice_structure.from_offset_and_linear_lattice(a.rows(), 1, offset, linlat)
        }
    }

    pub fn row_kernel(&self, a: Matrix<RS::Set>) -> LinearLattice<RS::Set> {
        let (_h, u, _u_det, pivs) = self.row_hermite_algorithm(a);
        LinearLatticeStructure::new(self.ring.clone()).from_basis(
            1,
            u.cols(),
            (pivs.len()..u.rows())
                .into_iter()
                .map(|r| u.submatrix(vec![r], (0..u.cols()).collect()))
                .collect(),
        )
    }

    pub fn col_kernel(&self, a: Matrix<RS::Set>) -> LinearLattice<RS::Set> {
        let (_h, u, _u_det, pivs) = self.col_hermite_algorithm(a);
        LinearLatticeStructure::new(self.ring.clone()).from_basis(
            u.rows(),
            1,
            (pivs.len()..u.cols())
                .into_iter()
                .map(|c| u.submatrix((0..u.rows()).collect(), vec![c]))
                .collect(),
        )
    }

    pub fn row_solve(
        &self,
        m: &Matrix<RS::Set>,
        y: impl Borrow<Matrix<RS::Set>>,
    ) -> Option<Matrix<RS::Set>> {
        match self.col_solve(&m.transpose_ref(), &y.borrow().transpose_ref()) {
            Some(x) => Some(x.transpose()),
            None => None,
        }
    }

    pub fn col_solve(
        &self,
        m: &Matrix<RS::Set>,
        y: impl Borrow<Matrix<RS::Set>>,
    ) -> Option<Matrix<RS::Set>> {
        assert_eq!(y.borrow().rows(), m.rows());
        assert_eq!(y.borrow().cols(), 1);
        //the kernel of ext_mat is related to the solution
        let ext_mat = Matrix::join_cols(m.rows(), vec![y.borrow(), m]);
        //we are looking for a point in the column kernel where the first coordinate is 1
        let col_ker = self.col_kernel(ext_mat);

        let first_coords: Vec<&RS::Set> = (0..LinearLatticeStructure::new(self.ring.clone())
            .rank(&col_ker))
            .map(|basis_num| {
                LinearLatticeStructure::new(self.ring.clone())
                    .basis_matrix_element(&col_ker, basis_num, 0, 0)
            })
            .collect();

        let (g, taps) = self.ring.xgcd_list(first_coords);

        if self.ring.is_unit(&g) {
            debug_assert!(self.ring.equal(&g, &self.ring.one()));
        }
        if self.ring.equal(&g, &self.ring.one()) {
            //there is a solution
            //it is given by -(sum(taps * col_ker.basis)) with the first coordinate (equal to 1) removed
            let mut ext_ans = self.zero(m.cols() + 1, 1);
            for basis_num in 0..LinearLatticeStructure::new(self.ring.clone()).rank(&col_ker) {
                self.add_mut(
                    &mut ext_ans,
                    &self.mul_scalar_ref(
                        &LinearLatticeStructure::new(self.ring.clone())
                            .basis_matrix(&col_ker, basis_num),
                        &taps[basis_num],
                    ),
                )
                .unwrap();
            }
            debug_assert!(self.ring.equal(ext_ans.at(0, 0).unwrap(), &self.ring.one()));
            let x = self.neg(ext_ans.submatrix((1..ext_ans.rows()).collect(), vec![0]));
            debug_assert!(self.equal(&self.mul(m, &x).unwrap(), y.borrow()));
            Some(x)
        } else {
            None //there is no solution
        }
    }

    pub fn row_solution_lattice(
        &self,
        m: &Matrix<RS::Set>,
        y: impl Borrow<Matrix<RS::Set>>,
    ) -> AffineLattice<RS::Set> {
        match self.row_solve(m, y) {
            Some(x) => AffineLatticeStructure::new(self.ring.clone())
                .from_offset_and_linear_lattice(1, m.rows(), x, self.row_kernel(m.clone())),
            None => AffineLatticeStructure::new(self.ring.clone()).empty(1, m.rows()),
        }
    }

    pub fn col_solution_lattice(
        &self,
        m: &Matrix<RS::Set>,
        y: impl Borrow<Matrix<RS::Set>>,
    ) -> AffineLattice<RS::Set> {
        match self.col_solve(m, y) {
            Some(x) => AffineLatticeStructure::new(self.ring.clone())
                .from_offset_and_linear_lattice(m.cols(), 1, x, self.col_kernel(m.clone())),
            None => AffineLatticeStructure::new(self.ring.clone()).empty(m.cols(), 1),
        }
    }

    //if A:=self return (H, U, u_det, pivots) such that
    //H is in row hermite normal form
    //U is invertible
    //H=UA
    //u det is the determinant of u. It is a unit
    //pivots[r] is the column of the rth pivot and pivots.len() == rank(A)
    pub fn row_hermite_algorithm(
        &self,
        mut m: Matrix<RS::Set>,
    ) -> (Matrix<RS::Set>, Matrix<RS::Set>, RS::Set, Vec<usize>) {
        //build up U by applying row opps to the identity as we go
        let mut u = self.ident(m.rows());
        let mut u_det = self.ring.one();
        let mut pivs = vec![];

        let (mut pr, mut pc) = (0, 0);
        'pivot_loop: while pr < m.rows() {
            //find the next pivot row
            //the next pivot row is the next row with a non-zero element below the previous pivot
            'next_pivot_loop: loop {
                if pc == m.cols() {
                    break 'pivot_loop;
                }

                for r in pr..m.rows() {
                    if !self.ring.equal(m.at(r, pc).unwrap(), &self.ring.zero()) {
                        break 'next_pivot_loop;
                    }
                }

                pc += 1;
            }
            pivs.push(pc);

            if pr + 1 < m.rows() {
                //reduce everything below the pivot to zero
                for r in pr + 1..m.rows() {
                    let a = m.at(pr, pc).unwrap();
                    let b = m.at(r, pc).unwrap();
                    //if a=0 and b=0 there is nothing to do. The reduction step would fail because d=0 and we divide by d, so just skip it in this case
                    if !self.ring.equal(a, &self.ring.zero())
                        || !self.ring.equal(b, &self.ring.zero())
                    {
                        let (d, x, y) = self.ring.xgcd(a, b);
                        debug_assert!(self.ring.equal(
                            &self.ring.add(&self.ring.mul(&x, a), &self.ring.mul(&y, b)),
                            &d
                        ));
                        // perform the following row opps on self
                        // / x  -b/d \
                        // \ y   a/d /
                        let row_opp = ElementaryOpp::new_row_opp(
                            self.ring.clone(),
                            ElementaryOppType::TwoInv {
                                i: pr,
                                j: r,
                                a: x,
                                b: y,
                                //TODO: compute b/d and a/d at the same time d is computed?
                                c: self.ring.neg(&self.ring.div(b, &d).unwrap()),
                                d: self.ring.div(a, &d).unwrap(),
                            },
                        );
                        //this will implicitly put the pivot into fav assoc form because that is what gcd does
                        row_opp.apply(&mut m);
                        row_opp.apply(&mut u);
                        self.ring.mul_mut(&mut u_det, &row_opp.det());
                    }
                }
            } else {
                //explicitly put the pivot into fav assoc form
                let (unit, _assoc) = self.ring.factor_fav_assoc(m.at(pr, pc).unwrap());
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring.clone(),
                    ElementaryOppType::UnitMul {
                        row: pr,
                        unit: self.ring.inv(&unit).unwrap(),
                    },
                );
                //this will implicitly put the pivot into fav assoc form because that is what the gcd returns
                row_opp.apply(&mut m);
                row_opp.apply(&mut u);
                self.ring.mul_mut(&mut u_det, &row_opp.det());
            }

            //should have eliminated everything below the pivot
            for r in pr + 1..m.rows() {
                debug_assert!(self.ring.equal(m.at(r, pc).unwrap(), &self.ring.zero()));
            }
            pr += 1;
        }

        if m.rows() <= 4 {
            debug_assert!(self.ring.equal(&self.det_naive(&u).unwrap(), &u_det));
        }

        (m, u, u_det, pivs)
    }

    pub fn col_hermite_algorithm(
        &self,
        a: Matrix<RS::Set>,
    ) -> (Matrix<RS::Set>, Matrix<RS::Set>, RS::Set, Vec<usize>) {
        let (rh, ru, u_det, pivs) = self.row_hermite_algorithm(a.transpose());
        (rh.transpose(), ru.transpose(), u_det, pivs)
    }

    fn det_hermite(&self, a: Matrix<RS::Set>) -> RS::Set {
        let n = a.rows();
        debug_assert_eq!(n, a.cols());
        let (h, _u, u_det, _pivs) = self.row_hermite_algorithm(a);
        //h = u * self, we know det(u), and h is upper triangular
        let mut h_det = self.ring.one();
        for i in 0..n {
            self.ring.mul_mut(&mut h_det, h.at(i, i).unwrap());
        }
        self.ring.div(&h_det, &u_det).unwrap()
    }

    pub fn det(&self, a: Matrix<RS::Set>) -> Result<RS::Set, MatOppErr> {
        let n = a.rows();
        if n != a.cols() {
            Err(MatOppErr::NotSquare)
        } else if n <= 3 {
            //for speed
            Ok(self.det_naive(&a).unwrap())
        } else {
            Ok(self.det_hermite(a))
        }
    }

    pub fn rank(&self, a: Matrix<RS::Set>) -> usize {
        let (_h, _u, _u_det, pivs) = self.row_hermite_algorithm(a);
        pivs.len()
    }

    //return (u, s, v, k) such that self = usv and s is in smith normal form (with diagonal entries their favorite associates) and u, v are invertible and k is the number of non-zero elements in the diagonal of s
    pub fn smith_algorithm(
        &self,
        mut m: Matrix<RS::Set>,
    ) -> (Matrix<RS::Set>, Matrix<RS::Set>, Matrix<RS::Set>, usize) {
        let mut u = self.ident(m.rows());
        let mut v = self.ident(m.cols());

        let mut n = 0;
        'inductive_loop: while n < m.rows() && n < m.cols() {
            //search for a non-zero element to make the new starting point for (n, n)
            //having a non-zero element is necessary later in the algorithm
            'search_for_nonzero_element: {
                if self.ring.equal(m.at(n, n).unwrap(), &self.ring.zero()) {
                    //search the first row to start with
                    for c in n + 1..m.cols() {
                        if !self.ring.equal(m.at(n, c).unwrap(), &self.ring.zero()) {
                            //swap column n and column c
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring.clone(),
                                ElementaryOppType::Swap(n, c),
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut v);

                            break 'search_for_nonzero_element;
                        }
                    }
                    //search all the rows below row n
                    for r in n + 1..m.rows() {
                        for c in n..m.cols() {
                            if !self.ring.equal(m.at(r, c).unwrap(), &self.ring.zero()) {
                                //swap column n and column c
                                let col_opp = ElementaryOpp::new_col_opp(
                                    self.ring.clone(),
                                    ElementaryOppType::Swap(n, c),
                                );
                                col_opp.apply(&mut m);
                                col_opp.apply(&mut v);

                                //swap row n and row r
                                let row_opp = ElementaryOpp::new_col_opp(
                                    self.ring.clone(),
                                    ElementaryOppType::Swap(n, r),
                                );
                                row_opp.apply(&mut m);
                                row_opp.apply(&mut u);

                                break 'search_for_nonzero_element;
                            }
                        }
                    }
                    //everything remaining is zero. The algorithm can terminate here
                    break 'inductive_loop;
                }
            }

            //turn (n, n) into its favorite associate
            let (unit, _assoc) = self.ring.factor_fav_assoc(m.at(n, n).unwrap());
            let row_opp = ElementaryOpp::new_row_opp(
                self.ring.clone(),
                ElementaryOppType::UnitMul {
                    row: n,
                    unit: self.ring.inv(&unit).unwrap(),
                },
            );
            row_opp.apply(&mut m);
            row_opp.apply(&mut u);

            let mut first = true;
            let mut all_divisible;
            'zero_first_row_and_column_loop: loop {
                //replace the first row (a0, a1, ..., ak) with (gcd, 0, ..., 0). Might mess up the first column in the process
                all_divisible = true;
                for c in n + 1..m.cols() {
                    let a = m.at(n, n).unwrap();
                    let b = m.at(n, c).unwrap();
                    match self.ring.div(b, a) {
                        Ok(q) => {
                            //b is a multiple of a
                            //replace (a, b) with (a, 0) by subtracting a multiple of a from b
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring.clone(),
                                ElementaryOppType::AddRowMul {
                                    i: c,
                                    j: n,
                                    x: self.ring.neg(&q),
                                },
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut v);
                        }
                        Err(RingDivisionError::NotDivisible) => {
                            all_divisible = false;
                            //b is not a multiple of a
                            //replace (a, b) with (gcd, 0)
                            let (d, x, y) = self.ring.xgcd(a, b);
                            debug_assert!(self.ring.equal(
                                &self.ring.add(&self.ring.mul(&x, a), &self.ring.mul(&y, b)),
                                &d
                            ));
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring.clone(),
                                ElementaryOppType::TwoInv {
                                    i: n,
                                    j: c,
                                    a: x,
                                    b: y,
                                    c: self.ring.neg(&self.ring.div(b, &d).unwrap()),
                                    d: self.ring.div(a, &d).unwrap(),
                                },
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut v);
                        }
                        Err(RingDivisionError::DivideByZero) => {
                            //swap a and b
                            //a=0 so this does have the effect of (a, b) -> (gcd(a, b), 0)
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring.clone(),
                                ElementaryOppType::Swap(n, c),
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut v);
                        }
                    }
                }
                if all_divisible && !first {
                    break 'zero_first_row_and_column_loop;
                }
                first = false;

                //replace the first column (a0, a1, ..., ak) with (gcd, 0, ..., 0). Might mess up the first row in the process
                all_divisible = true;
                for r in n + 1..m.rows() {
                    let a = m.at(n, n).unwrap();
                    let b = m.at(r, n).unwrap();
                    match self.ring.div(b, a) {
                        Ok(q) => {
                            //b is a multiple of a
                            //replace (a, b) with (a, 0) by subtracting a multiple of a from b
                            let col_opp = ElementaryOpp::new_row_opp(
                                self.ring.clone(),
                                ElementaryOppType::AddRowMul {
                                    i: r,
                                    j: n,
                                    x: self.ring.neg(&q),
                                },
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut u);
                        }
                        Err(RingDivisionError::NotDivisible) => {
                            all_divisible = false;
                            //b is not a multiple of a
                            //replace (a, b) with (gcd, 0)
                            let (d, x, y) = self.ring.xgcd(a, b);
                            debug_assert!(self.ring.equal(
                                &self.ring.add(&self.ring.mul(&x, a), &self.ring.mul(&y, b)),
                                &d
                            ));
                            let row_opp = ElementaryOpp::new_row_opp(
                                self.ring.clone(),
                                ElementaryOppType::TwoInv {
                                    i: n,
                                    j: r,
                                    a: x,
                                    b: y,
                                    c: self.ring.neg(&self.ring.div(b, &d).unwrap()),
                                    d: self.ring.div(a, &d).unwrap(),
                                },
                            );
                            row_opp.apply(&mut m);
                            row_opp.apply(&mut u);
                        }
                        Err(RingDivisionError::DivideByZero) => {
                            //swap a and b
                            //a=0 so this does have the effect of (a, b) -> (gcd(a, b), 0)
                            let col_opp = ElementaryOpp::new_row_opp(
                                self.ring.clone(),
                                ElementaryOppType::Swap(n, r),
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut u);
                        }
                    }
                }
                if all_divisible {
                    break 'zero_first_row_and_column_loop;
                }
            }
            //now the first row and the first column are all zero except the top left element at (n, n) which is non-zero
            debug_assert!(!self.ring.equal(m.at(n, n).unwrap(), &self.ring.zero()));
            //some more fiddling is needed now to make sure the top left element divides everything else
            for r in n + 1..m.rows() {
                //row(n) = row(n) + row(r)
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring.clone(),
                    ElementaryOppType::AddRowMul {
                        i: n,
                        j: r,
                        x: self.ring.one(),
                    },
                );
                row_opp.apply(&mut m);
                row_opp.apply(&mut u);

                //row(n) goes from (g, a1, a2, ..., an) to (gcd, 0, 0, ..., 0) by applying col opps
                for c in n + 1..m.cols() {
                    let a = m.at(n, n).unwrap();
                    let b = m.at(n, c).unwrap();
                    // println!("a={:?} b={:?}", a, b);
                    //if a=0 and b=0 then there is nothing to do and the following step would fail when dividing by g=0
                    if !self.ring().is_zero(a) || !self.ring.is_zero(b) {
                        //b might not be a multiple of a
                        //replace (a, b) with (gcd, 0) to fix this
                        let (g, x, y) = self.ring.xgcd(a, b);
                        debug_assert!(self.ring.equal(
                            &self.ring.add(&self.ring.mul(&x, a), &self.ring.mul(&y, b)),
                            &g
                        ));
                        let col_opp = ElementaryOpp::new_col_opp(
                            self.ring.clone(),
                            ElementaryOppType::TwoInv {
                                i: n,
                                j: c,
                                a: x,
                                b: y,
                                c: self.ring.neg(&self.ring.div(b, &g).unwrap()),
                                d: self.ring.div(a, &g).unwrap(),
                            },
                        );
                        col_opp.apply(&mut m);
                        col_opp.apply(&mut v);
                    }
                }
            }

            //fix the first column
            for fix_r in n + 1..m.rows() {
                let a = m.at(n, n).unwrap();
                let b = m.at(fix_r, n).unwrap();
                let q = self.ring.div(b, a).unwrap();
                let col_opp = ElementaryOpp::new_row_opp(
                    self.ring.clone(),
                    ElementaryOppType::AddRowMul {
                        i: fix_r,
                        j: n,
                        x: self.ring.neg(&q),
                    },
                );
                col_opp.apply(&mut m);
                col_opp.apply(&mut u);
            }

            if self.ring.equal(m.at(n, n).unwrap(), &self.ring.zero()) {
                //the bottom right submatrix is all zero
                break 'inductive_loop;
            }
            n += 1;
        }

        (u, m, v, n)
    }
}

impl<RS: EuclideanDivisionStructure + BezoutDomainStructure + FavoriteAssociateStructure>
    MatrixStructure<RS>
{
    //if A:=self return (H, U, pivots) such that
    //H is in row reduced hermite normal form
    //U is invertible
    //H=UA
    //pivots[r] is the column of the rth pivot and pivots.len() == rank(A)
    pub fn row_reduced_hermite_algorithm(
        &self,
        m: Matrix<RS::Set>,
    ) -> (Matrix<RS::Set>, Matrix<RS::Set>, Vec<usize>) {
        let (mut h, mut u, _u_det, pivs) = self.row_hermite_algorithm(m);

        for (pr, pc) in pivs.iter().enumerate() {
            for r in 0..pr {
                //reduce h[r, pc] so that it has norm less than h[pr, pc]
                let a = h.at(r, *pc).unwrap();
                let b = h.at(pr, *pc).unwrap();
                //a = b*q + r
                let q = self.ring.quo(a, b).unwrap();
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring.clone(),
                    ElementaryOppType::AddRowMul {
                        i: r,
                        j: pr,
                        x: self.ring.neg(&q),
                    },
                );
                row_opp.apply(&mut h);
                row_opp.apply(&mut u);
            }
        }

        (h, u, pivs)
    }

    pub fn row_reduced_hermite_normal_form(&self, m: Matrix<RS::Set>) -> Matrix<RS::Set> {
        self.row_reduced_hermite_algorithm(m).0
    }

    pub fn col_reduced_hermite_algorithm(
        &self,
        m: Matrix<RS::Set>,
    ) -> (Matrix<RS::Set>, Matrix<RS::Set>, Vec<usize>) {
        let (rh, ru, pivs) = self.row_reduced_hermite_algorithm(m.transpose());
        (rh.transpose(), ru.transpose(), pivs)
    }

    pub fn col_reduced_hermite_normal_form(&self, m: Matrix<RS::Set>) -> Matrix<RS::Set> {
        self.col_reduced_hermite_algorithm(m).0
    }

    pub fn inv(&self, a: Matrix<RS::Set>) -> Option<Matrix<RS::Set>> {
        let n = a.rows();
        assert_eq!(n, a.cols());
        let (h, u, _pivs) = self.row_reduced_hermite_algorithm(a);
        //h = u*a
        if self.equal(&h, &self.ident(n)) {
            Some(u)
        } else {
            None
        }
    }
}

impl<RS: GreatestCommonDivisorStructure> MatrixStructure<RS> {
    pub fn factor_primitive(&self, mut mat: Matrix<RS::Set>) -> Option<(RS::Set, Matrix<RS::Set>)> {
        let entries = mat.entries_list();
        let g = self.ring.gcd_list(entries);
        if self.ring.is_zero(&g) {
            None
        } else {
            for r in 0..mat.rows() {
                for c in 0..mat.cols() {
                    *mat.at_mut(r, c).unwrap() = self.ring.div(mat.at(r, c).unwrap(), &g).unwrap();
                }
            }
            Some((g, mat))
        }
    }

    pub fn primitive_part(&self, mat: Matrix<RS::Set>) -> Option<Matrix<RS::Set>> {
        match self.factor_primitive(mat) {
            Some((_unit, prim)) => Some(prim),
            None => None,
        }
    }
}

impl<FS: FieldOfFractionsStructure> MatrixStructure<FS>
where
    FS::RS: GreatestCommonDivisorStructure,
{
    pub fn factor_primitive_fof(
        &self,
        mat: &Matrix<FS::Set>,
    ) -> (FS::Set, Matrix<<FS::RS as Structure>::Set>) {
        let div = self.ring.base_ring_structure().lcm_list(
            mat.entries_list()
                .into_iter()
                .map(|c| self.ring.denominator(&c))
                .collect(),
        );

        let (mul, prim) = MatrixStructure::new(self.ring.base_ring_structure())
            .factor_primitive(mat.apply_map(|c| {
                self.ring
                    .as_base_ring(self.ring.mul(&self.ring.from_base_ring(div.clone()), c))
                    .unwrap()
            }))
            .unwrap();

        (
            self.ring
                .div(
                    &self.ring.from_base_ring(mul),
                    &self.ring.from_base_ring(div),
                )
                .unwrap(),
            prim,
        )
    }
}

impl<FS: FieldStructure> MatrixStructure<FS> {
    pub fn presentation_matrix(
        &self,
        m: Matrix<FS::Set>,
    ) -> Result<Matrix<Polynomial<FS::Set>>, MatOppErr> {
        let n = m.rows();
        if n != m.cols() {
            Err(MatOppErr::NotSquare)
        } else {
            let poly_ring = PolynomialStructure::new(self.ring.clone());
            let poly_mat_struct = MatrixStructure::new(poly_ring.clone().into());
            Ok(poly_mat_struct
                .add(
                    &m.apply_map(|x| Polynomial::constant(x.clone())),
                    &poly_mat_struct
                        .neg(poly_mat_struct.diag(&(0..n).map(|_i| poly_ring.var()).collect())),
                )
                .unwrap())
        }
    }

    pub fn minimal_polynomial(&self, m: Matrix<FS::Set>) -> Result<Polynomial<FS::Set>, MatOppErr> {
        match self.presentation_matrix(m) {
            Ok(pres_mat) => {
                let poly_ring = PolynomialStructure::new(self.ring.clone());
                let poly_mat_struct = MatrixStructure::new(poly_ring.into());
                let (_u, s, _v, k) = poly_mat_struct.smith_algorithm(pres_mat);
                debug_assert!(k > 0); //cant be all zero becasue we are taking SNF of a non-zero matrix
                Ok(s.at(k - 1, k - 1).unwrap().clone())
            }
            Err(MatOppErr::NotSquare) => Err(MatOppErr::NotSquare),
            Err(_) => panic!(),
        }
    }

    pub fn characteristic_polynomial(
        &self,
        m: Matrix<FS::Set>,
    ) -> Result<Polynomial<FS::Set>, MatOppErr> {
        match self.presentation_matrix(m) {
            Ok(pres_mat) => {
                let poly_ring = PolynomialStructure::new(self.ring.clone());
                let poly_mat_struct = MatrixStructure::new(poly_ring.clone().into());
                let (_u, s, _v, k) = poly_mat_struct.smith_algorithm(pres_mat);
                debug_assert!(k > 0); //cant be all zero becasue we are taking SNF of a non-zero matrix
                let mut char_poly = poly_ring.one();
                for i in 0..k {
                    poly_ring.mul_mut(&mut char_poly, s.at(i, i).unwrap())
                }
                Ok(char_poly)
            }
            Err(MatOppErr::NotSquare) => Err(MatOppErr::NotSquare),
            Err(_) => panic!(),
        }
    }
}

impl<FS: ComplexConjugateStructure> MatrixStructure<FS> {
    pub fn conjugate(&self, mat: &Matrix<FS::Set>) -> Matrix<FS::Set> {
        mat.apply_map(|x| self.ring().conjugate(x))
    }

    pub fn conjugate_transpose(&self, mat: &Matrix<FS::Set>) -> Matrix<FS::Set> {
        self.conjugate(mat).transpose()
    }
}

impl<FS: ComplexConjugateStructure + RingStructure> MatrixStructure<FS> {
    // dot product of a and conj(b)
    pub fn inner_product(&self, a: &Matrix<FS::Set>, b: &Matrix<FS::Set>) -> FS::Set {
        self.dot(a, &self.conjugate(b))
    }
}

impl<FS: ComplexConjugateStructure + FieldStructure> MatrixStructure<FS> {
    //return mat=LQ where L is lower triangular and Q is row-orthogonal (not orthonormal)
    pub fn gram_schmidt_row_orthogonalization_algorithm(
        &self,
        mut mat: Matrix<FS::Set>,
    ) -> (Matrix<FS::Set>, Matrix<FS::Set>) {
        #[cfg(debug_assertions)]
        let origional_mat = mat.clone();

        let mut lt = self.ident(mat.rows());
        for i in 0..mat.rows() {
            for j in 0..i {
                //col_i = col_i - (projection of col_j onto col_i)
                let lambda = self
                    .ring()
                    .div(
                        &self.inner_product(&mat.get_row(i), &mat.get_row(j)),
                        &self.inner_product(&mat.get_row(j), &mat.get_row(j)),
                    )
                    .unwrap();
                //col_i -= lambda col_j
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring.clone(),
                    ElementaryOppType::AddRowMul {
                        i: i,
                        j: j,
                        x: self.ring().neg(&lambda),
                    },
                );
                row_opp.apply(&mut lt);
                row_opp.apply(&mut mat);
            }
        }

        for i in 0..mat.rows() {
            for j in (i + 1)..mat.rows() {
                debug_assert!(self
                    .ring()
                    .is_zero(&self.inner_product(&mat.get_row(i), &mat.get_row(j)),));
            }
        }

        #[cfg(debug_assertions)]
        assert!(self.equal(&mat, &self.mul(&lt, &origional_mat).unwrap()));

        (lt, mat)
    }

    //return mat=QR where Q is col-orthogonal (not orthonormal) and R is upper triangular and
    pub fn gram_schmidt_col_orthogonalization_algorithm(
        &self,
        mat: Matrix<FS::Set>,
    ) -> (Matrix<FS::Set>, Matrix<FS::Set>) {
        let (l, q) = self.gram_schmidt_row_orthogonalization_algorithm(mat.transpose());
        (q.transpose(), l.transpose())
    }

    pub fn gram_schmidt_row_orthogonalization(&self, mat: Matrix<FS::Set>) -> Matrix<FS::Set> {
        self.gram_schmidt_row_orthogonalization_algorithm(mat).1
    }

    pub fn gram_schmidt_col_orthogonalization(&self, mat: Matrix<FS::Set>) -> Matrix<FS::Set> {
        self.gram_schmidt_col_orthogonalization_algorithm(mat).0
    }
}

impl<FS: ComplexConjugateStructure + PositiveRealNthRootStructure + FieldStructure>
    MatrixStructure<FS>
{
    //return L*mat=Q where L is lower triangular and Q is orthonormal
    pub fn lq_decomposition_algorithm(
        &self,
        mat: Matrix<FS::Set>,
    ) -> (Matrix<FS::Set>, Matrix<FS::Set>) {
        let (mut lt, mut mat) = self.gram_schmidt_row_orthogonalization_algorithm(mat);

        for r in 0..mat.rows() {
            let row = mat.get_row(r);
            let lensq = self.inner_product(&row, &row);
            let row_opp = ElementaryOpp::new_row_opp(
                self.ring.clone(),
                ElementaryOppType::UnitMul {
                    row: r,
                    unit: self
                        .ring()
                        .inv(&self.ring().nth_root(&lensq, 2).unwrap())
                        .unwrap(),
                },
            );
            row_opp.apply(&mut lt);
            row_opp.apply(&mut mat);
        }

        debug_assert!(self.equal(
            &self.ident(mat.rows()),
            &self
                .mul(&mat, &self.conjugate(&mat.transpose_ref()))
                .unwrap()
        ));

        (lt, mat)
    }

    //return mat*R=Q where Q is col-orthogonal (not orthonormal) and R is upper triangular
    pub fn qr_decomposition_algorithm(
        &self,
        mat: Matrix<FS::Set>,
    ) -> (Matrix<FS::Set>, Matrix<FS::Set>) {
        let (l, q) = self.lq_decomposition_algorithm(mat.transpose());
        (q.transpose(), l.transpose())
    }

    pub fn gram_schmidt_row_orthonormalization(&self, mat: Matrix<FS::Set>) -> Matrix<FS::Set> {
        self.lq_decomposition_algorithm(mat).1
    }

    pub fn gram_schmidt_col_orthonormalization(&self, mat: Matrix<FS::Set>) -> Matrix<FS::Set> {
        self.qr_decomposition_algorithm(mat).0
    }
}

#[derive(Debug, Clone)]
pub struct JordanBlock<FS: AlgebraicClosureStructure>
where
    PolynomialStructure<FS>: UniqueFactorizationStructure + Structure<Set = Polynomial<FS::Set>>,
{
    eigenvalue: <FS::ACFS as Structure>::Set,
    blocksize: usize,
}

impl<FS: AlgebraicClosureStructure> JordanBlock<FS>
where
    PolynomialStructure<FS>: UniqueFactorizationStructure + Structure<Set = Polynomial<FS::Set>>,
{
    pub fn matrix(&self, field: &FS) -> Matrix<<FS::ACFS as Structure>::Set> {
        let ac_field = field.algebraic_closure_field();
        Matrix::construct(self.blocksize, self.blocksize, |r, c| {
            if r == c {
                self.eigenvalue.clone()
            } else if r + 1 == c {
                ac_field.one()
            } else {
                ac_field.zero()
            }
        })
    }
}

#[derive(Debug, Clone)]
pub struct JordanNormalForm<FS: AlgebraicClosureStructure>
where
    PolynomialStructure<FS>: UniqueFactorizationStructure + Structure<Set = Polynomial<FS::Set>>,
{
    field: Rc<FS>,
    blocks: Vec<JordanBlock<FS>>,
}

impl<FS: AlgebraicClosureStructure> JordanNormalForm<FS>
where
    PolynomialStructure<FS>: UniqueFactorizationStructure + Structure<Set = Polynomial<FS::Set>>,
{
    pub fn matrix(&self) -> Matrix<<FS::ACFS as Structure>::Set> {
        let ac_field = self.field.algebraic_closure_field();
        let ac_mat_structure = MatrixStructure::new(ac_field.clone());
        ac_mat_structure.join_diag(
            self.blocks
                .iter()
                .map(|block| block.matrix(self.field.as_ref()))
                .collect(),
        )
    }
}

impl<FS: AlgebraicClosureStructure> MatrixStructure<FS>
where
    PolynomialStructure<FS>: UniqueFactorizationStructure + Structure<Set = Polynomial<FS::Set>>,
{
    pub fn eigenvalues_list(&self, mat: Matrix<FS::Set>) -> Vec<<FS::ACFS as Structure>::Set> {
        self.ring()
            .all_roots_list(&self.characteristic_polynomial(mat).unwrap())
            .unwrap()
    }

    pub fn eigenvalues_unique(&self, mat: Matrix<FS::Set>) -> Vec<<FS::ACFS as Structure>::Set> {
        self.ring()
            .all_roots_unique(&self.characteristic_polynomial(mat).unwrap())
            .unwrap()
    }

    pub fn eigenvalues_powers(
        &self,
        mat: Matrix<FS::Set>,
    ) -> Vec<(<FS::ACFS as Structure>::Set, usize)> {
        self.ring()
            .all_roots_powers(&self.characteristic_polynomial(mat).unwrap())
            .unwrap()
    }

    pub fn generalized_col_eigenspace(
        &self,
        mat: &Matrix<FS::Set>,
        eigenvalue: &<FS::ACFS as Structure>::Set,
        k: usize,
    ) -> LinearLattice<<FS::ACFS as Structure>::Set> {
        let n = mat.rows();
        assert_eq!(n, mat.cols());
        //compute ker((M - xI)^k)
        let ac_mat_structure = MatrixStructure::new(self.ring().algebraic_closure_field());
        ac_mat_structure.col_kernel(
            ac_mat_structure.nat_pow(
                &ac_mat_structure
                    .add(
                        &mat.apply_map(|x| self.ring().algebraic_closure_inclusion(x)),
                        &ac_mat_structure.neg(
                            ac_mat_structure.mul_scalar(ac_mat_structure.ident(n), eigenvalue),
                        ),
                    )
                    .unwrap(),
                &Natural::from(k),
            ),
        )
    }

    pub fn generalized_row_eigenspace(
        &self,
        mat: &Matrix<FS::Set>,
        eigenvalue: &<FS::ACFS as Structure>::Set,
        k: usize,
    ) -> LinearLattice<<FS::ACFS as Structure>::Set> {
        LinearLatticeStructure::new(self.ring().algebraic_closure_field())
            .transpose(&self.generalized_col_eigenspace(&mat.transpose_ref(), eigenvalue, k))
    }

    pub fn col_eigenspace(
        &self,
        mat: &Matrix<FS::Set>,
        eigenvalue: &<FS::ACFS as Structure>::Set,
    ) -> LinearLattice<<FS::ACFS as Structure>::Set> {
        self.generalized_col_eigenspace(mat, eigenvalue, 1)
    }

    pub fn row_eigenspace(
        &self,
        mat: &Matrix<FS::Set>,
        eigenvalue: &<FS::ACFS as Structure>::Set,
    ) -> LinearLattice<<FS::ACFS as Structure>::Set> {
        self.generalized_row_eigenspace(mat, eigenvalue, 1)
    }

    //return the jordan normal form F of the matrix M and a basis matrix B such that
    // B^-1 M B = J
    pub fn jordan_algorithm(
        &self,
        mat: &Matrix<FS::Set>,
    ) -> (JordanNormalForm<FS>, Matrix<<FS::ACFS as Structure>::Set>) {
        let n = mat.rows();
        assert_eq!(n, mat.cols());

        let ac_field = self.ring().algebraic_closure_field();
        let ac_mat_structure = MatrixStructure::new(ac_field.clone());
        let ac_linlat_structure = LinearLatticeStructure::new(ac_field.clone());

        let ac_mat = mat.apply_map(|x| self.ring().algebraic_closure_inclusion(x));

        let mut basis = vec![];
        let mut eigenvalues = vec![]; //store (gesp_basis, eigenvalue, multiplicity)
        for (eigenvalue, multiplicity) in self.eigenvalues_powers(mat.clone()) {
            let eigenspace = self.generalized_col_eigenspace(mat, &eigenvalue, multiplicity);
            debug_assert_eq!(ac_linlat_structure.rank(&eigenspace), multiplicity);
            basis.append(&mut ac_linlat_structure.basis_matrices(&eigenspace));
            eigenvalues.push((eigenvalue, multiplicity));
        }

        //b = direct sum of generalized eigenspace
        let gesp_basis = Matrix::join_cols(mat.rows(), basis);
        //b^-1 * mat * b = block diagonal of generalized eigenspaces
        let gesp_blocks_mat = ac_mat_structure
            .mul(
                &ac_mat_structure.inv(gesp_basis.clone()).unwrap(),
                &ac_mat_structure.mul(&ac_mat, &gesp_basis).unwrap(),
            )
            .unwrap();

        let mut idx_to_block = vec![];
        for (b, (_eval, mult)) in eigenvalues.iter().enumerate() {
            for _i in 0..*mult {
                idx_to_block.push(b);
            }
        }

        // println!("{:?}", idx_to_block);
        // println!("{:?}", eigenvalues);
        // ac_mat_structure.pprint(&gesp_blocks_mat);

        //extract the blocks from the block diagonal gesp_blocks_mat
        let mut gesp_blocks = vec![];
        let mut cum_mult = 0;
        for (eval, mult) in eigenvalues {
            gesp_blocks.push((
                eval,
                mult,
                Matrix::construct(mult, mult, |r, c| {
                    gesp_blocks_mat
                        .at(cum_mult + r, cum_mult + c)
                        .unwrap()
                        .clone()
                }),
            ));
            cum_mult += mult;
        }
        debug_assert_eq!(cum_mult, n);
        drop(gesp_blocks_mat);

        // Vec<(eval, multiplicity, Vec<Jordan Block>)>
        let jnf_info = gesp_blocks
            .into_iter()
            .map(|(eval, m, mat_t)| {
                debug_assert_eq!(mat_t.rows(), m);
                debug_assert_eq!(mat_t.cols(), m);
                // println!("eval = {:?} m={}", eval, m);
                // ac_mat_structure.pprint(&mat_t);
                //all eigenvalues of T are eval
                //let S = T - x I so that all eigenvlues of S are zero
                let mat_s = ac_mat_structure
                    .add(
                        &mat_t,
                        &ac_mat_structure
                            .mul_scalar(ac_mat_structure.ident(m), &ac_field.neg(&eval)),
                    )
                    .unwrap();

                let jb_basis = {
                    debug_assert!(m >= 1);

                    let mut mat_s_pows = vec![ac_mat_structure.ident(m), mat_s.clone()];
                    for _i in 0..(m - 1) {
                        mat_s_pows.push(
                            ac_mat_structure
                                .mul(mat_s_pows.last().unwrap(), &mat_s)
                                .unwrap(),
                        );
                    }
                    debug_assert!(
                        ac_mat_structure.equal(&ac_mat_structure.zero(m, m), &mat_s_pows[m])
                    );
                    // for (i, spow) in mat_s_pows.iter().enumerate() {
                    //     println!("s^{}", i);
                    //     ac_mat_structure.pprint(&spow);
                    // }
                    let mat_s_pow_kers = mat_s_pows
                        .into_iter()
                        .map(|s_mat_pow| ac_mat_structure.col_kernel(s_mat_pow))
                        .collect_vec();
                    // ker(S) in ker(S^2) in ker(S^3) in ...
                    // for (i, ker) in mat_s_pow_kers.iter().enumerate() {
                    //     println!("ker(s^{})", i);
                    //     ac_linlat_structure.pprint(&ker);
                    // }

                    let mut accounted = ac_linlat_structure.zero(m, 1);
                    let mut jordan_block_bases = vec![];
                    for k in (0..m).rev() {
                        //extend the basis by stuff in ker(S^{k+1}) but not in ker(S^k) and their images under S, and which are not already accounted for
                        // println!("k = {} {}", k + 1, k);
                        let ker_ext = ac_linlat_structure.from_basis(
                            m,
                            1,
                            ac_linlat_structure
                                .extension_basis(&mat_s_pow_kers[k], &mat_s_pow_kers[k + 1]),
                        );

                        let unaccounted_ker_ext_basis = ac_linlat_structure.extension_basis(
                            &ac_linlat_structure.intersect_pair(m, 1, &accounted, &ker_ext),
                            &ker_ext,
                        );
                        // let unaccounted_ker_ext =
                        //     ac_linlat_structure.from_basis(m, 1, unaccounted_ker_ext_basis.clone());

                        for ukeb in unaccounted_ker_ext_basis {
                            //one new jordan block for each ukeb
                            // println!("ukeb");
                            // ac_mat_structure.pprint(&ukeb);

                            let mut jb_basis = vec![ukeb];
                            for _i in 0..k {
                                let ukeb_img = ac_mat_structure
                                    .mul(&mat_s, &jb_basis.last().unwrap())
                                    .unwrap();
                                // println!("ukeb_img #{}", i);
                                // ac_mat_structure.pprint(&ukeb_img);
                                jb_basis.push(ukeb_img);
                            }

                            accounted = ac_linlat_structure.sum_pair(
                                m,
                                1,
                                &accounted,
                                &ac_linlat_structure.from_basis(m, 1, jb_basis.clone()),
                            );

                            jordan_block_bases.push(jb_basis.into_iter().rev().collect_vec());
                        }
                    }

                    // println!(
                    //     "jb sizes = {:?}",
                    //     jordan_block_bases.iter().map(|v| v.len()).collect_vec()
                    // );

                    jordan_block_bases

                    // Matrix::join_cols(m, jordan_block_bases.into_iter().flatten().collect_vec())
                };

                // println!("gesp_jordan_basis");
                // ac_mat_structure.pprint(&jb_basis);
                // ac_mat_structure.pprint(
                //     &ac_mat_structure
                //         .mul(
                //             &ac_mat_structure.inv(jb_basis.clone()).unwrap(),
                //             &ac_mat_structure.mul(&mat_t, &jb_basis).unwrap(),
                //         )
                //         .unwrap(),
                // );

                (eval, m, jb_basis)
            })
            .collect_vec();

        let mut jordan_blocks = vec![];
        let mut jnf_basis_rel_gesp_basis: Vec<Matrix<<FS::ACFS as Structure>::Set>> = vec![];
        for (eval, mult, blocks) in jnf_info {
            // println!("eval={:?}, mult={}", eval, mult);
            let mut eigenblock_basis = vec![];
            for mut block in blocks {
                jordan_blocks.push(JordanBlock {
                    eigenvalue: eval.clone(),
                    blocksize: block.len(),
                });
                eigenblock_basis.append(&mut block);
            }
            jnf_basis_rel_gesp_basis.push(Matrix::join_cols(mult, eigenblock_basis));
        }
        let jnf = JordanNormalForm {
            field: self.ring(),
            blocks: jordan_blocks,
        };
        let jordan_blocks_basis = ac_mat_structure.join_diag(jnf_basis_rel_gesp_basis);

        // ac_mat_structure.pprint(&jnf.matrix());

        // println!("jordan_blocks_basis");
        // ac_mat_structure.pprint(&jordan_blocks_basis);

        let jnf_basis = ac_mat_structure
            .mul(&gesp_basis, &jordan_blocks_basis)
            .unwrap();
        // println!("jnf_basis");
        // ac_mat_structure.pprint(&jnf_basis);

        //check that B^-1 M B = JNF
        debug_assert!(ac_mat_structure.equal(
            &ac_mat_structure
                .mul(
                    &ac_mat_structure.inv(jnf_basis.clone()).unwrap(),
                    &ac_mat_structure.mul(&ac_mat, &jnf_basis).unwrap(),
                )
                .unwrap(),
            &jnf.matrix()
        ));
        // println!("jnf");
        // ac_mat_structure.pprint(&jnf);

        // todo!()

        (jnf, jnf_basis)
    }

    pub fn jordan_normal_form(
        &self,
        mat: &Matrix<FS::Set>,
    ) -> Matrix<<FS::ACFS as Structure>::Set> {
        self.jordan_algorithm(mat).0.matrix()
    }

    //TODO: find basis which make two matricies similar if one exists
    /*
    def similar_basis(self, other):
    #find a basis in which self looks like other
    #equivelently, find P such that P^-1*self*P == other
    if type(self) == type(other) == Matrix:
        if self.n == other.n:
            if self.jcf_spec() == other.jcf_spec():
                #need to find a jcf basis for self and other such that the jcf matricies are identical (including order (thats the only hard part))
                #NOTE - by the implementation of the algorithm used, each eigen block will be consistently ordered - largest first
                #HOWEVER, the order of the eigen block is still unknown (and an order cant be imposed in the algorhtm becasue in general, arbitrary number things cant be ordered in a consistent way)
                self_jcf_info = self.jcf_info()
                other_jcf_info = other.jcf_info()
                #rewrite these in terms of {e_val : info}
                self_jcf_info = {info["ev"] : info for info in self_jcf_info}
                other_jcf_info = {info["ev"] : info for info in other_jcf_info}
                assert self_jcf_info.keys() == other_jcf_info.keys()
                keys = list(self_jcf_info.keys()) #decide a consistent order here
                #reorder the info
                self_jcf_info = [self_jcf_info[ev] for ev in keys]
                other_jcf_info = [other_jcf_info[ev] for ev in keys]
                #now both info lists have the eigen values in the same order
                #as well as all blocks within each eigen block being in the right order
                self_jcf_basis = Matrix.jcf_info_to_jcf_basis(self_jcf_info)
                other_jcf_basis = Matrix.jcf_info_to_jcf_basis(other_jcf_info)
                return self_jcf_basis * other_jcf_basis ** -1
            else:
                raise Exception("Matricies are not similar so cant find a basis in which one looks like the other")
    raise NotImplementedError
    */
}

impl<R: StructuredType> StructuredType for Matrix<R>
where
    R::Structure: Structure,
{
    type Structure = MatrixStructure<R::Structure>;

    fn structure() -> Rc<Self::Structure> {
        MatrixStructure::new(R::structure()).into()
    }
}

impl<R: StructuredType> Matrix<R>
where
    R::Structure: DisplayableStructure,
{
    pub fn pprint(&self) {
        Self::structure().pprint(self)
    }
}

impl<R: StructuredType> PartialEq for Matrix<R>
where
    R::Structure: RingStructure,
{
    fn eq(&self, other: &Self) -> bool {
        Self::structure().equal(self, other)
    }
}

impl<R: StructuredType> Eq for Matrix<R> where R::Structure: RingStructure {}

impl<R: StructuredType> Matrix<R>
where
    R::Structure: RingStructure,
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
        Self::structure().neg_mut(self)
    }

    pub fn neg(&self) -> Self {
        Self::structure().neg(self.clone())
    }

    pub fn mul(a: &Self, b: &Self) -> Result<Self, MatOppErr> {
        Self::structure().mul(a, b)
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
        Self::structure().trace(&self)
    }
}

impl<R: StructuredType> Matrix<R>
where
    R::Structure: BezoutDomainStructure,
{
    pub fn row_span(&self) -> LinearLattice<R> {
        Self::structure().row_span(self.clone())
    }

    pub fn col_span(&self) -> LinearLattice<R> {
        Self::structure().col_span(self.clone())
    }

    pub fn row_affine_span(&self) -> AffineLattice<R> {
        Self::structure().row_affine_span(self.clone())
    }

    pub fn col_affine_span(&self) -> AffineLattice<R> {
        Self::structure().col_affine_span(self.clone())
    }

    pub fn row_kernel(&self) -> LinearLattice<R> {
        Self::structure().row_kernel(self.clone())
    }

    pub fn col_kernel(&self) -> LinearLattice<R> {
        Self::structure().col_kernel(self.clone())
    }

    pub fn row_solve(&self, y: impl Borrow<Self>) -> Option<Self> {
        Self::structure().row_solve(self, y)
    }

    pub fn col_solve(&self, y: impl Borrow<Self>) -> Option<Self> {
        Self::structure().col_solve(self, y)
    }

    pub fn row_solution_lattice(&self, y: impl Borrow<Self>) -> AffineLattice<R> {
        Self::structure().row_solution_lattice(self, y)
    }

    pub fn col_solution_lattice(&self, y: impl Borrow<Self>) -> AffineLattice<R> {
        Self::structure().col_solution_lattice(self, y)
    }

    pub fn row_hermite_algorithm(&self) -> (Self, Self, R, Vec<usize>) {
        Self::structure().row_hermite_algorithm(self.clone())
    }

    pub fn col_hermite_algorithm(&self) -> (Self, Self, R, Vec<usize>) {
        Self::structure().col_hermite_algorithm(self.clone())
    }

    pub fn det(&self) -> Result<R, MatOppErr> {
        Self::structure().det(self.clone())
    }

    pub fn rank(&self) -> usize {
        Self::structure().rank(self.clone())
    }

    pub fn smith_algorithm(&self) -> (Self, Self, Self, usize) {
        Self::structure().smith_algorithm(self.clone())
    }
}

impl<R: StructuredType> Matrix<R>
where
    R::Structure: EuclideanDivisionStructure + BezoutDomainStructure + FavoriteAssociateStructure,
{
    pub fn row_reduced_hermite_algorithm(&self) -> (Self, Self, Vec<usize>) {
        Self::structure().row_reduced_hermite_algorithm(self.clone())
    }

    pub fn row_reduced_hermite_normal_form(&self) -> Self {
        Self::structure().row_reduced_hermite_normal_form(self.clone())
    }

    pub fn col_reduced_hermite_algorithm(&self) -> (Self, Self, Vec<usize>) {
        Self::structure().col_reduced_hermite_algorithm(self.clone())
    }

    pub fn col_reduced_hermite_normal_form(&self) -> Self {
        Self::structure().col_reduced_hermite_normal_form(self.clone())
    }

    pub fn inv(&self) -> Option<Matrix<R>> {
        Self::structure().inv(self.clone())
    }
}

impl<R: StructuredType> Matrix<R>
where
    R::Structure: GreatestCommonDivisorStructure,
{
    pub fn factor_primitive(self) -> Option<(R, Matrix<R>)> {
        Self::structure().factor_primitive(self)
    }

    pub fn primitive_part(self) -> Option<Matrix<R>> {
        Self::structure().primitive_part(self)
    }
}

impl<F: StructuredType> Matrix<F>
where
    F::Structure: FieldOfFractionsStructure,
    <F::Structure as FieldOfFractionsStructure>::RS: GreatestCommonDivisorStructure,
{
    pub fn factor_primitive_fof(
        &self,
    ) -> (
        F,
        Matrix<<<F::Structure as FieldOfFractionsStructure>::RS as Structure>::Set>,
    ) {
        Self::structure().factor_primitive_fof(self)
    }
}

impl<F: StructuredType> Matrix<F>
where
    F::Structure: FieldStructure,
{
    pub fn presentation_matrix(&self) -> Result<Matrix<Polynomial<F>>, MatOppErr> {
        Self::structure().presentation_matrix(self.clone())
    }

    pub fn minimal_polynomial(&self) -> Result<Polynomial<F>, MatOppErr> {
        Self::structure().minimal_polynomial(self.clone())
    }

    pub fn characteristic_polynomial(&self) -> Result<Polynomial<F>, MatOppErr> {
        Self::structure().characteristic_polynomial(self.clone())
    }
}

impl<F: StructuredType> Matrix<F>
where
    F::Structure: ComplexConjugateStructure + FieldStructure,
{
    pub fn gram_schmidt_row_orthogonalization_algorithm(self) -> (Matrix<F>, Matrix<F>) {
        Self::structure().gram_schmidt_row_orthogonalization_algorithm(self)
    }

    pub fn gram_schmidt_col_orthogonalization_algorithm(self) -> (Matrix<F>, Matrix<F>) {
        Self::structure().gram_schmidt_col_orthogonalization_algorithm(self)
    }

    pub fn gram_schmidt_row_orthogonalization(self) -> Matrix<F> {
        Self::structure().gram_schmidt_row_orthogonalization(self)
    }

    pub fn gram_schmidt_col_orthogonalization(self) -> Matrix<F> {
        Self::structure().gram_schmidt_col_orthogonalization(self)
    }
}

impl<F: StructuredType> Matrix<F>
where
    F::Structure: ComplexConjugateStructure + PositiveRealNthRootStructure + FieldStructure,
{
    pub fn lq_decomposition_algorithm(self) -> (Matrix<F>, Matrix<F>) {
        Self::structure().lq_decomposition_algorithm(self)
    }

    pub fn qr_decomposition_algorithm(self) -> (Matrix<F>, Matrix<F>) {
        Self::structure().qr_decomposition_algorithm(self)
    }

    pub fn gram_schmidt_row_orthonormalization(self) -> Matrix<F> {
        Self::structure().gram_schmidt_row_orthonormalization(self)
    }

    pub fn gram_schmidt_col_orthonormalization(self) -> Matrix<F> {
        Self::structure().gram_schmidt_col_orthonormalization(self)
    }
}

impl<F: StructuredType> Matrix<F>
where
    F::Structure: AlgebraicClosureStructure,
    PolynomialStructure<F::Structure>:
        UniqueFactorizationStructure + Structure<Set = Polynomial<F>>,
{
    pub fn eigenvalues_list(
        self,
    ) -> Vec<<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set> {
        Self::structure().eigenvalues_list(self)
    }

    pub fn eigenvalues_unique(
        self,
    ) -> Vec<<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set> {
        Self::structure().eigenvalues_unique(self)
    }

    pub fn eigenvalues_powers(
        self,
    ) -> Vec<(
        <<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set,
        usize,
    )> {
        Self::structure().eigenvalues_powers(self)
    }

    pub fn generalized_col_eigenspace(
        &self,
        eigenvalue: &<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set,
        k: usize,
    ) -> LinearLattice<<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set> {
        Self::structure().generalized_col_eigenspace(self, eigenvalue, k)
    }

    pub fn generalized_row_eigenspace(
        &self,
        eigenvalue: &<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set,
        k: usize,
    ) -> LinearLattice<<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set> {
        Self::structure().generalized_row_eigenspace(self, eigenvalue, k)
    }

    pub fn col_eigenspace(
        &self,
        eigenvalue: &<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set,
    ) -> LinearLattice<<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set> {
        Self::structure().col_eigenspace(self, eigenvalue)
    }

    pub fn row_eigenspace(
        &self,
        eigenvalue: &<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set,
    ) -> LinearLattice<<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set> {
        Self::structure().row_eigenspace(self, eigenvalue)
    }

    //return the jordan normal form F of the matrix M and a basis matrix B such that
    // B^-1 M B = J
    pub fn jordan_algorithm(
        &self,
    ) -> (
        JordanNormalForm<F::Structure>,
        Matrix<<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set>,
    ) {
        Self::structure().jordan_algorithm(self)
    }

    pub fn jordan_normal_form(
        &self,
    ) -> Matrix<<<F::Structure as AlgebraicClosureStructure>::ACFS as Structure>::Set> {
        Self::structure().jordan_normal_form(self)
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use malachite_nz::integer::Integer;
    use malachite_q::Rational;

    use crate::{
        number::algebraic::isolated_roots::{complex::ComplexAlgebraic, real::RealAlgebraic},
        ring_structure::cannonical::{EuclideanDivisionDomain, FavoriteAssociateDomain, Ring},
    };

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
        match m.check_invariants() {
            Ok(()) => panic!(),
            Err(_) => {}
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

    #[test]
    fn hermite_algorithm() {
        for a in vec![
            Matrix::from_rows(vec![
                vec![
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(6),
                    Integer::from(2),
                ],
                vec![
                    Integer::from(5),
                    Integer::from(6),
                    Integer::from(1),
                    Integer::from(6),
                ],
                vec![
                    Integer::from(8),
                    Integer::from(3),
                    Integer::from(1),
                    Integer::from(1),
                ],
            ]),
            Matrix::from_rows(vec![
                vec![
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(6),
                    Integer::from(2),
                ],
                vec![
                    Integer::from(5),
                    Integer::from(6),
                    Integer::from(1),
                    Integer::from(6),
                ],
                vec![
                    Integer::from(8),
                    Integer::from(3),
                    Integer::from(1),
                    Integer::from(1),
                ],
            ])
            .transpose(),
            Matrix::from_rows(vec![
                vec![
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(5),
                    Integer::from(2),
                ],
                vec![
                    Integer::from(5),
                    Integer::from(6),
                    Integer::from(11),
                    Integer::from(6),
                ],
                vec![
                    Integer::from(8),
                    Integer::from(3),
                    Integer::from(11),
                    Integer::from(1),
                ],
            ]),
            Matrix::from_rows(vec![
                vec![Integer::from(0), Integer::from(4), Integer::from(4)],
                vec![Integer::from(1), Integer::from(6), Integer::from(12)],
                vec![Integer::from(1), Integer::from(4), Integer::from(16)],
            ]),
        ] {
            println!();
            println!("hermite reduced row algorithm for");
            a.pprint();
            let (h, u, pivs) = a.clone().row_reduced_hermite_algorithm();
            println!("H =");
            h.pprint();
            println!("pivs = {:?}", pivs);
            println!("U =");
            u.pprint();
            assert_eq!(h, Matrix::mul(&u, &a).unwrap());

            //trace the boundary of zeros and check that everything under is zero
            let mut rz = 0;
            for cz in 0..h.cols() {
                match pivs.get(rz) {
                    Some(cp) => {
                        if cp == &cz {
                            rz += 1;
                        }
                    }
                    None => {}
                }
                for r in rz..h.rows() {
                    assert_eq!(h.at(r, cz).unwrap(), &Integer::zero());
                }
            }

            //check pivot columns
            for (pr, pc) in pivs.iter().enumerate() {
                assert!(h.at(pr, *pc).unwrap() != &Integer::zero());
                for r in 0..h.rows() {
                    if r > pr {
                        assert_eq!(h.at(r, *pc).unwrap(), &Integer::zero());
                    } else if r == pr {
                        let (_unit, assoc) = Integer::factor_fav_assoc(h.at(r, *pc).unwrap());
                        assert_eq!(&assoc, h.at(r, *pc).unwrap());
                    } else {
                        assert!(
                            Integer::norm(h.at(r, *pc).unwrap())
                                < Integer::norm(h.at(pr, *pc).unwrap())
                        );
                    }
                }
            }

            println!();
            println!("hermite reduced col algorithm for");
            a.pprint();
            let (h, u, pivs) = a.clone().col_reduced_hermite_algorithm();
            println!("H =");
            h.pprint();
            println!("pivs = {:?}", pivs);
            println!("U =");
            u.pprint();

            //trace the boundary of zeros and check that everything to the right is zero
            let mut cz = 0;
            for rz in 0..h.rows() {
                match pivs.get(cz) {
                    Some(rp) => {
                        if rp == &rz {
                            cz += 1;
                        }
                    }
                    None => {}
                }
                for c in cz..h.cols() {
                    assert_eq!(h.at(rz, c).unwrap(), &Integer::zero());
                }
            }

            //check the pivot rows
            assert_eq!(h, Matrix::mul(&a, &u).unwrap());
            for (pc, pr) in pivs.iter().enumerate() {
                assert!(h.at(*pr, pc).unwrap() != &Integer::zero());
                for c in 0..h.cols() {
                    if c > pc {
                        assert_eq!(h.at(*pr, c).unwrap(), &Integer::zero());
                    } else if c == pc {
                        let (_unit, assoc) = Integer::factor_fav_assoc(h.at(*pr, c).unwrap());
                        assert_eq!(&assoc, h.at(*pr, c).unwrap());
                    } else {
                        assert!(
                            Integer::norm(h.at(*pr, c).unwrap())
                                < Integer::norm(h.at(*pr, pc).unwrap())
                        );
                    }
                }
            }
        }

        {
            //integer reduced hermite normal form is unique, so we can fully check an example
            let a = Matrix::<Integer>::from_rows(vec![
                vec![
                    Integer::from(2),
                    Integer::from(3),
                    Integer::from(6),
                    Integer::from(2),
                ],
                vec![
                    Integer::from(5),
                    Integer::from(6),
                    Integer::from(1),
                    Integer::from(6),
                ],
                vec![
                    Integer::from(8),
                    Integer::from(3),
                    Integer::from(1),
                    Integer::from(1),
                ],
            ]);

            let expected_h = Matrix::from_rows(vec![
                vec![
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(50),
                    Integer::from(-11),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(3),
                    Integer::from(28),
                    Integer::from(-2),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(61),
                    Integer::from(-13),
                ],
            ]);

            let expected_u = Matrix::from_rows(vec![
                vec![Integer::from(9), Integer::from(-5), Integer::from(1)],
                vec![Integer::from(5), Integer::from(-2), Integer::from(0)],
                vec![Integer::from(11), Integer::from(-6), Integer::from(1)],
            ]);

            let (h, u, pivs) = a.clone().row_reduced_hermite_algorithm();

            assert_eq!(h, expected_h);
            assert_eq!(u, expected_u);
            assert_eq!(pivs, vec![0, 1, 2]);
        }

        {
            //this one used to cause a dividion by zero error when replacing (a, b) with (gcd, 0) when (a, b) = (0, 0)
            let a = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(1)],
            ]);

            let expected_h = Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(1)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
            ]);

            let expected_u = Matrix::from_rows(vec![
                vec![
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(-1),
                    Integer::from(0),
                    Integer::from(0),
                ],
            ]);

            let (h, u, pivs) = a.clone().row_reduced_hermite_algorithm();

            assert_eq!(h, expected_h);
            assert_eq!(u, expected_u);
            assert_eq!(pivs, vec![0, 1, 2]);
        }
    }

    #[test]
    fn invert() {
        let a = Matrix::<Rational>::from_rows(vec![
            vec![Rational::from(2), Rational::from(4), Rational::from(4)],
            vec![Rational::from(-6), Rational::from(6), Rational::from(12)],
            vec![Rational::from(10), Rational::from(7), Rational::from(17)],
        ]);
        a.pprint();
        println!("{:?}", a.rank());
        let b = Matrix::inv(&a).unwrap();
        b.pprint();
        debug_assert_eq!(Matrix::mul(&a, &b).unwrap(), Matrix::ident(3));
        debug_assert_eq!(Matrix::mul(&b, &a).unwrap(), Matrix::ident(3));
    }

    #[test]
    fn smith_algorithm() {
        {
            let a = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(2), Integer::from(4), Integer::from(4)],
                vec![Integer::from(-6), Integer::from(6), Integer::from(12)],
                vec![Integer::from(10), Integer::from(4), Integer::from(16)],
            ]);
            let (u, s, v, k) = a.clone().smith_algorithm();
            assert_eq!(s, Matrix::mul(&Matrix::mul(&u, &a).unwrap(), &v).unwrap());
            assert_eq!(k, 3);
            assert_eq!(
                s,
                Matrix::from_rows(vec![
                    vec![Integer::from(2), Integer::from(0), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(2), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(0), Integer::from(156)]
                ])
            );

            let a = Matrix::<Integer>::from_rows(vec![
                vec![
                    Integer::from(-6),
                    Integer::from(111),
                    Integer::from(-36),
                    Integer::from(6),
                ],
                vec![
                    Integer::from(5),
                    Integer::from(-672),
                    Integer::from(210),
                    Integer::from(74),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(-255),
                    Integer::from(81),
                    Integer::from(24),
                ],
                vec![
                    Integer::from(-7),
                    Integer::from(255),
                    Integer::from(-81),
                    Integer::from(-10),
                ],
            ]);
            let (u, s, v, k) = a.clone().smith_algorithm();
            assert_eq!(s, Matrix::mul(&Matrix::mul(&u, &a).unwrap(), &v).unwrap());
            assert_eq!(k, 3);
            assert_eq!(
                s,
                Matrix::from_rows(vec![
                    vec![
                        Integer::from(1),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(3),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(21),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ]
                ])
            );
        }

        {
            //used to cause a divide by zero
            let a = Matrix::<Integer>::from_rows(vec![
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                ],
            ]);
            let (_u, _s, _v, k) = a.clone().smith_algorithm();
            assert_eq!(k, 1);
        }
    }
    #[test]
    fn min_and_char_polys() {
        {
            let a = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(0), Integer::from(4), Integer::from(4)],
                vec![Integer::from(1), Integer::from(4), Integer::from(16)],
            ]);
            let (_u, s, _v, _k) = a.smith_algorithm();

            assert_eq!(
                s,
                Matrix::from_rows(vec![
                    vec![Integer::from(1), Integer::from(0), Integer::from(0),],
                    vec![Integer::from(0), Integer::from(4), Integer::from(0)],
                ])
            );
        }

        {
            let a = Matrix::<Rational>::from_rows(vec![
                vec![Rational::from(0), Rational::from(0), Rational::from(0)],
                vec![Rational::from(0), Rational::from(0), Rational::from(1)],
                vec![Rational::from(0), Rational::from(0), Rational::from(0)],
            ]);
            let min_p = a.clone().minimal_polynomial().unwrap();
            let char_p = a.clone().characteristic_polynomial().unwrap();
            assert_eq!(
                &min_p,
                &Polynomial::from_coeffs(vec![
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(1)
                ])
            );
            assert_eq!(
                &char_p,
                &Polynomial::from_coeffs(vec![
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(1)
                ])
            );
        }

        {
            let a = Matrix::<Rational>::from_rows(vec![
                vec![
                    Rational::from(4),
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(0),
                ],
                vec![
                    Rational::from(0),
                    Rational::from(4),
                    Rational::from(0),
                    Rational::from(0),
                ],
                vec![
                    Rational::from(0),
                    Rational::from(1),
                    Rational::from(4),
                    Rational::from(0),
                ],
                vec![
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(1),
                    Rational::from(4),
                ],
            ]);
            let min_p = a.clone().minimal_polynomial().unwrap();
            let char_p = a.clone().characteristic_polynomial().unwrap();
            assert_eq!(
                &min_p,
                &Polynomial::from_coeffs(vec![Rational::from(-4), Rational::from(1),])
                    .nat_pow(&Natural::from(3u8))
            );
            assert_eq!(
                &char_p,
                &Polynomial::from_coeffs(vec![Rational::from(-4), Rational::from(1),])
                    .nat_pow(&Natural::from(4u8))
            );
        }
    }

    #[test]
    fn span_and_kernel_rank() {
        let mat = Matrix::<Integer>::from_rows(vec![
            vec![Integer::from(1), Integer::from(1), Integer::from(1)],
            vec![Integer::from(1), Integer::from(2), Integer::from(1)],
            vec![Integer::from(1), Integer::from(1), Integer::from(1)],
            vec![Integer::from(1), Integer::from(1), Integer::from(1)],
        ]);

        assert_eq!(mat.clone().row_span().rank(), 2);
        assert_eq!(mat.clone().col_span().rank(), 2);
        assert_eq!(mat.clone().row_kernel().rank(), 2);
        assert_eq!(mat.clone().col_kernel().rank(), 1);
    }

    #[test]
    fn affine_span() {
        {
            //row affine span
            let lat1 = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(1), Integer::from(1)],
                vec![Integer::from(3), Integer::from(1)],
                vec![Integer::from(2), Integer::from(3)],
            ])
            .row_affine_span();

            let lat2 = AffineLattice::from_offset_and_linear_lattice(
                1,
                2,
                Matrix::from_rows(vec![vec![Integer::from(2), Integer::from(3)]]),
                Matrix::from_rows(vec![
                    vec![Integer::from(1), Integer::from(2)],
                    vec![Integer::from(-1), Integer::from(2)],
                ])
                .row_span(),
            );

            lat1.pprint();
            lat2.pprint();

            assert_eq!(lat1, lat2);
        }

        {
            //column affine span
            let lat1 = Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(1), Integer::from(3), Integer::from(2)],
                vec![Integer::from(1), Integer::from(1), Integer::from(3)],
            ])
            .col_affine_span();

            let lat2 = AffineLattice::from_offset_and_linear_lattice(
                2,
                1,
                Matrix::from_rows(vec![vec![Integer::from(2)], vec![Integer::from(3)]]),
                Matrix::from_rows(vec![
                    vec![Integer::from(1), Integer::from(-1)],
                    vec![Integer::from(2), Integer::from(2)],
                ])
                .col_span(),
            );

            lat1.pprint();
            lat2.pprint();

            assert_eq!(&lat1, &lat2);
        }
    }

    #[test]
    fn span_and_kernel_points() {
        let mat = Matrix::<Integer>::from_rows(vec![
            vec![
                Integer::from(1),
                Integer::from(1),
                Integer::from(0),
                Integer::from(0),
            ],
            vec![
                Integer::from(1),
                Integer::from(1),
                Integer::from(0),
                Integer::from(0),
            ],
            vec![
                Integer::from(0),
                Integer::from(0),
                Integer::from(3),
                Integer::from(5),
            ],
        ]);

        println!("matrix");
        mat.pprint();

        let k = mat.col_kernel();
        println!("kernel");
        k.pprint();

        assert!(k.contains_point(Matrix::from_rows(vec![
            vec![Integer::from(-1)],
            vec![Integer::from(1)],
            vec![Integer::from(5)],
            vec![Integer::from(-3)]
        ])));
    }

    #[test]
    fn rational_gram_schmidt() {
        let mat = Matrix::<Rational>::from_rows(vec![
            vec![Rational::from(1), Rational::from(-1), Rational::from(3)],
            vec![Rational::from(1), Rational::from(0), Rational::from(5)],
            vec![Rational::from(1), Rational::from(2), Rational::from(6)],
        ]);

        let mat_expected_gs = Matrix::from_rows(vec![
            vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("-4/3").unwrap(),
                Rational::from_str("-3/7").unwrap(),
            ],
            vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("-1/3").unwrap(),
                Rational::from_str("9/14").unwrap(),
            ],
            vec![
                Rational::from_str("1").unwrap(),
                Rational::from_str("5/3").unwrap(),
                Rational::from_str("-3/14").unwrap(),
            ],
        ]);

        mat.pprint();
        mat_expected_gs.pprint();

        let (mat_actual_gs, mat_actual_ut) =
            mat.clone().gram_schmidt_col_orthogonalization_algorithm();

        mat_actual_gs.pprint();
        mat_actual_ut.pprint();
        Matrix::mul(&mat, &mat_actual_ut).unwrap().pprint();

        assert_eq!(mat.gram_schmidt_col_orthogonalization(), mat_expected_gs);
    }

    #[test]
    fn complex_gram_schmidt() {
        let i = &ComplexAlgebraic::i().into_ring();

        let mat = Matrix::<ComplexAlgebraic>::from_rows(vec![
            vec![(1 + 0 * i).into_set(), (1 * i).into_set()],
            vec![(1 + 0 * i).into_set(), (1 + 0 * i).into_set()],
        ]);
        mat.pprint();
        mat.gram_schmidt_col_orthogonalization().pprint();

        let mat = Matrix::<ComplexAlgebraic>::from_rows(vec![
            vec![
                (-2 + 2 * i).into_set(),
                (7 + 3 * i).into_set(),
                (7 + 3 * i).into_set(),
            ],
            vec![
                (3 + 3 * i).into_set(),
                (-2 + 4 * i).into_set(),
                (6 + 2 * i).into_set(),
            ],
            vec![
                (2 + 2 * i).into_set(),
                (8 + 0 * i).into_set(),
                (-9 + 1 * i).into_set(),
            ],
        ]);
        mat.pprint();
        mat.clone().gram_schmidt_col_orthogonalization().pprint();
    }

    #[test]
    fn complex_gram_schmidt_normalized() {
        let one = &RealAlgebraic::one().into_ring();

        let mat = Matrix::<RealAlgebraic>::from_rows(vec![
            vec![(1 * one).into_set(), (1 * one).into_set()],
            vec![(1 * one).into_set(), (2 * one).into_set()],
        ]);
        mat.pprint();
        mat.gram_schmidt_col_orthonormalization().pprint();

        let i = &ComplexAlgebraic::i().into_ring();
        let mat = Matrix::<ComplexAlgebraic>::from_rows(vec![
            vec![(-2 + 2 * i).into_set(), (-9 + 1 * i).into_set()],
            vec![(3 + 3 * i).into_set(), (-2 + 4 * i).into_set()],
        ]);
        mat.pprint();
        mat.clone().gram_schmidt_col_orthonormalization().pprint();
    }

    #[test]
    fn jordan_normal_form() {
        let mat = Matrix::from_rows(vec![
            vec![
                Rational::from(3),
                Rational::from(1),
                Rational::from(0),
                Rational::from(1),
            ],
            vec![
                Rational::from(-1),
                Rational::from(5),
                Rational::from(4),
                Rational::from(1),
            ],
            vec![
                Rational::from(0),
                Rational::from(0),
                Rational::from(2),
                Rational::from(0),
            ],
            vec![
                Rational::from(0),
                Rational::from(0),
                Rational::from(0),
                Rational::from(4),
            ],
        ]);

        mat.pprint();
        for root in MatrixStructure::new(Rational::structure()).eigenvalues_list(mat.clone()) {
            println!("{}", root);
        }

        let (j, b) = MatrixStructure::new(Rational::structure()).jordan_algorithm(&mat);
        println!("{:?}", j);
        j.matrix().pprint();
        b.pprint();

        let mat = Matrix::<Rational>::from_rows(vec![vec![1, 0, 0], vec![0, 0, -1], vec![0, 2, 0]]);

        mat.pprint();
        for root in MatrixStructure::new(Rational::structure()).eigenvalues_list(mat.clone()) {
            println!("{}", root);
        }

        let (j, b) = MatrixStructure::new(Rational::structure()).jordan_algorithm(&mat);
        println!("{:?}", j);
        j.matrix().pprint();
        b.pprint();
    }
}

/*
//for LLL testing from wikipeida
let mat = Matrix::<ComplexAlgebraic>::from_rows(vec![
            vec![
                (-2 + 2 * i).into_set(),
                (7 + 3 * i).into_set(),
                (7 + 3 * i).into_set(),
                (-5 + 4 * i).into_set(),
            ],
            vec![
                (3 + 3 * i).into_set(),
                (-2 + 4 * i).into_set(),
                (6 + 2 * i).into_set(),
                (-1 + 4 * i).into_set(),
            ],
            vec![
                (2 + 2 * i).into_set(),
                (8 + 0 * i).into_set(),
                (-9 + 1 * i).into_set(),
                (-7 + 5 * i).into_set(),
            ],
            vec![
                (8 + 2 * i).into_set(),
                (-9 + 0 * i).into_set(),
                (6 + 3 * i).into_set(),
                (-4 + 4 * i).into_set(),
            ],
        ]);
        */
