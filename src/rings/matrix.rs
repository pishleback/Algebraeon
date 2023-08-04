#![allow(dead_code)]

use std::borrow::Borrow;

use super::{lattice::*, nzq::*, poly::*, ring::*};

pub const ZZ_MAT: MatrixStructure<IntegerRing> = MatrixStructure { ring: &ZZ };
pub const QQ_MAT: MatrixStructure<RationalField> = MatrixStructure { ring: &QQ };

#[derive(Debug)]
pub enum MatOppErr {
    DimMissmatch,
    InvalidIndex,
    NotSquare,
}

#[derive(Debug, Clone)]
pub struct Matrix<ElemT: Clone> {
    dim1: usize,
    dim2: usize,
    transpose: bool,
    elems: Vec<ElemT>, //length self.rows * self.cols. row r and column c is index c + r * self.cols
}

impl<ElemT: Clone> Matrix<ElemT> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        if self.elems.len() != self.dim1 * self.dim2 {
            return Err("matrix entries has the wrong length");
        }
        Ok(())
    }

    pub fn full(rows: usize, cols: usize, elem: &ElemT) -> Self {
        let mut elems = Vec::with_capacity(rows * cols);
        for _i in 0..rows * cols {
            elems.push(elem.clone());
        }
        Self {
            dim1: rows,
            dim2: cols,
            transpose: false,
            elems,
        }
    }

    pub fn construct(rows: usize, cols: usize, make_entry: impl Fn(usize, usize) -> ElemT) -> Self {
        let mut elems = Vec::with_capacity(rows * cols);
        for idx in 0..rows * cols {
            let (r, c) = (idx / cols, idx % cols); //idx_to_rc for transpose=false
            elems.push(make_entry(r, c).clone());
        }
        Self {
            dim1: rows,
            dim2: cols,
            transpose: false,
            elems,
        }
    }

    pub fn from_rows(rows_elems: Vec<Vec<ElemT>>) -> Self {
        let rows = rows_elems.len();
        assert!(rows >= 1);
        let cols = rows_elems[0].len();
        for r in 1..rows {
            assert_eq!(rows_elems[r].len(), cols);
        }
        Self::construct(rows, cols, |r, c| rows_elems[r][c].clone())
    }

    pub fn from_cols(cols_elems: Vec<Vec<ElemT>>) -> Self {
        Self::from_rows(cols_elems).transpose()
    }

    fn rc_to_idx(&self, r: usize, c: usize) -> usize {
        match self.transpose {
            false => c + r * self.dim2,
            true => r + c * self.dim2,
        }
    }

    pub fn at(&self, r: usize, c: usize) -> Result<&ElemT, MatOppErr> {
        if r >= self.rows() {
            Err(MatOppErr::InvalidIndex)
        } else if c >= self.cols() {
            Err(MatOppErr::InvalidIndex)
        } else {
            let idx = self.rc_to_idx(r, c);
            Ok(&self.elems[idx])
        }
    }

    pub fn at_mut(&mut self, r: usize, c: usize) -> Result<&mut ElemT, MatOppErr> {
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
            elems,
        }
    }

    pub fn apply_map<NewElemT: Clone>(&self, f: impl Fn(&ElemT) -> NewElemT) -> Matrix<NewElemT> {
        Matrix {
            dim1: self.dim1,
            dim2: self.dim2,
            transpose: self.transpose,
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
    }

    pub fn join_rows<MatT: Borrow<Matrix<ElemT>>>(cols: usize, mats: Vec<MatT>) -> Matrix<ElemT> {
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

    pub fn join_cols<MatT: Borrow<Matrix<ElemT>>>(rows: usize, mats: Vec<MatT>) -> Matrix<ElemT> {
        let mut t_mats = vec![];
        for mat in mats {
            t_mats.push(mat.borrow().clone().transpose());
        }
        let joined = Self::join_rows(rows, t_mats.iter().collect());
        joined.transpose()
    }
}

impl<ElemT: Clone + PartialEq> PartialEq for Matrix<ElemT> {
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

impl<ElemT: ComRing + PartialEq + Eq> Eq for Matrix<ElemT> {}

pub struct MatrixStructure<'a, R: ComRing> {
    ring: &'a R,
}

impl<'a, R: ComRing> MatrixStructure<'a, R> {
    pub fn new(ring: &'a R) -> Self {
        Self { ring }
    }
}

impl<'a, R: ComRing> MatrixStructure<'a, R> {
    pub fn pprint(&self, mat: &Matrix<R::ElemT>) {
        let mut str_rows = vec![];
        for r in 0..mat.rows() {
            str_rows.push(vec![]);
            for c in 0..mat.cols() {
                str_rows[r].push(self.ring.to_string(mat.at(r, c).unwrap()));
            }
        }
        let cols_widths: Vec<usize> = (0..mat.cols())
            .map(|c| {
                (0..mat.rows())
                    .map(|r| str_rows[r][c].len())
                    .fold(0usize, |a, b| a.max(b))
            })
            .collect();

        for r in 0..mat.rows() {
            for c in 0..mat.cols() {
                while str_rows[r][c].len() < cols_widths[c] {
                    str_rows[r][c].push(' ');
                }
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

    pub fn zero(&self, rows: usize, cols: usize) -> Matrix<R::ElemT> {
        Matrix::construct(rows, cols, |_r, _c| self.ring.zero())
    }

    pub fn ident(&self, n: usize) -> Matrix<R::ElemT> {
        Matrix::construct(n, n, |r, c| {
            if r == c {
                self.ring.one()
            } else {
                self.ring.zero()
            }
        })
    }

    pub fn diag(&self, diag: &Vec<R::ElemT>) -> Matrix<R::ElemT> {
        Matrix::construct(diag.len(), diag.len(), |r, c| {
            if r == c {
                diag[r].clone()
            } else {
                self.ring.zero()
            }
        })
    }
}

impl<'a, R: ComRing> MatrixStructure<'a, R> {
    pub fn add_mut(&self, a: &mut Matrix<R::ElemT>, b: &Matrix<R::ElemT>) -> Result<(), MatOppErr> {
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
        mut a: Matrix<R::ElemT>,
        b: Matrix<R::ElemT>,
    ) -> Result<Matrix<R::ElemT>, MatOppErr> {
        match self.add_mut(&mut a, &b) {
            Ok(()) => Ok(a),
            Err(e) => Err(e),
        }
    }

    pub fn add_ref(
        &self,
        mut a: Matrix<R::ElemT>,
        b: &Matrix<R::ElemT>,
    ) -> Result<Matrix<R::ElemT>, MatOppErr> {
        match self.add_mut(&mut a, b) {
            Ok(()) => Ok(a),
            Err(e) => Err(e),
        }
    }

    pub fn add_refs(
        &self,
        a: &Matrix<R::ElemT>,
        b: &Matrix<R::ElemT>,
    ) -> Result<Matrix<R::ElemT>, MatOppErr> {
        let mut new_a = a.clone();
        match self.add_mut(&mut new_a, &b) {
            Ok(()) => Ok(new_a),
            Err(e) => Err(e),
        }
    }

    pub fn neg_mut(&self, a: &mut Matrix<R::ElemT>) {
        for r in 0..a.rows() {
            for c in 0..a.cols() {
                let neg_elem = self.ring.neg_ref(a.at(r, c).unwrap());
                *a.at_mut(r, c).unwrap() = neg_elem;
            }
        }
    }

    pub fn neg(&self, mut a: Matrix<R::ElemT>) -> Matrix<R::ElemT> {
        self.neg_mut(&mut a);
        a
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

    pub fn mul_refs(
        &self,
        a: &Matrix<R::ElemT>,
        b: &Matrix<R::ElemT>,
    ) -> Result<Matrix<R::ElemT>, MatOppErr> {
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
                        &self.ring.mul_refs(a.at(r, m).unwrap(), b.at(m, c).unwrap()),
                    );
                }
            }
        }
        Ok(s)
    }

    pub fn mul_scalar(&self, mut a: Matrix<R::ElemT>, scalar: &R::ElemT) -> Matrix<R::ElemT> {
        for r in 0..a.rows() {
            for c in 0..a.cols() {
                self.ring.mul_mut(a.at_mut(r, c).unwrap(), scalar);
            }
        }
        a
    }

    pub fn mul_scalar_ref(&self, a: &Matrix<R::ElemT>, scalar: &R::ElemT) -> Matrix<R::ElemT> {
        self.mul_scalar(a.clone(), scalar)
    }

    pub fn det_naive(&self, a: &Matrix<R::ElemT>) -> Result<R::ElemT, MatOppErr> {
        let n = a.dim1;
        if n != a.dim2 {
            Err(MatOppErr::NotSquare)
        } else {
            let mut det = self.ring.zero();
            for perm in super::super::sets::permutations::all_perms(n) {
                let mut prod = self.ring.one();
                for k in 0..n {
                    self.ring
                        .mul_mut(&mut prod, a.at(k, perm.call(k).unwrap()).unwrap());
                }
                if !perm.sign() {
                    self.ring.neg_mut(&mut prod);
                }
                self.ring.add_mut(&mut det, &prod);
            }
            Ok(det)
        }
    }
}

#[derive(Debug)]
enum ElementaryOppType<R: ComRing> {
    //swap distinct rows
    Swap(usize, usize),
    //multiply a row by a unit
    UnitMul {
        row: usize,
        unit: R::ElemT,
    },
    //row(i) -> row(i) + x*row(j)
    AddRowMul {
        i: usize,
        j: usize,
        x: R::ElemT,
    },
    //apply invertible row operations to two rows
    // /a b\
    // \c d/
    //such that ad-bc is a unit
    TwoInv {
        i: usize,
        j: usize,
        a: R::ElemT,
        b: R::ElemT,
        c: R::ElemT,
        d: R::ElemT,
    },
}

struct ElementaryOpp<'a, R: ComRing> {
    ring: &'a R,
    transpose: bool, //false = row opp, true = column opp
    opp: ElementaryOppType<R>,
}

impl<'a, R: PrincipalIdealDomain> ElementaryOpp<'a, R> {
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
                if !self.ring.is_unit(unit.clone()) {
                    return Err("can only multiply a row by a unit");
                }
            }
            ElementaryOppType::TwoInv { i, j, a, b, c, d } => {
                if i == j {
                    return Err("rows must be distinct");
                }
                let m = Matrix::<R::ElemT> {
                    dim1: 2,
                    dim2: 2,
                    transpose: false,
                    elems: vec![a.clone(), b.clone(), c.clone(), d.clone()],
                };
                if !self
                    .ring
                    .is_unit(MatrixStructure { ring: self.ring }.det_naive(&m).unwrap())
                {
                    return Err("can only apply an invertible row opperation to two rows");
                }
            }
        }
        Ok(())
    }

    fn new_row_opp(ring: &'a R, opp: ElementaryOppType<R>) -> Self {
        Self {
            ring,
            transpose: false,
            opp,
        }
    }

    fn new_col_opp(ring: &'a R, opp: ElementaryOppType<R>) -> Self {
        Self {
            ring,
            transpose: true,
            opp,
        }
    }

    fn det(&self) -> R::ElemT {
        match &self.opp {
            ElementaryOppType::Swap(_i, _j) => self.ring.neg(self.ring.one()),
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
            } => self.ring.add(
                self.ring.mul_refs(a, d),
                self.ring.neg(self.ring.mul_refs(b, c)),
            ),
        }
    }

    fn apply(&self, m: &mut Matrix<R::ElemT>) {
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
                    let offset = self.ring.mul_refs(m.at(*j, col).unwrap(), x);
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
                        self.ring.mul_refs(c, m.at(*i, col).unwrap()),
                        self.ring.mul_refs(d, m.at(*j, col).unwrap()),
                    );
                    // row(i) = a*row(i) + b*row(j)
                    *m.at_mut(*i, col).unwrap() = self.ring.add(
                        self.ring.mul_refs(a, m.at(*i, col).unwrap()),
                        self.ring.mul_refs(b, m.at(*j, col).unwrap()),
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

impl<'a, R: PrincipalIdealDomain> MatrixStructure<'a, R> {
    pub fn row_span(&self, a: Matrix<R::ElemT>) -> LinearLattice<R::ElemT> {
        LinearLatticeStructure::new(self.ring).from_span(
            1,
            a.cols(),
            (0..a.rows())
                .map(|r| a.submatrix(vec![r], (0..a.cols()).collect()))
                .collect(),
        )
    }

    pub fn col_span(&self, a: Matrix<R::ElemT>) -> LinearLattice<R::ElemT> {
        LinearLatticeStructure::new(self.ring).from_span(
            a.rows(),
            1,
            (0..a.cols())
                .map(|c| a.submatrix((0..a.rows()).collect(), vec![c]))
                .collect(),
        )
    }

    pub fn row_kernel(&self, a: Matrix<R::ElemT>) -> LinearLattice<R::ElemT> {
        let (_h, u, _u_det, pivs) = self.row_hermite_algorithm(a);
        LinearLatticeStructure::new(self.ring).from_basis(
            1,
            u.cols(),
            (pivs.len()..u.rows())
                .into_iter()
                .map(|r| u.submatrix(vec![r], (0..u.cols()).collect()))
                .collect(),
        )
    }

    pub fn col_kernel(&self, a: Matrix<R::ElemT>) -> LinearLattice<R::ElemT> {
        let (_h, u, _u_det, pivs) = self.col_hermite_algorithm(a);
        LinearLatticeStructure::new(self.ring).from_basis(
            u.rows(),
            1,
            (pivs.len()..u.cols())
                .into_iter()
                .map(|c| u.submatrix((0..u.rows()).collect(), vec![c]))
                .collect(),
        )
    }

    pub fn row_solve<VecT: Borrow<Matrix<R::ElemT>>>(
        &self,
        m: &Matrix<R::ElemT>,
        y: VecT,
    ) -> Option<Matrix<R::ElemT>> {
        match self.col_solve(&m.transpose_ref(), &y.borrow().transpose_ref()) {
            Some(x) => Some(x.transpose()),
            None => None,
        }
    }

    pub fn col_solve<VecT: Borrow<Matrix<R::ElemT>>>(
        &self,
        m: &Matrix<R::ElemT>,
        y: VecT,
    ) -> Option<Matrix<R::ElemT>> {
        assert_eq!(y.borrow().rows(), m.rows());
        assert_eq!(y.borrow().cols(), 1);
        //the kernel of ext_mat is related to the solution
        let ext_mat = Matrix::join_cols(m.rows(), vec![y.borrow(), m]);
        //we are looking for a point in the column kernel where the first coordinate is 1
        let col_ker = self.col_kernel(ext_mat);

        let first_coords: Vec<&R::ElemT> = (0..LinearLatticeStructure::new(self.ring)
            .rank(&col_ker))
            .map(|basis_num| {
                LinearLatticeStructure::new(self.ring)
                    .basis_matrix_element(&col_ker, basis_num, 0, 0)
            })
            .collect();

        let (g, taps) = self.ring.xgcd_list(first_coords);

        if self.ring.is_unit(g.clone()) {
            debug_assert_eq!(g, self.ring.one());
        }
        if g == self.ring.one() {
            //there is a solution
            //it is given by -(sum(taps * col_ker.basis)) with the first coordinate (equal to 1) removed
            let mut ext_ans = self.zero(m.cols() + 1, 1);
            for basis_num in 0..LinearLatticeStructure::new(self.ring).rank(&col_ker) {
                self.add_mut(
                    &mut ext_ans,
                    &self.mul_scalar_ref(
                        &LinearLatticeStructure::new(self.ring).basis_matrix(&col_ker, basis_num),
                        &taps[basis_num],
                    ),
                )
                .unwrap();
            }
            debug_assert_eq!(ext_ans.at(0, 0).unwrap(), &self.ring.one());
            let x = self.neg(ext_ans.submatrix((1..ext_ans.rows()).collect(), vec![0]));
            debug_assert_eq!(&self.mul_refs(m, &x).unwrap(), y.borrow());
            Some(x)
        } else {
            None //there is no solution
        }
    }

    pub fn row_solution_lattice<VecT: Borrow<Matrix<R::ElemT>>>(
        &self,
        m: &Matrix<R::ElemT>,
        y: VecT,
    ) -> AffineLattice<R::ElemT> {
        match self.row_solve(m, y) {
            Some(x) => AffineLatticeStructure::new(self.ring).from_offset_and_linear_lattice(
                1,
                m.rows(),
                x,
                self.row_kernel(m.clone()),
            ),
            None => AffineLatticeStructure::new(self.ring).empty(1, m.rows()),
        }
    }

    pub fn col_solution_lattice<VecT: Borrow<Matrix<R::ElemT>>>(
        &self,
        m: &Matrix<R::ElemT>,
        y: VecT,
    ) -> AffineLattice<R::ElemT> {
        match self.col_solve(m, y) {
            Some(x) => AffineLatticeStructure::new(self.ring).from_offset_and_linear_lattice(
                m.cols(),
                1,
                x,
                self.col_kernel(m.clone()),
            ),
            None => AffineLatticeStructure::new(self.ring).empty(m.cols(), 1),
        }
    }

    //TODO: replace with over a pid
    //if A:=self return (H, U, u_det, pivots) such that
    //H is in row hermite normal form
    //U is invertible
    //H=UA
    //u det is the determinant of u. It is a unit
    //pivots[r] is the column of the rth pivot and pivots.len() == rank(A)
    pub fn row_hermite_algorithm(
        &self,
        mut m: Matrix<R::ElemT>,
    ) -> (Matrix<R::ElemT>, Matrix<R::ElemT>, R::ElemT, Vec<usize>) {
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
                    if m.at(r, pc).unwrap() != &self.ring.zero() {
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
                    if a != &self.ring.zero() || b != &self.ring.zero() {
                        let (d, x, y) = self.ring.xgcd(a.clone(), b.clone());
                        debug_assert_eq!(
                            self.ring
                                .add(self.ring.mul_refs(&x, a), self.ring.mul_refs(&y, b)),
                            d
                        );
                        // perform the following row opps on self
                        // / x  -b/d \
                        // \ y   a/d /
                        let row_opp = ElementaryOpp::new_row_opp(
                            self.ring,
                            ElementaryOppType::TwoInv {
                                i: pr,
                                j: r,
                                a: x,
                                b: y,
                                //TODO: compute b/d and a/d at the same time d is computed?
                                c: self.ring.neg(self.ring.div(b.clone(), d.clone()).unwrap()),
                                d: self.ring.div(a.clone(), d.clone()).unwrap(),
                            },
                        );
                        //this will implicitly put the pivot into fav assoc form because that is what the gcd returns
                        row_opp.apply(&mut m);
                        row_opp.apply(&mut u);
                        self.ring.mul_mut(&mut u_det, &row_opp.det());
                    }
                }
            } else {
                //explicitly put the pivot into fav assoc form
                let (unit, _assoc) = self.ring.factor_fav_assoc_ref(m.at(pr, pc).unwrap());
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring,
                    ElementaryOppType::UnitMul {
                        row: pr,
                        unit: self.ring.inv(unit).unwrap(),
                    },
                );
                //this will implicitly put the pivot into fav assoc form because that is what the gcd returns
                row_opp.apply(&mut m);
                row_opp.apply(&mut u);
                self.ring.mul_mut(&mut u_det, &row_opp.det());
            }

            //should have eliminated everything below the pivot
            for r in pr + 1..m.rows() {
                debug_assert_eq!(m.at(r, pc).unwrap(), &self.ring.zero());
            }
            pr += 1;
        }

        if m.rows() <= 4 {
            debug_assert_eq!(self.det_naive(&u).unwrap(), u_det);
        }

        (m, u, u_det, pivs)
    }

    pub fn col_hermite_algorithm(
        &self,
        a: Matrix<R::ElemT>,
    ) -> (Matrix<R::ElemT>, Matrix<R::ElemT>, R::ElemT, Vec<usize>) {
        let (rh, ru, u_det, pivs) = self.row_hermite_algorithm(a.transpose());
        (rh.transpose(), ru.transpose(), u_det, pivs)
    }

    pub fn det(&self, a: Matrix<R::ElemT>) -> Result<R::ElemT, MatOppErr> {
        let n = a.rows();
        if n != a.cols() {
            Err(MatOppErr::NotSquare)
        } else {
            let (h, _u, u_det, _pivs) = self.row_hermite_algorithm(a);
            //h = u * self, we know det(u), and h is upper triangular
            let mut h_det = self.ring.one();
            for i in 0..n {
                self.ring.mul_mut(&mut h_det, h.at(i, i).unwrap());
            }
            Ok(self.ring.div(h_det, u_det).unwrap())
        }
    }

    pub fn rank(&self, a: Matrix<R::ElemT>) -> usize {
        let (_h, _u, _u_det, pivs) = self.row_hermite_algorithm(a);
        pivs.len()
    }

    //return (u, s, v, k) such that self = usv and s is in smith normal form and u, v are invertible and k is the number of non-zero elements in the diagonal of s
    pub fn smith_algorithm(
        &self,
        mut m: Matrix<R::ElemT>,
    ) -> (Matrix<R::ElemT>, Matrix<R::ElemT>, Matrix<R::ElemT>, usize) {
        let mut u = self.ident(m.rows());
        let mut v = self.ident(m.cols());

        let mut n = 0;
        'inductive_loop: while n < m.rows() && n < m.cols() {
            //search for a non-zero element to make the new starting point for (n, n)
            //having a non-zero element is necessary later in the algorithm
            'search_for_nonzero_element: {
                if m.at(n, n).unwrap() == &self.ring.zero() {
                    //search the first row to start with
                    for c in n + 1..m.cols() {
                        if m.at(n, c).unwrap() != &self.ring.zero() {
                            //swap column n and column c
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring,
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
                            if m.at(r, c).unwrap() != &self.ring.zero() {
                                //swap column n and column c
                                let col_opp = ElementaryOpp::new_col_opp(
                                    self.ring,
                                    ElementaryOppType::Swap(n, c),
                                );
                                col_opp.apply(&mut m);
                                col_opp.apply(&mut v);

                                //swap row n and row r
                                let row_opp = ElementaryOpp::new_col_opp(
                                    self.ring,
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
            let (unit, _assoc) = self.ring.factor_fav_assoc_ref(m.at(n, n).unwrap());
            let row_opp = ElementaryOpp::new_row_opp(
                self.ring,
                ElementaryOppType::UnitMul {
                    row: n,
                    unit: self.ring.inv(unit).unwrap(),
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
                    match self.ring.div_refs(b, a) {
                        Ok(q) => {
                            //b is a multiple of a
                            //replace (a, b) with (a, 0) by subtracting a multiple of a from b
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring,
                                ElementaryOppType::AddRowMul {
                                    i: c,
                                    j: n,
                                    x: self.ring.neg(q),
                                },
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut v);
                        }
                        Err(RingDivisionError::NotDivisible) => {
                            all_divisible = false;
                            //b is not a multiple of a
                            //replace (a, b) with (gcd, 0)
                            let (d, x, y) = self.ring.xgcd(a.clone(), b.clone());
                            debug_assert_eq!(
                                self.ring
                                    .add(self.ring.mul_refs(&x, a), self.ring.mul_refs(&y, b)),
                                d
                            );
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring,
                                ElementaryOppType::TwoInv {
                                    i: n,
                                    j: c,
                                    a: x,
                                    b: y,
                                    c: self.ring.neg(self.ring.div(b.clone(), d.clone()).unwrap()),
                                    d: self.ring.div(a.clone(), d.clone()).unwrap(),
                                },
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut v);
                        }
                        Err(RingDivisionError::DivideByZero) => {
                            //swap a and b
                            //a=0 so this does have the effect of (a, b) -> (gcd(a, b), 0)
                            let col_opp = ElementaryOpp::new_col_opp(
                                self.ring,
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
                    match self.ring.div_refs(b, a) {
                        Ok(q) => {
                            //b is a multiple of a
                            //replace (a, b) with (a, 0) by subtracting a multiple of a from b
                            let col_opp = ElementaryOpp::new_row_opp(
                                self.ring,
                                ElementaryOppType::AddRowMul {
                                    i: r,
                                    j: n,
                                    x: self.ring.neg(q),
                                },
                            );
                            col_opp.apply(&mut m);
                            col_opp.apply(&mut u);
                        }
                        Err(RingDivisionError::NotDivisible) => {
                            all_divisible = false;
                            //b is not a multiple of a
                            //replace (a, b) with (gcd, 0)
                            let (d, x, y) = self.ring.xgcd(a.clone(), b.clone());
                            debug_assert_eq!(
                                self.ring
                                    .add(self.ring.mul_refs(&x, a), self.ring.mul_refs(&y, b)),
                                d
                            );
                            let row_opp = ElementaryOpp::new_row_opp(
                                self.ring,
                                ElementaryOppType::TwoInv {
                                    i: n,
                                    j: r,
                                    a: x,
                                    b: y,
                                    c: self.ring.neg(self.ring.div(b.clone(), d.clone()).unwrap()),
                                    d: self.ring.div(a.clone(), d.clone()).unwrap(),
                                },
                            );
                            row_opp.apply(&mut m);
                            row_opp.apply(&mut u);
                        }
                        Err(RingDivisionError::DivideByZero) => {
                            //swap a and b
                            //a=0 so this does have the effect of (a, b) -> (gcd(a, b), 0)
                            let col_opp = ElementaryOpp::new_row_opp(
                                self.ring,
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
            debug_assert_ne!(m.at(n, n).unwrap(), &self.ring.zero());
            //some more fiddling is needed now to make the top left element divides everything else
            for r in n + 1..m.rows() {
                //row(n) = row(n) + row(r)
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring,
                    ElementaryOppType::AddRowMul {
                        i: n,
                        j: r,
                        x: self.ring.one(),
                    },
                );
                row_opp.apply(&mut m);
                row_opp.apply(&mut u);

                //row(n) goes from (g, a1, a2, ..., an) to (gcd, 0, 0, ..., 0)
                for c in n + 1..m.cols() {
                    let a = m.at(n, n).unwrap();
                    let b = m.at(n, c).unwrap();
                    //if a=0 and b=0 then there is nothing to do and the following step would fail when dividing by g=0
                    //b might not be a multiple of a
                    //replace (a, b) with (gcd, 0) to fix this
                    let (g, x, y) = self.ring.xgcd(a.clone(), b.clone());
                    debug_assert_eq!(
                        self.ring
                            .add(self.ring.mul_refs(&x, a), self.ring.mul_refs(&y, b)),
                        g
                    );
                    let col_opp = ElementaryOpp::new_col_opp(
                        self.ring,
                        ElementaryOppType::TwoInv {
                            i: n,
                            j: c,
                            a: x,
                            b: y,
                            c: self.ring.neg(self.ring.div(b.clone(), g.clone()).unwrap()),
                            d: self.ring.div(a.clone(), g.clone()).unwrap(),
                        },
                    );
                    col_opp.apply(&mut m);
                    col_opp.apply(&mut v);
                }

                //fix the first column
                for fix_r in n + 1..m.rows() {
                    let a = m.at(n, n).unwrap();
                    let b = m.at(fix_r, n).unwrap();
                    let q = self.ring.div_refs(b, a).unwrap();
                    let col_opp = ElementaryOpp::new_row_opp(
                        self.ring,
                        ElementaryOppType::AddRowMul {
                            i: fix_r,
                            j: n,
                            x: self.ring.neg(q),
                        },
                    );
                    col_opp.apply(&mut m);
                    col_opp.apply(&mut u);
                }
            }

            if m.at(n, n).unwrap() == &self.ring.zero() {
                //the bottom right submatrix is all zero
                break 'inductive_loop;
            }
            n += 1;
        }

        (u, m, v, n)
    }
}

impl<'a, R: EuclideanDomain + FavoriteAssociate> MatrixStructure<'a, R> {
    //if A:=self return (H, U, pivots) such that
    //H is in row reduced hermite normal form
    //U is invertible
    //H=UA
    //pivots[r] is the column of the rth pivot and pivots.len() == rank(A)
    pub fn row_reduced_hermite_algorithm(
        &self,
        m: Matrix<R::ElemT>,
    ) -> (Matrix<R::ElemT>, Matrix<R::ElemT>, Vec<usize>) {
        let (mut h, mut u, _u_det, pivs) = self.row_hermite_algorithm(m);

        for (pr, pc) in pivs.iter().enumerate() {
            for r in 0..pr {
                //reduce h[r, pc] so that it has norm less than h[pr, pc]
                let a = h.at(r, *pc).unwrap();
                let b = h.at(pr, *pc).unwrap();
                //a = b*q + r
                let q = self.ring.quo_refs(a, b).unwrap();
                let row_opp = ElementaryOpp::new_row_opp(
                    self.ring,
                    ElementaryOppType::AddRowMul {
                        i: r,
                        j: pr,
                        x: self.ring.neg(q),
                    },
                );
                row_opp.apply(&mut h);
                row_opp.apply(&mut u);
            }
        }

        (h, u, pivs)
    }

    pub fn col_reduced_hermite_algorithm(
        &self,
        m: Matrix<R::ElemT>,
    ) -> (Matrix<R::ElemT>, Matrix<R::ElemT>, Vec<usize>) {
        let (rh, ru, pivs) = self.row_reduced_hermite_algorithm(m.transpose());
        (rh.transpose(), ru.transpose(), pivs)
    }
}

impl<'a, F: Field> MatrixStructure<'a, F> {
    pub fn presentation_matrix(
        &self,
        m: Matrix<F::ElemT>,
    ) -> Result<Matrix<Polynomial<F::ElemT>>, MatOppErr> {
        let n = m.rows();
        if n != m.cols() {
            Err(MatOppErr::NotSquare)
        } else {
            let poly_ring = PolynomialRing::new(self.ring);
            let poly_mat_struct = MatrixStructure { ring: &poly_ring };
            Ok(poly_mat_struct
                .add(
                    m.apply_map(|x| poly_ring.constant(x.clone())),
                    poly_mat_struct
                        .neg(poly_mat_struct.diag(&(0..n).map(|_i| poly_ring.var()).collect())),
                )
                .unwrap())
        }
    }

    pub fn minimal_polynomial(
        &self,
        m: Matrix<F::ElemT>,
    ) -> Result<Polynomial<F::ElemT>, MatOppErr> {
        match self.presentation_matrix(m) {
            Ok(pres_mat) => {
                let poly_ring = PolynomialRing::new(self.ring);
                let poly_mat_struct = MatrixStructure::new(&poly_ring);
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
        m: Matrix<F::ElemT>,
    ) -> Result<Polynomial<F::ElemT>, MatOppErr> {
        match self.presentation_matrix(m) {
            Ok(pres_mat) => {
                let poly_ring = PolynomialRing::new(self.ring);
                let poly_mat_struct = MatrixStructure::new(&poly_ring);
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

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;
    use malachite_q::Rational;

    use super::super::nzq::*;
    use super::*;

    #[test]
    fn test_join_rows() {
        let top = Matrix::from_rows(vec![vec![
            Integer::from(1),
            Integer::from(2),
            Integer::from(3),
        ],
        vec![
            Integer::from(4),
            Integer::from(5),
            Integer::from(6),
        ]]);
        let bot = Matrix::from_rows(vec![vec![
            Integer::from(7),
            Integer::from(8),
            Integer::from(9),
        ]]);

        let both = Matrix::from_rows(vec![
            vec![Integer::from(1), Integer::from(2), Integer::from(3)],
            vec![Integer::from(4), Integer::from(5), Integer::from(6)],
            vec![Integer::from(7), Integer::from(8), Integer::from(9)],
        ]);

        println!("top");
        ZZ_MAT.pprint(&top);
        println!("bot");
        ZZ_MAT.pprint(&bot);
        println!("both");
        ZZ_MAT.pprint(&both);

        let ans = Matrix::join_rows(3, vec![top, bot]);
        println!("ans");
        ZZ_MAT.pprint(&ans);

        assert_eq!(ans, both);
    }

    #[test]
    fn invariants() {
        let m = Matrix {
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

        let m = Matrix {
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
        let a = Matrix {
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

        let b = Matrix {
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
            let mut a = Matrix {
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

            let b = Matrix {
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

            let c = Matrix {
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

            ZZ_MAT.add_mut(&mut a, &b).unwrap();

            assert_eq!(a, c);
        }

        {
            let mut a = Matrix {
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

            let b = Matrix {
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

            let c = Matrix {
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

            ZZ_MAT.add_mut(&mut a, &b).unwrap();

            assert_eq!(a, c);
        }

        {
            let mut a = Matrix {
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

            let b = Matrix {
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

            match ZZ_MAT.add_mut(&mut a, &b) {
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

            assert_eq!(ZZ_MAT.add(a, b).unwrap(), c);
        }
    }

    #[test]
    fn mul() {
        {
            let a = Matrix {
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

            let b = Matrix {
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

            let c = Matrix {
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

            assert_eq!(ZZ_MAT.mul_refs(&a, &b).unwrap(), c);
        }
    }

    #[test]
    fn det_naive() {
        let m = Matrix::from_rows(vec![
            vec![Integer::from(1), Integer::from(3), Integer::from(2)],
            vec![Integer::from(-3), Integer::from(-1), Integer::from(-3)],
            vec![Integer::from(2), Integer::from(3), Integer::from(1)],
        ]);
        assert_eq!(ZZ_MAT.det_naive(&m).unwrap(), Integer::from(-15));
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
            ZZ_MAT.pprint(&a);
            let (h, u, pivs) = ZZ_MAT.row_reduced_hermite_algorithm(a.clone());
            println!("H =");
            ZZ_MAT.pprint(&h);
            println!("pivs = {:?}", pivs);
            println!("U =");
            ZZ_MAT.pprint(&u);
            assert_eq!(h, ZZ_MAT.mul_refs(&u, &a).unwrap());

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
                    assert_eq!(h.at(r, cz).unwrap(), &ZZ.zero());
                }
            }

            //check pivot columns
            for (pr, pc) in pivs.iter().enumerate() {
                assert!(h.at(pr, *pc).unwrap() != &ZZ.zero());
                for r in 0..h.rows() {
                    if r > pr {
                        assert_eq!(h.at(r, *pc).unwrap(), &ZZ.zero());
                    } else if r == pr {
                        let (_unit, assoc) = ZZ.factor_fav_assoc(h.at(r, *pc).unwrap().clone());
                        assert_eq!(&assoc, h.at(r, *pc).unwrap());
                    } else {
                        assert!(ZZ.norm(h.at(r, *pc).unwrap()) < ZZ.norm(h.at(pr, *pc).unwrap()));
                    }
                }
            }

            println!();
            println!("hermite reduced col algorithm for");
            ZZ_MAT.pprint(&a);
            let (h, u, pivs) = ZZ_MAT.col_reduced_hermite_algorithm(a.clone());
            println!("H =");
            ZZ_MAT.pprint(&h);
            println!("pivs = {:?}", pivs);
            println!("U =");
            ZZ_MAT.pprint(&u);

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
                    assert_eq!(h.at(rz, c).unwrap(), &ZZ.zero());
                }
            }

            //check the pivot rows
            assert_eq!(h, ZZ_MAT.mul_refs(&a, &u).unwrap());
            for (pc, pr) in pivs.iter().enumerate() {
                assert!(h.at(*pr, pc).unwrap() != &ZZ.zero());
                for c in 0..h.cols() {
                    if c > pc {
                        assert_eq!(h.at(*pr, c).unwrap(), &ZZ.zero());
                    } else if c == pc {
                        let (_unit, assoc) = ZZ.factor_fav_assoc(h.at(*pr, c).unwrap().clone());
                        assert_eq!(&assoc, h.at(*pr, c).unwrap());
                    } else {
                        assert!(ZZ.norm(h.at(*pr, c).unwrap()) < ZZ.norm(h.at(*pr, pc).unwrap()));
                    }
                }
            }
        }

        {
            //integer reduced hermite normal form is unique, so we can fully check an example
            let a = Matrix::from_rows(vec![
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

            let (h, u, pivs) = ZZ_MAT.row_reduced_hermite_algorithm(a.clone());

            assert_eq!(h, expected_h);
            assert_eq!(u, expected_u);
            assert_eq!(pivs, vec![0, 1, 2]);
        }

        {
            //this one used to cause a dividion by zero error when replacing (a, b) with (gcd, 0) when (a, b) = (0, 0)
            let a = Matrix::from_rows(vec![
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

            let (h, u, pivs) = ZZ_MAT.row_reduced_hermite_algorithm(a.clone());

            assert_eq!(h, expected_h);
            assert_eq!(u, expected_u);
            assert_eq!(pivs, vec![0, 1, 2]);
        }
    }

    #[test]
    fn smith_algorithm() {
        {
            let a = Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(4), Integer::from(4)],
                vec![Integer::from(-6), Integer::from(6), Integer::from(12)],
                vec![Integer::from(10), Integer::from(4), Integer::from(16)],
            ]);
            let (u, s, v, k) = ZZ_MAT.smith_algorithm(a.clone());
            assert_eq!(
                s,
                ZZ_MAT
                    .mul_refs(&ZZ_MAT.mul_refs(&u, &a).unwrap(), &v)
                    .unwrap()
            );
            assert_eq!(k, 3);
            assert_eq!(
                s,
                Matrix::from_rows(vec![
                    vec![Integer::from(2), Integer::from(0), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(2), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(0), Integer::from(156)]
                ])
            );

            let a = Matrix::from_rows(vec![
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
            let (u, s, v, k) = ZZ_MAT.smith_algorithm(a.clone());
            assert_eq!(
                s,
                ZZ_MAT
                    .mul_refs(&ZZ_MAT.mul_refs(&u, &a).unwrap(), &v)
                    .unwrap()
            );
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
            let a = Matrix::from_rows(vec![
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
            let (_u, _s, _v, k) = ZZ_MAT.smith_algorithm(a.clone());
            assert_eq!(k, 1);
        }
    }
    #[test]
    fn min_and_char_polys() {
        {
            let a = Matrix::from_rows(vec![
                vec![Integer::from(0), Integer::from(4), Integer::from(4)],
                vec![Integer::from(1), Integer::from(4), Integer::from(16)],
            ]);
            let (_u, s, _v, _k) = ZZ_MAT.smith_algorithm(a);

            assert_eq!(
                s,
                Matrix::from_rows(vec![
                    vec![Integer::from(1), Integer::from(0), Integer::from(0),],
                    vec![Integer::from(0), Integer::from(4), Integer::from(0)],
                ])
            );
        }

        {
            let a = Matrix::from_rows(vec![
                vec![Rational::from(0), Rational::from(0), Rational::from(0)],
                vec![Rational::from(0), Rational::from(0), Rational::from(1)],
                vec![Rational::from(0), Rational::from(0), Rational::from(0)],
            ]);
            let min_p = QQ_MAT.minimal_polynomial(a.clone()).unwrap();
            let char_p = QQ_MAT.characteristic_polynomial(a.clone()).unwrap();
            assert_eq!(
                min_p,
                QQ_POLY.from_coeffs(vec![
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(1)
                ])
            );
            assert_eq!(
                char_p,
                QQ_POLY.from_coeffs(vec![
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(1)
                ])
            );
        }
    }

    #[test]
    fn span_and_kernel_rank() {
        let mat = Matrix::from_rows(vec![
            vec![Integer::from(1), Integer::from(1), Integer::from(1)],
            vec![Integer::from(1), Integer::from(2), Integer::from(1)],
            vec![Integer::from(1), Integer::from(1), Integer::from(1)],
            vec![Integer::from(1), Integer::from(1), Integer::from(1)],
        ]);

        assert_eq!(ZZ_LINLAT.rank(&ZZ_MAT.row_span(mat.clone())), 2);
        assert_eq!(ZZ_LINLAT.rank(&ZZ_MAT.col_span(mat.clone())), 2);
        assert_eq!(ZZ_LINLAT.rank(&ZZ_MAT.row_kernel(mat.clone())), 2);
        assert_eq!(ZZ_LINLAT.rank(&ZZ_MAT.col_kernel(mat.clone())), 1);
    }

    #[test]
    fn span_and_kernel_points() {
        let mat = Matrix::from_rows(vec![
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
        ZZ_MAT.pprint(&mat);

        let k = ZZ_MAT.col_kernel(mat);
        println!("kernel");
        ZZ_LINLAT.pprint(&k);

        assert!(ZZ_LINLAT.contains_point(
            &k,
            Matrix::from_rows(vec![
                vec![Integer::from(-1)],
                vec![Integer::from(1)],
                vec![Integer::from(5)],
                vec![Integer::from(-3)]
            ])
        ));
    }
}
