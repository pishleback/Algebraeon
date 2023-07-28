#![allow(dead_code)]

use std::borrow::Borrow;

use super::{lattice::*, poly::*, ring::*};

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

pub fn join_rows<R: ComRing, MatT: Borrow<Matrix<R>>>(cols: usize, mats: Vec<MatT>) -> Matrix<R> {
    let mut rows = 0;
    for mat in &mats {
        assert_eq!(cols, mat.borrow().cols());
        rows += mat.borrow().rows();
    }
    let mut joined = Matrix::zero(rows, cols);
    let mut row_offset = 0;
    for mat in &mats {
        for r in 0..mat.borrow().rows() {
            for c in 0..cols {
                *joined.at_mut(row_offset + r, c).unwrap() = mat.borrow().at(r, c).unwrap().clone();
            }
        }
        row_offset += mat.borrow().rows();
    }
    joined
}

pub fn join_cols<R: ComRing, MatT: Borrow<Matrix<R>>>(rows: usize, mats: Vec<MatT>) -> Matrix<R> {
    let mut t_mats = vec![];
    for mat in mats {
        t_mats.push(mat.borrow().clone().transpose());
    }
    let joined = join_rows(rows, t_mats.iter().collect());
    joined.transpose()
}

impl<R: ComRing> Matrix<R> {
    fn check_invariants(&self) -> Result<(), &'static str> {
        if self.elems.len() != self.dim1 * self.dim2 {
            return Err("matrix entries has the wrong length");
        }
        Ok(())
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

    pub fn from_cols(cols_elems: Vec<Vec<R>>) -> Self {
        Self::from_rows(cols_elems).transpose()
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

    pub fn diag(diag: &Vec<R>) -> Self {
        let n = diag.len();
        let mut elems = Vec::with_capacity(n * n);
        for r in 0..n {
            for c in 0..n {
                match r == c {
                    true => elems.push(diag[r].clone()),
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

    pub fn apply_map<S: ComRing>(&self, f: impl Fn(&R) -> S) -> Matrix<S> {
        Matrix {
            dim1: self.dim1,
            dim2: self.dim2,
            transpose: self.transpose,
            elems: self.elems.iter().map(|x| f(x)).collect(),
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
}

impl<R: ComRing> Matrix<R> {
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
            if self.rows() == 1 {
                print!("( ");
            } else if r == 0 {
                print!("/ ");
            } else if r == self.rows() - 1 {
                print!("\\ ");
            } else {
                print!("| ");
            }
            for c in 0..self.cols() {
                if c != 0 {
                    print!("    ");
                }
                print!("{}", str_rows[r][c]);
            }
            if self.rows() == 1 {
                print!(" )");
            } else if r == 0 {
                print!(" \\");
            } else if r == self.rows() - 1 {
                print!(" /");
            } else {
                print!(" |");
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
            for c in 0..cols {
                for r in 0..rows {
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

    pub fn neg_mut(&mut self) {
        for r in 0..self.rows() {
            for c in 0..self.cols() {
                let neg_elem = self.at(r, c).unwrap().neg_ref();
                *self.at_mut(r, c).unwrap() = neg_elem;
            }
        }
    }

    pub fn neg(mut self) -> Self {
        self.neg_mut();
        self
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

    pub fn mul_scalar(mut self, scalar: &R) -> Matrix<R> {
        for r in 0..self.rows() {
            for c in 0..self.cols() {
                self.at_mut(r, c).unwrap().mul_mut(scalar);
            }
        }
        self
    }

    pub fn mul_scalar_ref(&self, scalar: &R) -> Matrix<R> {
        self.clone().mul_scalar(scalar)
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

#[derive(Debug)]
enum ElementaryOppType<R: ComRing> {
    //swap distinct rows
    Swap(usize, usize),
    //multiply a row by a unit
    UnitMul {
        row: usize,
        unit: R,
    },
    //row(i) -> row(i) + x*row(j)
    AddRowMul {
        i: usize,
        j: usize,
        x: R,
    },
    //apply invertible row operations to two rows
    // /a b\
    // \c d/
    //such that ad-bc is a unit
    TwoInv {
        i: usize,
        j: usize,
        a: R,
        b: R,
        c: R,
        d: R,
    },
}

struct ElementaryOpp<R: ComRing> {
    transpose: bool, //false = row opp, true = column opp
    opp: ElementaryOppType<R>,
}

impl<R: PrincipalIdealDomain> ElementaryOpp<R> {
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
                if !unit.clone().is_unit() {
                    return Err("can only multiply a row by a unit");
                }
            }
            ElementaryOppType::TwoInv { i, j, a, b, c, d } => {
                if i == j {
                    return Err("rows must be distinct");
                }
                let m = Matrix {
                    dim1: 2,
                    dim2: 2,
                    transpose: false,
                    elems: vec![a.clone(), b.clone(), c.clone(), d.clone()],
                };
                if !m.det_naive().unwrap().is_unit() {
                    return Err("can only apply an invertible row opperation to two rows");
                }
            }
        }
        Ok(())
    }

    fn new_row_opp(opp: ElementaryOppType<R>) -> Self {
        Self {
            transpose: false,
            opp,
        }
    }

    fn new_col_opp(opp: ElementaryOppType<R>) -> Self {
        Self {
            transpose: true,
            opp,
        }
    }

    fn det(&self) -> R {
        match &self.opp {
            ElementaryOppType::Swap(_i, _j) => R::one().neg(),
            ElementaryOppType::UnitMul { row: _row, unit } => unit.clone(),
            ElementaryOppType::AddRowMul {
                i: _i,
                j: _j,
                x: _x,
            } => R::one(),
            ElementaryOppType::TwoInv {
                i: _i,
                j: _j,
                a,
                b,
                c,
                d,
            } => R::add(R::mul_refs(a, d), R::mul_refs(b, c).neg()),
        }
    }

    fn apply(&self, m: &mut Matrix<R>) {
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
                    let offset = R::mul_refs(m.at(*j, col).unwrap(), x);
                    m.at_mut(*i, col).unwrap().add_mut(&offset)
                }
            }
            // /u 0\
            // \0 1/
            ElementaryOppType::UnitMul { row, unit } => {
                for col in 0..m.cols() {
                    m.at_mut(*row, col).unwrap().mul_mut(unit)
                }
            }
            // /a b\
            // \c d/
            ElementaryOppType::TwoInv { i, j, a, b, c, d } => {
                for col in 0..m.cols() {
                    // tmp = c*row(i) + d*row(j)
                    let tmp = R::add(
                        R::mul_refs(c, m.at(*i, col).unwrap()),
                        R::mul_refs(d, m.at(*j, col).unwrap()),
                    );
                    // row(i) = a*row(i) + b*row(j)
                    *m.at_mut(*i, col).unwrap() = R::add(
                        R::mul_refs(a, m.at(*i, col).unwrap()),
                        R::mul_refs(b, m.at(*j, col).unwrap()),
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

impl<R: PrincipalIdealDomain> Matrix<R> {
    pub fn row_span(&self) -> LinearLattice<R> {
        LinearLattice::from_span(
            1,
            self.cols(),
            (0..self.rows())
                .map(|r| self.submatrix(vec![r], (0..self.cols()).collect()))
                .collect(),
        )
    }

    pub fn col_span(&self) -> LinearLattice<R> {
        LinearLattice::from_span(
            self.rows(),
            1,
            (0..self.cols())
                .map(|c| self.submatrix((0..self.rows()).collect(), vec![c]))
                .collect(),
        )
    }

    pub fn row_kernel(self) -> LinearLattice<R> {
        let (_h, u, _u_det, pivs) = self.row_hermite_algorithm();
        LinearLattice::from_basis(
            1,
            u.cols(),
            (pivs.len()..u.rows())
                .into_iter()
                .map(|r| u.submatrix(vec![r], (0..u.cols()).collect()))
                .collect(),
        )
    }

    pub fn col_kernel(self) -> LinearLattice<R> {
        let (_h, u, pivs) = self.col_hermite_algorithm();
        LinearLattice::from_basis(
            u.rows(),
            1,
            (pivs.len()..u.cols())
                .into_iter()
                .map(|c| u.submatrix((0..u.rows()).collect(), vec![c]))
                .collect(),
        )
    }

    pub fn row_solve<VecT: Borrow<Matrix<R>>>(&self, y: VecT) -> Option<Matrix<R>> {
        match self.transpose_ref().col_solve(&y.borrow().transpose_ref()) {
            Some(x) => Some(x.transpose()),
            None => None,
        }
    }

    pub fn col_solve<VecT: Borrow<Matrix<R>>>(&self, y: VecT) -> Option<Matrix<R>> {
        assert_eq!(y.borrow().rows(), self.rows());
        assert_eq!(y.borrow().cols(), 1);
        //the kernel of ext_mat is related to the solution
        let ext_mat = join_cols(self.rows(), vec![y.borrow(), self]);
        //we are looking for a point in the column kernel where the first coordinate is 1
        let col_ker = ext_mat.col_kernel();

        let first_coords: Vec<&R> = (0..col_ker.rank())
            .map(|basis_num| col_ker.basis_matrix_element(basis_num, 0, 0))
            .collect();

        let (g, taps) = R::xgcd_list(first_coords);

        if g.clone().is_unit() {
            debug_assert_eq!(g, R::one());
        }
        if g == R::one() {
            //there is a solution
            //it is given by -(sum(taps * col_ker.basis)) with the first coordinate (equal to 1) removed
            let mut ext_ans = Matrix::zero(self.cols() + 1, 1);
            for basis_num in 0..col_ker.rank() {
                ext_ans
                    .add_mut(
                        &col_ker
                            .basis_matrix(basis_num)
                            .mul_scalar_ref(&taps[basis_num]),
                    )
                    .unwrap();
            }
            debug_assert_eq!(ext_ans.at(0, 0).unwrap(), &R::one());
            let x = ext_ans
                .submatrix((1..ext_ans.rows()).collect(), vec![0])
                .neg();
            debug_assert_eq!(&Self::mul_refs(self, &x).unwrap(), y.borrow());
            Some(x)
        } else {
            None //there is no solution
        }
    }

    pub fn row_solution_lattice<VecT: Borrow<Matrix<R>>>(&self, y: VecT) -> AffineLattice<R> {
        match self.row_solve(y) {
            Some(x) => AffineLattice::from_offset_and_linear_lattice(
                1,
                self.rows(),
                x,
                self.clone().row_kernel(),
            ),
            None => AffineLattice::empty(1, self.rows()),
        }
    }

    pub fn col_solution_lattice<VecT: Borrow<Matrix<R>>>(&self, y: VecT) -> AffineLattice<R> {
        match self.col_solve(y) {
            Some(x) => AffineLattice::from_offset_and_linear_lattice(
                self.cols(),
                1,
                x,
                self.clone().col_kernel(),
            ),
            None => AffineLattice::empty(self.cols(), 1),
        }
    }

    //TODO: replace with over a pid
    //if A:=self return (H, U, u_det, pivots) such that
    //H is in row hermite normal form
    //U is invertible
    //H=UA
    //u det is the determinant of u. It is a unit
    //pivots[r] is the column of the rth pivot and pivots.len() == rank(A)
    pub fn row_hermite_algorithm(mut self) -> (Self, Self, R, Vec<usize>) {
        //build up U by applying row opps to the identity as we go
        let mut u = Self::ident(self.rows());
        let mut u_det = R::one();
        let mut pivs = vec![];

        let (mut pr, mut pc) = (0, 0);
        'pivot_loop: while pr < self.rows() {
            //find the next pivot row
            //the next pivot row is the next row with a non-zero element below the previous pivot
            'next_pivot_loop: loop {
                if pc == self.cols() {
                    break 'pivot_loop;
                }

                for r in pr..self.rows() {
                    if self.at(r, pc).unwrap() != &R::zero() {
                        break 'next_pivot_loop;
                    }
                }

                pc += 1;
            }
            pivs.push(pc);

            if pr + 1 < self.rows() {
                //reduce everything below the pivot to zero
                for r in pr + 1..self.rows() {
                    let a = self.at(pr, pc).unwrap();
                    let b = self.at(r, pc).unwrap();
                    //if a=0 and b=0 there is nothing to do. The reduction step would fail because d=0 and we divide by d, so just skip it in this case
                    if a != &R::zero() || b != &R::zero() {
                        let (d, x, y) = R::xgcd(a.clone(), b.clone());
                        debug_assert_eq!(R::add(R::mul_refs(&x, a), R::mul_refs(&y, b)), d);
                        // perform the following row opps on self
                        // / x  -b/d \
                        // \ y   a/d /
                        let row_opp = ElementaryOpp::new_row_opp(ElementaryOppType::TwoInv {
                            i: pr,
                            j: r,
                            a: x,
                            b: y,
                            //TODO: compute b/d and a/d at the same time d is computed?
                            c: R::div(b.clone(), d.clone()).unwrap().neg(),
                            d: R::div(a.clone(), d.clone()).unwrap(),
                        });
                        //this will implicitly put the pivot into fav assoc form because that is what the gcd returns
                        row_opp.apply(&mut self);
                        row_opp.apply(&mut u);
                        u_det.mul_mut(&row_opp.det());
                    }
                }
            } else {
                //explicitly put the pivot into fav assoc form
                let (unit, _assoc) = self.at(pr, pc).unwrap().factor_fav_assoc_ref();
                let row_opp = ElementaryOpp::new_row_opp(ElementaryOppType::UnitMul {
                    row: pr,
                    unit: unit.inv().unwrap(),
                });
                //this will implicitly put the pivot into fav assoc form because that is what the gcd returns
                row_opp.apply(&mut self);
                row_opp.apply(&mut u);
                u_det.mul_mut(&row_opp.det());
            }

            //should have eliminated everything below the pivot
            for r in pr + 1..self.rows() {
                debug_assert_eq!(self.at(r, pc).unwrap(), &R::zero());
            }
            pr += 1;
        }

        if self.rows() <= 4 {
            debug_assert_eq!(u.det_naive().unwrap(), u_det);
        }

        (self, u, u_det, pivs)
    }

    pub fn col_hermite_algorithm(self) -> (Self, Self, Vec<usize>) {
        let (rh, ru, _u_det, pivs) = self.transpose().row_hermite_algorithm();
        (rh.transpose(), ru.transpose(), pivs)
    }

    pub fn det(self) -> Result<R, MatOppErr> {
        let n = self.rows();
        if n != self.cols() {
            Err(MatOppErr::NotSquare)
        } else {
            let (h, u, u_det, pivs) = self.row_hermite_algorithm();
            //h = u * self, we know det(u), and h is upper triangular
            let mut h_det = R::one();
            for i in 0..n {
                h_det.mul_mut(h.at(i, i).unwrap());
            }
            Ok(R::div(h_det, u_det).unwrap())
        }
    }

    pub fn rank(self) -> usize {
        let (_h, _u, _u_det, pivs) = self.row_hermite_algorithm();
        pivs.len()
    }

    //return (u, s, v, k) such that self = usv and s is in smith normal form and u, v are invertible and k is the number of non-zero elements in the diagonal of s
    pub fn smith_algorithm(mut self) -> (Self, Self, Self, usize) {
        let mut u = Self::ident(self.rows());
        let mut v = Self::ident(self.cols());

        let mut n = 0;
        'inductive_loop: while n < self.rows() && n < self.cols() {
            //search for a non-zero element to make the new starting point for (n, n)
            //having a non-zero element is necessary later in the algorithm
            'search_for_nonzero_element: {
                if self.at(n, n).unwrap() == &R::zero() {
                    //search the first row to start with
                    for c in n + 1..self.cols() {
                        if self.at(n, c).unwrap() != &R::zero() {
                            //swap column n and column c
                            let col_opp = ElementaryOpp::new_col_opp(ElementaryOppType::Swap(n, c));
                            col_opp.apply(&mut self);
                            col_opp.apply(&mut v);

                            break 'search_for_nonzero_element;
                        }
                    }
                    //search all the rows below row n
                    for r in n + 1..self.rows() {
                        for c in n..self.cols() {
                            if self.at(r, c).unwrap() != &R::zero() {
                                //swap column n and column c
                                let col_opp =
                                    ElementaryOpp::new_col_opp(ElementaryOppType::Swap(n, c));
                                col_opp.apply(&mut self);
                                col_opp.apply(&mut v);

                                //swap row n and row r
                                let row_opp =
                                    ElementaryOpp::new_col_opp(ElementaryOppType::Swap(n, r));
                                row_opp.apply(&mut self);
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
            let (unit, _assoc) = self.at(n, n).unwrap().factor_fav_assoc_ref();
            let row_opp = ElementaryOpp::new_row_opp(ElementaryOppType::UnitMul {
                row: n,
                unit: unit.inv().unwrap(),
            });
            row_opp.apply(&mut self);
            row_opp.apply(&mut u);

            let mut first = true;
            let mut all_divisible;
            'zero_first_row_and_column_loop: loop {
                //replace the first row (a0, a1, ..., ak) with (gcd, 0, ..., 0). Might mess up the first column in the process
                all_divisible = true;
                for c in n + 1..self.cols() {
                    let a = self.at(n, n).unwrap();
                    let b = self.at(n, c).unwrap();
                    match R::div_refs(b, a) {
                        Ok(q) => {
                            //b is a multiple of a
                            //replace (a, b) with (a, 0) by subtracting a multiple of a from b
                            let col_opp =
                                ElementaryOpp::new_col_opp(ElementaryOppType::AddRowMul {
                                    i: c,
                                    j: n,
                                    x: q.neg(),
                                });
                            col_opp.apply(&mut self);
                            col_opp.apply(&mut v);
                        }
                        Err(RingDivisionError::NotDivisible) => {
                            all_divisible = false;
                            //b is not a multiple of a
                            //replace (a, b) with (gcd, 0)
                            let (d, x, y) = R::xgcd(a.clone(), b.clone());
                            debug_assert_eq!(R::add(R::mul_refs(&x, a), R::mul_refs(&y, b)), d);
                            let col_opp = ElementaryOpp::new_col_opp(ElementaryOppType::TwoInv {
                                i: n,
                                j: c,
                                a: x,
                                b: y,
                                c: R::div(b.clone(), d.clone()).unwrap().neg(),
                                d: R::div(a.clone(), d.clone()).unwrap(),
                            });
                            col_opp.apply(&mut self);
                            col_opp.apply(&mut v);
                        }
                        Err(RingDivisionError::DivideByZero) => {
                            //swap a and b
                            //a=0 so this does have the effect of (a, b) -> (gcd(a, b), 0)
                            let col_opp = ElementaryOpp::new_col_opp(ElementaryOppType::Swap(n, c));
                            col_opp.apply(&mut self);
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
                for r in n + 1..self.rows() {
                    let a = self.at(n, n).unwrap();
                    let b = self.at(r, n).unwrap();
                    match R::div_refs(b, a) {
                        Ok(q) => {
                            //b is a multiple of a
                            //replace (a, b) with (a, 0) by subtracting a multiple of a from b
                            let col_opp =
                                ElementaryOpp::new_row_opp(ElementaryOppType::AddRowMul {
                                    i: r,
                                    j: n,
                                    x: q.neg(),
                                });
                            col_opp.apply(&mut self);
                            col_opp.apply(&mut u);
                        }
                        Err(RingDivisionError::NotDivisible) => {
                            all_divisible = false;
                            //b is not a multiple of a
                            //replace (a, b) with (gcd, 0)
                            let (d, x, y) = R::xgcd(a.clone(), b.clone());
                            debug_assert_eq!(R::add(R::mul_refs(&x, a), R::mul_refs(&y, b)), d);
                            let row_opp = ElementaryOpp::new_row_opp(ElementaryOppType::TwoInv {
                                i: n,
                                j: r,
                                a: x,
                                b: y,
                                c: R::div(b.clone(), d.clone()).unwrap().neg(),
                                d: R::div(a.clone(), d.clone()).unwrap(),
                            });
                            row_opp.apply(&mut self);
                            row_opp.apply(&mut u);
                        }
                        Err(RingDivisionError::DivideByZero) => {
                            //swap a and b
                            //a=0 so this does have the effect of (a, b) -> (gcd(a, b), 0)
                            let col_opp = ElementaryOpp::new_row_opp(ElementaryOppType::Swap(n, r));
                            col_opp.apply(&mut self);
                            col_opp.apply(&mut u);
                        }
                    }
                }
                if all_divisible {
                    break 'zero_first_row_and_column_loop;
                }
            }
            //now the first row and the first column are all zero except the top left element at (n, n) which is non-zero
            debug_assert_ne!(self.at(n, n).unwrap(), &R::zero());
            //some more fiddling is needed now to make the top left element divides everything else
            for r in n + 1..self.rows() {
                //row(n) = row(n) + row(r)
                let row_opp = ElementaryOpp::new_row_opp(ElementaryOppType::AddRowMul {
                    i: n,
                    j: r,
                    x: R::one(),
                });
                row_opp.apply(&mut self);
                row_opp.apply(&mut u);

                //row(n) goes from (g, a1, a2, ..., an) to (gcd, 0, 0, ..., 0)
                for c in n + 1..self.cols() {
                    let a = self.at(n, n).unwrap();
                    let b = self.at(n, c).unwrap();
                    //if a=0 and b=0 then there is nothing to do and the following step would fail when dividing by g=0
                    //b might not be a multiple of a
                    //replace (a, b) with (gcd, 0) to fix this
                    let (g, x, y) = R::xgcd(a.clone(), b.clone());
                    debug_assert_eq!(R::add(R::mul_refs(&x, a), R::mul_refs(&y, b)), g);
                    let col_opp = ElementaryOpp::new_col_opp(ElementaryOppType::TwoInv {
                        i: n,
                        j: c,
                        a: x,
                        b: y,
                        c: R::div(b.clone(), g.clone()).unwrap().neg(),
                        d: R::div(a.clone(), g.clone()).unwrap(),
                    });
                    col_opp.apply(&mut self);
                    col_opp.apply(&mut v);
                }

                //fix the first column
                for fix_r in n + 1..self.rows() {
                    let a = self.at(n, n).unwrap();
                    let b = self.at(fix_r, n).unwrap();
                    let q = R::div_refs(b, a).unwrap();
                    let col_opp = ElementaryOpp::new_row_opp(ElementaryOppType::AddRowMul {
                        i: fix_r,
                        j: n,
                        x: q.neg(),
                    });
                    col_opp.apply(&mut self);
                    col_opp.apply(&mut u);
                }
            }

            if self.at(n, n).unwrap() == &R::zero() {
                //the bottom right submatrix is all zero
                break 'inductive_loop;
            }
            n += 1;
        }

        (u, self, v, n)
    }
}

impl<R: EuclideanDomain + FavoriteAssociate> Matrix<R> {
    //if A:=self return (H, U, pivots) such that
    //H is in row reduced hermite normal form
    //U is invertible
    //H=UA
    //pivots[r] is the column of the rth pivot and pivots.len() == rank(A)
    pub fn row_reduced_hermite_algorithm(self) -> (Self, Self, Vec<usize>) {
        let (mut h, mut u, _u_det, pivs) = self.row_hermite_algorithm();

        for (pr, pc) in pivs.iter().enumerate() {
            for r in 0..pr {
                //reduce h[r, pc] so that it has norm less than h[pr, pc]
                let a = h.at(r, *pc).unwrap();
                let b = h.at(pr, *pc).unwrap();
                //a = b*q + r
                let q = R::quo_refs(a, b).unwrap();
                let row_opp = ElementaryOpp::new_row_opp(ElementaryOppType::AddRowMul {
                    i: r,
                    j: pr,
                    x: q.neg(),
                });
                row_opp.apply(&mut h);
                row_opp.apply(&mut u);
            }
        }

        (h, u, pivs)
    }

    pub fn col_reduced_hermite_algorithm(self) -> (Self, Self, Vec<usize>) {
        let (rh, ru, pivs) = self.transpose().row_reduced_hermite_algorithm();
        (rh.transpose(), ru.transpose(), pivs)
    }
}

impl<F: Field> Matrix<F> {
    pub fn presentation_matrix(self) -> Result<Matrix<Polynomial<F>>, MatOppErr> {
        let n = self.rows();
        if n != self.cols() {
            Err(MatOppErr::NotSquare)
        } else {
            Ok(Matrix::add(
                self.apply_map(|x| Polynomial::from(x.clone())),
                Matrix::diag(&(0..n).map(|_i| Polynomial::var()).collect()).neg(),
            )
            .unwrap())
        }
    }

    pub fn minimal_polynomial(self) -> Result<Polynomial<F>, MatOppErr> {
        match self.presentation_matrix() {
            Ok(pres_mat) => {
                let (_u, s, _v, k) = pres_mat.smith_algorithm();
                debug_assert!(k > 0); //cant be all zero becasue we are taking SNF of a non-zero matrix
                Ok(s.at(k - 1, k - 1).unwrap().clone())
            }
            Err(MatOppErr::NotSquare) => Err(MatOppErr::NotSquare),
            Err(_) => panic!(),
        }
    }

    pub fn characteristic_polynomial(self) -> Result<Polynomial<F>, MatOppErr> {
        match self.presentation_matrix() {
            Ok(pres_mat) => {
                let (_u, s, _v, k) = pres_mat.smith_algorithm();
                debug_assert!(k > 0); //cant be all zero becasue we are taking SNF of a non-zero matrix
                let mut char_poly = Polynomial::one();
                for i in 0..k {
                    char_poly.mul_mut(s.at(i, i).unwrap())
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

    use super::*;

    #[test]
    fn invariants() {
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
            assert_eq!(h, Matrix::<Integer>::mul_refs(&u, &a).unwrap());

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
                        let (_unit, assoc) = h.at(r, *pc).unwrap().clone().factor_fav_assoc();
                        assert_eq!(&assoc, h.at(r, *pc).unwrap());
                    } else {
                        assert!(h.at(r, *pc).unwrap().norm() < h.at(pr, *pc).unwrap().norm());
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
            assert_eq!(h, Matrix::<Integer>::mul_refs(&a, &u).unwrap());
            for (pc, pr) in pivs.iter().enumerate() {
                assert!(h.at(*pr, pc).unwrap() != &Integer::zero());
                for c in 0..h.cols() {
                    if c > pc {
                        assert_eq!(h.at(*pr, c).unwrap(), &Integer::zero());
                    } else if c == pc {
                        let (_unit, assoc) = h.at(*pr, c).unwrap().clone().factor_fav_assoc();
                        assert_eq!(&assoc, h.at(*pr, c).unwrap());
                    } else {
                        assert!(h.at(*pr, c).unwrap().norm() < h.at(*pr, pc).unwrap().norm());
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

            let (h, u, pivs) = a.clone().row_reduced_hermite_algorithm();

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

            let (h, u, pivs) = a.clone().row_reduced_hermite_algorithm();

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
            let (u, s, v, k) = a.clone().smith_algorithm();
            assert_eq!(
                s,
                Matrix::mul_refs(&Matrix::mul_refs(&u, &a).unwrap(), &v).unwrap()
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
            let (u, s, v, k) = a.clone().smith_algorithm();
            assert_eq!(
                s,
                Matrix::mul_refs(&Matrix::mul_refs(&u, &a).unwrap(), &v).unwrap()
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
            let (_u, _s, _v, k) = a.clone().smith_algorithm();
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
            let a = Matrix::from_rows(vec![
                vec![Rational::from(0), Rational::from(0), Rational::from(0)],
                vec![Rational::from(0), Rational::from(0), Rational::from(1)],
                vec![Rational::from(0), Rational::from(0), Rational::from(0)],
            ]);
            let min_p = a.clone().minimal_polynomial().unwrap();
            let char_p = a.clone().characteristic_polynomial().unwrap();
            assert_eq!(
                min_p,
                Polynomial::new(vec![
                    Rational::from(0),
                    Rational::from(0),
                    Rational::from(1)
                ])
            );
            assert_eq!(
                char_p,
                Polynomial::new(vec![
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

        assert_eq!(mat.clone().row_span().rank(), 2);
        assert_eq!(mat.clone().col_span().rank(), 2);
        assert_eq!(mat.clone().row_kernel().rank(), 2);
        assert_eq!(mat.clone().col_kernel().rank(), 1);
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
}
