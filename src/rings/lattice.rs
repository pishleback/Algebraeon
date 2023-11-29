use std::borrow::Borrow;

use super::matrix::*;
use super::nzq::*;
use super::ring::*;

pub const ZZ_LINLAT: LinearLatticeStructure<IntegerRing> = LinearLatticeStructure { ring: &ZZ };
pub const QQ_LINLAT: LinearLatticeStructure<RationalField> = LinearLatticeStructure { ring: &QQ };
pub const ZZ_AFFLAT: AffineLatticeStructure<IntegerRing> = AffineLatticeStructure { ring: &ZZ };
pub const QQ_AFFLAT: AffineLatticeStructure<RationalField> = AffineLatticeStructure { ring: &QQ };

//return a metamatrix whose rows are a basis for the joint row span of all the passed metamatricies
fn metamatrix_row_sum<'a, R: PrincipalIdealDomain, MetaMatT: Borrow<Matrix<R::ElemT>>>(
    ring: &'a R,
    cols: usize,
    metamats: Vec<MetaMatT>,
) -> Matrix<R::ElemT> {
    let mat_struct = MatrixStructure::new(ring);
    for metamat in &metamats {
        assert_eq!(metamat.borrow().cols(), cols);
    }
    let joined_metamat = Matrix::join_rows(cols, metamats);
    let (h, _u, _u_det, pivs) = mat_struct.row_hermite_algorithm(joined_metamat);
    h.submatrix((0..pivs.len()).collect(), (0..cols).collect()) //return the top non-zero and linearly independent rows from h
}

//return a metamatrix whose rows are a basis for the intersection of the row spans of the passed metamatricies
fn metamatrix_row_intersection<'a, R: PrincipalIdealDomain, MetaMatT: Borrow<Matrix<R::ElemT>>>(
    ring: &'a R,
    cols: usize,
    mut metamat1: MetaMatT,
    mut metamat2: MetaMatT,
) -> Matrix<R::ElemT> {
    let mat_struct = MatrixStructure::new(ring);
    assert_eq!(metamat1.borrow().cols(), cols);
    assert_eq!(metamat2.borrow().cols(), cols);
    //metamats should have linearly independent rows
    debug_assert_eq!(
        mat_struct.rank(metamat1.borrow().clone()),
        metamat1.borrow().rows()
    );
    debug_assert_eq!(
        mat_struct.rank(metamat2.borrow().clone()),
        metamat2.borrow().rows()
    );
    if metamat1.borrow().rows() > metamat2.borrow().rows() {
        //optimize for when we take linear combinations of rows from metamat1 rather than metamat2 later
        (metamat1, metamat2) = (metamat2, metamat1);
    }
    let joined_metamat = Matrix::join_rows(
        cols,
        vec![
            metamat1.borrow(),
            &mat_struct.neg(metamat2.borrow().clone()),
        ],
    );
    //the row kernel of joined_metamat tells us about which linear combinations of rows of metamat1 are equal to which linear combinations of rows of metamat2
    let row_ker = mat_struct.row_kernel(joined_metamat);
    //the rows in row_ker are in two halves
    //the first represents a linear combination of metamat1 rows
    //the second represents a linear combination of metamat2 rows
    //take without loss of generality all the linear combinations of metamat1 rows
    //the rows of the resulting matric are linearly independent
    //    because projection from (linear combinations of rows of metamat1 and metamat2) to (linear combinations of rows of metamat1) is injective
    //    because if a linear combination of rows of metamat1 and metamat2 is such that the metamat1 part is zero, then a linear combination of rows of metamat2 is zero
    //    but metamat2 is linearly independent, so the whole linear combination is zero
    Matrix::join_rows(
        cols,
        (0..LinearLatticeStructure::new(ring).rank(&row_ker))
            .map(|i| {
                mat_struct
                    .mul_refs(
                        &LinearLatticeStructure::new(ring)
                            .basis_matrix(&row_ker, i)
                            .submatrix(vec![0], (0..metamat1.borrow().rows()).collect()),
                        metamat1.borrow(),
                    )
                    .unwrap()
            })
            .collect(),
    )
}

#[derive(Debug, Clone)]
pub struct LinearLattice<ElemT: Clone> {
    //matrix whose rows are a basis of the linear lattice
    //NOT necessarily in row hermite normal form
    metamatrix: Matrix<ElemT>,
    //each row represents a matrix of this shape
    rows: usize,
    cols: usize,
}

impl<ElemT: Clone> LinearLattice<ElemT> {
    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }
}

//from matrix coords to meta row index
fn rc_to_idx(rows: usize, cols: usize, r: usize, c: usize) -> usize {
    if rows <= r || cols <= c {
        panic!();
    }
    c + r * cols
}

//from meta row index to matrix coords
fn idx_to_rc(rows: usize, cols: usize, idx: usize) -> (usize, usize) {
    if rows * cols <= idx {
        panic!();
    }
    (idx / cols, idx % cols)
}

fn mats_to_rows<ElemT: Clone, MatT: Borrow<Matrix<ElemT>>>(
    rows: usize,
    cols: usize,
    mats: Vec<MatT>,
) -> Matrix<ElemT> {
    for mat in &mats {
        assert_eq!(mat.borrow().rows(), rows);
        assert_eq!(mat.borrow().cols(), cols);
    }
    Matrix::construct(mats.len(), rows * cols, |r, c| {
        let (mr, mc) = idx_to_rc(rows, cols, c);
        mats[r].borrow().at(mr, mc).unwrap().clone()
    })
}

#[derive(Debug, Clone)]
pub struct LinearLatticeStructure<'a, R: PrincipalIdealDomain> {
    ring: &'a R,
}

impl<'a, R: PrincipalIdealDomain> LinearLatticeStructure<'a, R> {
    pub fn new(ring: &'a R) -> Self {
        Self { ring }
    }
}

impl<'a, R: PrincipalIdealDomain> LinearLatticeStructure<'a, R> {
    pub fn check_invariants(&self, lat: &LinearLattice<R::ElemT>) -> Result<(), &'static str> {
        if lat.rows * lat.cols != lat.metamatrix.cols() {
            return Err("the number of colnums of the meta_matrix should be rows*cols");
        }
        if MatrixStructure::new(self.ring).rank(lat.metamatrix.clone()) != lat.metamatrix.rows() {
            return Err("the rows of meta_matrix should be linearly independent");
        }
        Ok(())
    }

    pub fn from_span<MatT: Borrow<Matrix<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        mats: Vec<MatT>,
    ) -> LinearLattice<R::ElemT> {
        let spanning_meta_matrix = mats_to_rows(rows, cols, mats);
        let (h, _u, _u_det, pivs) =
            MatrixStructure::new(self.ring).row_hermite_algorithm(spanning_meta_matrix);
        let metamatrix = h.submatrix((0..pivs.len()).collect(), (0..rows * cols).collect());
        let lattice = LinearLattice {
            metamatrix,
            rows,
            cols,
        };
        debug_assert!(self.check_invariants(&lattice).is_ok());
        lattice
    }

    pub fn from_basis<MatT: Borrow<Matrix<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        mats: Vec<MatT>,
    ) -> LinearLattice<R::ElemT> {
        let metamatrix = mats_to_rows(rows, cols, mats);
        let lattice = LinearLattice {
            metamatrix: metamatrix,
            rows,
            cols,
        };
        debug_assert!(self.check_invariants(&lattice).is_ok());
        lattice
    }

    pub fn rank(&self, lat: &LinearLattice<R::ElemT>) -> usize {
        lat.metamatrix.rows()
    }

    fn basis_row(&self, lat: &LinearLattice<R::ElemT>, basis_num: usize) -> Matrix<R::ElemT> {
        lat.metamatrix
            .submatrix(vec![basis_num], (0..lat.metamatrix.cols()).collect())
    }

    pub fn basis_matrix(&self, lat: &LinearLattice<R::ElemT>, r: usize) -> Matrix<R::ElemT> {
        if self.rank(lat) <= r {
            panic!();
        }
        Matrix::construct(lat.rows, lat.cols, |mr, mc| {
            lat.metamatrix
                .at(r, rc_to_idx(lat.rows, lat.cols, mr, mc))
                .unwrap()
                .clone()
        })
    }

    pub fn basis_matrix_element<'b>(
        &self,
        lat: &'b LinearLattice<R::ElemT>,
        basis_num: usize,
        r: usize,
        c: usize,
    ) -> &'b R::ElemT {
        lat.metamatrix
            .at(basis_num, rc_to_idx(lat.rows, lat.cols, r, c))
            .unwrap()
    }

    fn contains_row<MatT: Borrow<Matrix<R::ElemT>>>(
        &self,
        lat: &LinearLattice<R::ElemT>,
        mat_as_row: MatT,
    ) -> bool {
        match MatrixStructure::new(self.ring).row_solve(&lat.metamatrix, mat_as_row) {
            Some(_taps) => true,
            None => false,
        }
    }

    pub fn contains_point<MatT: Borrow<Matrix<R::ElemT>>>(
        &self,
        lat: &LinearLattice<R::ElemT>,
        mat: MatT,
    ) -> bool {
        self.contains_row(lat, mats_to_rows(lat.rows, lat.cols, vec![mat]))
    }

    //is lat a subset of self?
    pub fn contains_sublattice<LatT: Borrow<LinearLattice<R::ElemT>>>(
        &self,
        lat: &LinearLattice<R::ElemT>,
        sublat: LatT,
    ) -> bool {
        assert_eq!(lat.metamatrix.cols(), sublat.borrow().metamatrix.cols());
        for basis_num in 0..self.rank(sublat.borrow()) {
            if !self.contains_row(lat, self.basis_row(sublat.borrow(), basis_num)) {
                return false;
            }
        }
        true
    }

    pub fn eq(&self, lat1: &LinearLattice<R::ElemT>, lat2: &LinearLattice<R::ElemT>) -> bool {
        self.contains_sublattice(lat1, lat2) && self.contains_sublattice(lat2, lat1)
    }

    pub fn sum<LatT: Borrow<LinearLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> LinearLattice<R::ElemT> {
        for lat in &lats {
            assert_eq!(lat.borrow().rows, rows);
            assert_eq!(lat.borrow().cols, cols);
        }
        LinearLattice {
            rows,
            cols,
            metamatrix: metamatrix_row_sum(
                self.ring,
                rows * cols,
                lats.iter().map(|lat| &lat.borrow().metamatrix).collect(),
            ),
        }
    }

    pub fn sum_pair<LatT: Borrow<LinearLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> LinearLattice<R::ElemT> {
        self.sum(rows, cols, vec![lat1, lat2])
    }

    pub fn intersect<LatT: Borrow<LinearLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> LinearLattice<R::ElemT> {
        if lats.len() == 0 {
            LinearLattice {
                rows,
                cols,
                metamatrix: MatrixStructure::new(self.ring).ident(rows * cols),
            }
        } else if lats.len() == 1 {
            lats[0].borrow().clone()
        } else {
            let mut int_lat = self.intersect_pair(rows, cols, lats[0].borrow(), lats[1].borrow());
            for i in 2..lats.len() {
                int_lat = self.intersect_pair(rows, cols, &int_lat, lats[i].borrow());
            }
            int_lat
        }
    }

    pub fn intersect_pair<LatT: Borrow<LinearLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> LinearLattice<R::ElemT> {
        assert_eq!(lat1.borrow().rows, rows);
        assert_eq!(lat1.borrow().cols, cols);
        assert_eq!(lat2.borrow().rows, rows);
        assert_eq!(lat2.borrow().cols, cols);
        let intersection_lattice = LinearLattice {
            rows,
            cols,
            metamatrix: metamatrix_row_intersection(
                self.ring,
                rows * cols,
                &lat1.borrow().metamatrix,
                &lat2.borrow().metamatrix,
            ),
        };
        debug_assert!(self.check_invariants(&intersection_lattice).is_ok());
        intersection_lattice
    }

    pub fn as_hyperplane_intersection(
        &self,
        lat: &LinearLattice<R::ElemT>,
    ) -> Vec<LinearLattice<R::ElemT>> {
        //extend the basis of lat to a full basis
        let mut extended_basis = vec![];
        let (rows, cols) = (lat.rows, lat.cols);
        let mut extended_lat = lat.clone();
        for i in 0..rows {
            for j in 0..cols {
                let v = Matrix::construct(rows, cols, |r, c| {
                    if i == r && j == c {
                        self.ring.one()
                    } else {
                        self.ring.zero()
                    }
                });
                let new_extended_lat = self.sum_pair(
                    rows,
                    cols,
                    &extended_lat,
                    &self.from_basis(rows, cols, vec![v.clone()]),
                );

                debug_assert!(
                    self.rank(&new_extended_lat) == self.rank(&extended_lat)
                        || self.rank(&new_extended_lat) == self.rank(&extended_lat) + 1
                );

                if self.rank(&new_extended_lat) == self.rank(&extended_lat) + 1 {
                    extended_lat = new_extended_lat;
                    extended_basis.push(v);
                }
            }
        }
        debug_assert_eq!(self.rank(&extended_lat), rows * cols);

        //now there is one hyperplane for each subset of extended_basis which omits one element
        (0..extended_basis.len())
            .map(|i| {
                self.sum_pair(
                    rows,
                    cols,
                    lat,
                    &self.from_basis(
                        rows,
                        cols,
                        extended_basis
                            .iter()
                            .enumerate()
                            .filter(|(j, _b)| *j != i)
                            .map(|(_j, b)| b)
                            .collect(),
                    ),
                )
            })
            .collect()
    }
}

impl<'a, R: PrincipalIdealDomain> LinearLatticeStructure<'a, R> {
    pub fn pprint(&self, lat: &LinearLattice<R::ElemT>) {
        println!("Start Linear Lattice");
        for r in 0..lat.metamatrix.rows() {
            MatrixStructure::new(self.ring).pprint(&self.basis_matrix(lat, r));
        }
        println!("End Linear Lattice");
    }
}

#[derive(Debug, Clone)]
pub enum AffineLatticeElements<ElemT: Clone> {
    Empty(),
    NonEmpty {
        offset: Matrix<ElemT>,        //offset.rows == 1 and offset.cols == self.cols
        linlat: LinearLattice<ElemT>, //linlat.rows == self.rows and linlat.cols == self.cols
    },
}

#[derive(Debug, Clone)]
pub struct AffineLattice<ElemT: Clone> {
    rows: usize,
    cols: usize,
    elems: AffineLatticeElements<ElemT>,
}

impl<ElemT: Clone> AffineLattice<ElemT> {
    pub fn elems(&self) -> &AffineLatticeElements<ElemT> {
        &self.elems
    }
}

pub struct AffineLatticeStructure<'a, R: PrincipalIdealDomain> {
    ring: &'a R,
}

impl<'a, R: PrincipalIdealDomain> AffineLatticeStructure<'a, R> {
    pub fn new(ring: &'a R) -> Self {
        Self { ring }
    }
}

impl<'a, R: PrincipalIdealDomain> AffineLatticeStructure<'a, R> {
    pub fn check_invariants(&self, lat: &AffineLattice<R::ElemT>) -> Result<(), &'static str> {
        match &lat.elems {
            AffineLatticeElements::Empty() => {}
            AffineLatticeElements::NonEmpty { offset, linlat } => {
                match LinearLatticeStructure::new(self.ring).check_invariants(linlat) {
                    Ok(()) => {
                        if offset.rows() != lat.rows {
                            return Err("offset rows doesnt match self rows");
                        }
                        if offset.cols() != lat.cols {
                            return Err("offset columns doesnt match self columns");
                        }
                        if linlat.rows() != lat.rows {
                            return Err("linlat rows doesnt match self rows");
                        }
                        if linlat.cols() != lat.cols {
                            return Err("linlat columns doesnt match self columns");
                        }
                    }
                    Err(msg) => {
                        return Err(msg);
                    }
                }
            }
        }
        Ok(())
    }

    pub fn rows(&self, lat: &AffineLattice<R::ElemT>) -> usize {
        lat.rows
    }

    pub fn cols(&self, lat: &AffineLattice<R::ElemT>) -> usize {
        lat.cols
    }

    pub fn rank(&self, lat: &AffineLattice<R::ElemT>) -> Option<usize> {
        match &lat.elems {
            AffineLatticeElements::Empty() => None,
            AffineLatticeElements::NonEmpty {
                offset: _offset,

                linlat,
            } => Some(LinearLatticeStructure::new(self.ring).rank(linlat)),
        }
    }

    pub fn empty(&self, rows: usize, cols: usize) -> AffineLattice<R::ElemT> {
        AffineLattice {
            rows,
            cols,
            elems: AffineLatticeElements::Empty(),
        }
    }

    pub fn from_offset_and_linear_lattice(
        &self,
        rows: usize,
        cols: usize,
        offset: Matrix<R::ElemT>,
        linlat: LinearLattice<R::ElemT>,
    ) -> AffineLattice<R::ElemT> {
        assert_eq!(offset.rows(), rows);
        assert_eq!(offset.cols(), cols);
        assert_eq!(linlat.rows(), rows);
        assert_eq!(linlat.cols(), cols);
        let afflat = AffineLattice {
            rows,
            cols,
            elems: AffineLatticeElements::NonEmpty { offset, linlat },
        };
        debug_assert!(self.check_invariants(&afflat).is_ok());
        afflat
    }

    pub fn contains_point<MatT: Borrow<Matrix<R::ElemT>>>(
        &self,
        lat: &AffineLattice<R::ElemT>,
        mat: MatT,
    ) -> bool {
        let mat_struct = MatrixStructure::new(self.ring);
        match &lat.elems {
            AffineLatticeElements::Empty() => false,
            AffineLatticeElements::NonEmpty { offset, linlat } => {
                LinearLatticeStructure::new(self.ring).contains_point(
                    linlat,
                    mat_struct
                        .add_ref(mat_struct.neg(offset.clone()), mat.borrow())
                        .unwrap(),
                )
            }
        }
    }

    //is other a subset of self?
    pub fn contains_sublattice<LatT: Borrow<AffineLattice<R::ElemT>>>(
        &self,
        lat: &AffineLattice<R::ElemT>,
        other: LatT,
    ) -> bool {
        let mat_struct = MatrixStructure::new(self.ring);
        let linlat_struct = LinearLatticeStructure::new(self.ring);
        match &other.borrow().elems {
            AffineLatticeElements::Empty() => true,
            AffineLatticeElements::NonEmpty {
                offset: other_offset,
                linlat: other_linlat,
            } => match &lat.elems {
                AffineLatticeElements::Empty() => false,
                AffineLatticeElements::NonEmpty {
                    offset: _self_offset,
                    linlat: _self_linlat,
                } => {
                    for bn in 0..linlat_struct.rank(other_linlat) {
                        if !self.contains_point(
                            lat,
                            mat_struct
                                .add_ref(linlat_struct.basis_matrix(other_linlat, bn), other_offset)
                                .unwrap(),
                        ) {
                            return false;
                        }
                    }
                    true
                }
            },
        }
    }

    pub fn eq(&self, lat1: &AffineLattice<R::ElemT>, lat2: &AffineLattice<R::ElemT>) -> bool {
        self.contains_sublattice(lat1, lat2) && self.contains_sublattice(lat2, lat1)
    }

    pub fn sum<LatT: Borrow<AffineLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> AffineLattice<R::ElemT> {
        let mut sum_offset = MatrixStructure::new(self.ring).zero(rows, cols);
        let mut sum_linlats = vec![];
        for lat in &lats {
            assert_eq!(self.rows(lat.borrow()), rows);
            assert_eq!(self.cols(lat.borrow()), cols);
            match &lat.borrow().elems {
                AffineLatticeElements::Empty() => {
                    return self.empty(rows, cols);
                }
                AffineLatticeElements::NonEmpty { offset, linlat } => {
                    MatrixStructure::new(self.ring)
                        .add_mut(&mut sum_offset, offset)
                        .unwrap();
                    sum_linlats.push(linlat);
                }
            }
        }
        AffineLattice {
            rows: rows,
            cols: cols,
            elems: AffineLatticeElements::NonEmpty {
                offset: sum_offset,
                linlat: LinearLatticeStructure::new(self.ring).sum(rows, cols, sum_linlats),
            },
        }
    }

    pub fn sum_pair<LatT: Borrow<AffineLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> AffineLattice<R::ElemT> {
        self.sum(rows, cols, vec![lat1, lat2])
    }

    pub fn intersect<LatT: Borrow<AffineLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> AffineLattice<R::ElemT> {
        if lats.len() == 0 {
            AffineLattice {
                rows,
                cols,
                elems: AffineLatticeElements::NonEmpty {
                    offset: MatrixStructure::new(self.ring).zero(rows, cols),
                    linlat: LinearLattice {
                        rows,
                        cols,
                        metamatrix: MatrixStructure::new(self.ring).ident(rows * cols),
                    },
                },
            }
        } else if lats.len() == 1 {
            lats[0].borrow().clone()
        } else {
            let mut int_lat = self.intersect_pair(rows, cols, lats[0].borrow(), lats[1].borrow());
            for i in 2..lats.len() {
                int_lat = self.intersect_pair(rows, cols, &int_lat, lats[i].borrow());
            }
            int_lat
        }
    }

    fn offset_linlat_to_metamat(
        &self,
        rows: usize,
        cols: usize,
        offset: &Matrix<R::ElemT>,
        linlat: &LinearLattice<R::ElemT>,
    ) -> Matrix<R::ElemT> {
        let mut metamat = MatrixStructure::new(self.ring).zero(
            1 + LinearLatticeStructure::new(self.ring).rank(linlat),
            1 + rows * cols,
        );
        *metamat.at_mut(0, 0).unwrap() = self.ring.one();
        for idx in 0..rows * cols {
            let (r, c) = idx_to_rc(rows, cols, idx);
            println!("rows={} cols={} r={} c={} idx={}", rows, cols, r, c, idx);
            *metamat.at_mut(0, 1 + idx).unwrap() = offset.at(r, c).unwrap().clone();
        }
        for bn in 0..LinearLatticeStructure::new(self.ring).rank(linlat) {
            for idx in 0..rows * cols {
                let (r, c) = idx_to_rc(rows, cols, idx);
                *metamat.at_mut(0 + 1 + bn, 1 + idx).unwrap() =
                    LinearLatticeStructure::new(self.ring)
                        .basis_matrix_element(linlat, bn, r, c)
                        .clone();
            }
        }
        metamat
    }

    pub fn intersect_pair<LatT: Borrow<AffineLattice<R::ElemT>>>(
        &self,
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> AffineLattice<R::ElemT> {
        assert_eq!(self.rows(lat1.borrow()), rows);
        assert_eq!(self.cols(lat1.borrow()), cols);
        assert_eq!(self.rows(lat2.borrow()), rows);
        assert_eq!(self.cols(lat2.borrow()), cols);
        match &lat1.borrow().elems {
            AffineLatticeElements::Empty() => self.empty(rows, cols),
            AffineLatticeElements::NonEmpty {
                offset: offset1,
                linlat: linlat1,
            } => match &lat2.borrow().elems {
                AffineLatticeElements::Empty() => self.empty(rows, cols),
                AffineLatticeElements::NonEmpty {
                    offset: offset2,
                    linlat: linlat2,
                } => {
                    //model an affine lattice as the intersection of a linear lattice (henceforth called the hyperlattice) living in one higher dimension with the plane (1, *, *, ..., *)
                    //take the intersection of the linear lattices and row reduce and get something like
                    // a * * * *
                    // 0 0 b * *
                    // 0 0 0 c *
                    //if a=1 then the rest of the top row is the offset and the bottom right submatrix is the basis of the linear lattice

                    let metamat1 = self.offset_linlat_to_metamat(rows, cols, offset1, linlat1);
                    let metamat2 = self.offset_linlat_to_metamat(rows, cols, offset2, linlat2);
                    let int_metamat =
                        metamatrix_row_intersection(self.ring, 1 + rows * cols, metamat1, metamat2);

                    if int_metamat.rows() == 0 {
                        //the hyperlattice is just the origin, so the coresponding affine lattice - the intersection with the plane (1, *, ..., *) - is empty.
                        self.empty(rows, cols)
                    } else {
                        let (int_metamat_h, _u, _u_det, pivs) =
                            MatrixStructure::new(self.ring).row_hermite_algorithm(int_metamat);
                        MatrixStructure::new(self.ring).pprint(&int_metamat_h);
                        if self.ring.is_unit(int_metamat_h.at(0, 0).unwrap().clone()) {
                            debug_assert!(self
                                .ring
                                .equal(int_metamat_h.at(0, 0).unwrap(), &self.ring.one()));
                        }
                        if self
                            .ring
                            .equal(int_metamat_h.at(0, 0).unwrap(), &self.ring.one())
                        {
                            let mut int_offset = MatrixStructure::new(self.ring).zero(rows, cols);
                            for idx in 0..rows * cols {
                                let (r, c) = idx_to_rc(rows, cols, idx);
                                *int_offset.at_mut(r, c).unwrap() =
                                    int_metamat_h.at(0, 1 + idx).unwrap().clone();
                            }
                            let int_basis_mats = (0..pivs.len() - 1)
                                .map(|bn| {
                                    debug_assert!(self.ring.equal(
                                        int_metamat_h.at(1 + bn, 0).unwrap(),
                                        &self.ring.zero()
                                    ));
                                    let mut basis_mat =
                                        MatrixStructure::new(self.ring).zero(rows, cols);
                                    for idx in 0..rows * cols {
                                        let (r, c) = idx_to_rc(rows, cols, idx);
                                        *basis_mat.at_mut(r, c).unwrap() =
                                            int_metamat_h.at(1 + bn, 1 + idx).unwrap().clone();
                                    }
                                    basis_mat
                                })
                                .collect();
                            MatrixStructure::new(self.ring).pprint(&int_offset);
                            for basis_mat in &int_basis_mats {
                                MatrixStructure::new(self.ring).pprint(basis_mat);
                            }
                            self.from_offset_and_linear_lattice(
                                rows,
                                cols,
                                int_offset,
                                LinearLatticeStructure::new(self.ring).from_basis(
                                    rows,
                                    cols,
                                    int_basis_mats,
                                ),
                            )
                        } else {
                            //the hyperlattice does not intersect the plane (1, *, ..., *) because int_metamat_h(0, 0) is not a unit
                            self.empty(rows, cols)
                        }
                    }
                }
            },
        }
    }
}

impl<'a, R: PrincipalIdealDomain> AffineLatticeStructure<'a, R> {
    pub fn pprint(&self, lat: &AffineLattice<R::ElemT>) {
        println!("Start Affine Lattice");
        match &lat.elems {
            AffineLatticeElements::Empty() => println!("Empty"),
            AffineLatticeElements::NonEmpty { offset, linlat } => {
                println!("Offset");
                MatrixStructure::new(self.ring).pprint(&offset);
                LinearLatticeStructure::new(self.ring).pprint(&linlat);
            }
        }
        println!("End Affine Lattice");
    }
}

#[cfg(test)]
mod tests {
    use malachite_nz::integer::Integer;

    use super::*;

    #[test]
    fn linear_lattice_invariant() {
        let lattice = LinearLattice {
            metamatrix: Matrix::from_rows(vec![
                vec![
                    Integer::from(0),
                    Integer::from(3),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(2),
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                ],
            ]),
            rows: 2,
            cols: 2,
        };
        ZZ_LINLAT.check_invariants(&lattice).unwrap();

        let lattice = LinearLattice {
            metamatrix: Matrix::from_rows(vec![
                vec![Integer::from(0), Integer::from(3), Integer::from(0)],
                vec![Integer::from(2), Integer::from(0), Integer::from(1)],
            ]),
            rows: 2,
            cols: 2,
        };
        assert!(ZZ_LINLAT.check_invariants(&lattice).is_err());

        let lattice = LinearLattice {
            metamatrix: Matrix::from_rows(vec![
                vec![
                    Integer::from(6),
                    Integer::from(0),
                    Integer::from(3),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(2),
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                ],
            ]),
            rows: 2,
            cols: 2,
        };
        assert!(ZZ_LINLAT.check_invariants(&lattice).is_err());
    }

    #[test]
    fn index_conversions() {
        let rows = 5;
        let cols = 7;

        for idx in 0..rows * cols {
            let (r, c) = idx_to_rc(rows, cols, idx);
            assert_eq!(idx, rc_to_idx(rows, cols, r, c));
        }
        for r in 0..rows {
            for c in 0..cols {
                let idx = rc_to_idx(rows, cols, r, c);
                assert_eq!((r, c), idx_to_rc(rows, cols, idx));
            }
        }
    }

    #[test]
    fn containment() {
        let lattice = ZZ_LINLAT.from_span(
            2,
            2,
            vec![
                &Matrix::from_rows(vec![
                    vec![Integer::from(0), Integer::from(3)],
                    vec![Integer::from(0), Integer::from(0)],
                ]),
                &Matrix::from_rows(vec![
                    vec![Integer::from(2), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(1)],
                ]),
            ],
        );

        assert_eq!(
            true,
            ZZ_LINLAT.contains_point(
                &lattice,
                &Matrix::from_rows(vec![
                    vec![Integer::from(2), Integer::from(3)],
                    vec![Integer::from(0), Integer::from(1)],
                ])
            )
        );

        assert_eq!(
            false,
            ZZ_LINLAT.contains_point(
                &lattice,
                &Matrix::from_rows(vec![
                    vec![Integer::from(2), Integer::from(4)],
                    vec![Integer::from(0), Integer::from(1)],
                ])
            )
        );

        assert_eq!(
            false,
            ZZ_LINLAT.contains_point(
                &lattice,
                &Matrix::from_rows(vec![
                    vec![Integer::from(2), Integer::from(3)],
                    vec![Integer::from(1), Integer::from(1)],
                ])
            )
        );

        assert!(!ZZ_LINLAT.eq(
            &ZZ_LINLAT.from_span(
                2,
                3,
                vec![
                    &Matrix::from_rows(vec![
                        vec![Integer::from(0), Integer::from(2), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                    &Matrix::from_rows(vec![
                        vec![Integer::from(0), Integer::from(4), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                ],
            ),
            &ZZ_LINLAT.from_span(
                2,
                3,
                vec![
                    &Matrix::from_rows(vec![
                        vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                    &Matrix::from_rows(vec![
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                ],
            )
        ));
    }

    #[test]
    fn linear_lattice_sum_and_intersection() {
        {
            //standard basis sum and intersection
            let a = Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
            ]);

            let b = Matrix::from_rows(vec![
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(1)],
            ]);

            let c = Matrix::from_rows(vec![
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(0)],
            ]);

            let d = Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(1)],
            ]);

            println!();
            println!("a");
            ZZ_MAT.pprint(&a);
            println!("b");
            ZZ_MAT.pprint(&b);
            println!("a & b");
            let int = ZZ_LINLAT.intersect_pair(
                3,
                1,
                ZZ_MAT.col_span(a.clone()),
                ZZ_MAT.col_span(b.clone()),
            );
            ZZ_LINLAT.pprint(&int);
            println!("a + b");
            let sum =
                ZZ_LINLAT.sum_pair(3, 1, ZZ_MAT.col_span(a.clone()), ZZ_MAT.col_span(b.clone()));
            ZZ_LINLAT.pprint(&sum);

            assert!(ZZ_LINLAT.eq(&int, &ZZ_MAT.col_span(c)));
            assert!(ZZ_LINLAT.eq(&sum, &ZZ_MAT.col_span(d)));
        }

        {
            //sum and intersection as gcd and lcm
            let a = Matrix::from_rows(vec![
                vec![Integer::from(3), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(5), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(7)],
            ]);

            let b = Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(4), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(14)],
            ]);

            let c = Matrix::from_rows(vec![
                vec![Integer::from(6), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(20), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(14)],
            ]);

            let d = Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(0), Integer::from(1), Integer::from(0)],
                vec![Integer::from(0), Integer::from(0), Integer::from(7)],
            ]);

            println!();
            println!("a");
            ZZ_MAT.pprint(&a);
            println!("b");
            ZZ_MAT.pprint(&b);
            println!("a & b");
            let int = ZZ_LINLAT.intersect_pair(
                3,
                1,
                ZZ_MAT.col_span(a.clone()),
                ZZ_MAT.col_span(b.clone()),
            );
            ZZ_LINLAT.pprint(&int);
            println!("a + b");
            let sum =
                ZZ_LINLAT.sum_pair(3, 1, ZZ_MAT.col_span(a.clone()), ZZ_MAT.col_span(b.clone()));
            ZZ_LINLAT.pprint(&sum);

            assert!(ZZ_LINLAT.eq(&int, &ZZ_MAT.col_span(c)));
            assert!(ZZ_LINLAT.eq(&sum, &ZZ_MAT.col_span(d)));
        }

        {
            //triple intersection
            let a = Matrix::from_rows(vec![
                vec![
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
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                ],
            ]);

            let b = Matrix::from_rows(vec![
                vec![
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(1),
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

            let c = Matrix::from_rows(vec![
                vec![
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                    Integer::from(0),
                ],
            ]);

            let int = ZZ_LINLAT.intersect(
                4,
                1,
                vec![ZZ_MAT.col_span(a), ZZ_MAT.col_span(b), ZZ_MAT.col_span(c)],
            );

            assert!(ZZ_LINLAT.eq(
                &int,
                &ZZ_MAT.col_span(Matrix::from_rows(vec![
                    vec![
                        Integer::from(1),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                    vec![
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0),
                        Integer::from(0)
                    ],
                ]))
            ));
        }

        {
            //complex example
            let a = Matrix::from_rows(vec![
                vec![Integer::from(3), Integer::from(9), Integer::from(27)],
                vec![Integer::from(-4), Integer::from(6), Integer::from(-100)],
                vec![Integer::from(2), Integer::from(8), Integer::from(7)],
            ]);

            let b = Matrix::from_rows(vec![
                vec![Integer::from(12), Integer::from(-1), Integer::from(18)],
                vec![Integer::from(-5), Integer::from(12), Integer::from(-24)],
                vec![Integer::from(1), Integer::from(2), Integer::from(14)],
            ]);

            let c = Matrix::from_rows(vec![
                vec![
                    Integer::from(21),
                    Integer::from(3852),
                    Integer::from(3315300),
                ],
                vec![
                    Integer::from(-252),
                    Integer::from(-46214),
                    Integer::from(-39775000),
                ],
                vec![
                    Integer::from(-42),
                    Integer::from(-4454),
                    Integer::from(-3833450),
                ],
            ]);

            let d = Matrix::from_rows(vec![
                vec![Integer::from(1), Integer::from(0), Integer::from(0)],
                vec![Integer::from(-12), Integer::from(1), Integer::from(0)],
                vec![Integer::from(-2), Integer::from(325), Integer::from(1)],
            ]);

            println!();
            println!("a");
            ZZ_MAT.pprint(&a);
            println!("b");
            ZZ_MAT.pprint(&b);
            println!("a & b");
            let int = ZZ_LINLAT.intersect_pair(
                3,
                1,
                ZZ_MAT.col_span(a.clone()),
                ZZ_MAT.col_span(b.clone()),
            );
            ZZ_LINLAT.pprint(&int);
            println!("a + b");
            let sum =
                ZZ_LINLAT.sum_pair(3, 1, ZZ_MAT.col_span(a.clone()), ZZ_MAT.col_span(b.clone()));
            ZZ_LINLAT.pprint(&sum);

            assert!(ZZ_LINLAT.eq(&int, &ZZ_MAT.col_span(c)));
            assert!(ZZ_LINLAT.eq(&sum, &ZZ_MAT.col_span(d)));
        }
    }

    #[test]
    fn affine_lattice_invariants() {
        let afflat = AffineLattice {
            rows: 2,
            cols: 2,
            elems: AffineLatticeElements::Empty(),
        };
        assert!(ZZ_AFFLAT.check_invariants(&afflat).is_ok());

        let afflat = AffineLattice {
            rows: 2,
            cols: 2,
            elems: AffineLatticeElements::NonEmpty {
                linlat: ZZ_LINLAT.from_basis(
                    2,
                    2,
                    vec![
                        Matrix::from_rows(vec![
                            vec![Integer::from(1), Integer::from(0)],
                            vec![Integer::from(0), Integer::from(1)],
                        ]),
                        Matrix::from_rows(vec![
                            vec![Integer::from(0), Integer::from(1)],
                            vec![Integer::from(1), Integer::from(0)],
                        ]),
                    ],
                ),
                offset: Matrix::from_rows(vec![
                    vec![Integer::from(1), Integer::from(0)],
                    vec![Integer::from(0), Integer::from(0)],
                ]),
            },
        };
        assert!(ZZ_AFFLAT.check_invariants(&afflat).is_ok());
    }

    fn affine_lattice_sum_and_intersection() {
        let a1 = Matrix::from_rows(vec![
            vec![Integer::from(3), Integer::from(1), Integer::from(0)],
            vec![Integer::from(3), Integer::from(1), Integer::from(0)],
            vec![Integer::from(3), Integer::from(1), Integer::from(1)],
        ]);

        let y1 = Matrix::from_rows(vec![
            vec![Integer::from(1)],
            vec![Integer::from(1)],
            vec![Integer::from(1)],
        ]);

        let a2 = Matrix::from_rows(vec![
            vec![Integer::from(3), Integer::from(5), Integer::from(0)],
            vec![Integer::from(3), Integer::from(5), Integer::from(0)],
            vec![Integer::from(0), Integer::from(0), Integer::from(1)],
        ]);

        let y2 = Matrix::from_rows(vec![
            vec![Integer::from(1)],
            vec![Integer::from(1)],
            vec![Integer::from(1)],
        ]);

        let alat1 = ZZ_MAT.col_solution_lattice(&a1, y1);
        let alat2 = ZZ_MAT.col_solution_lattice(&a2, y2);

        ZZ_AFFLAT.pprint(&alat1);
        println!();
        ZZ_AFFLAT.pprint(&alat2);

        let alat3 = ZZ_AFFLAT.sum(3, 1, vec![alat1, alat2]);
        println!();
        ZZ_AFFLAT.pprint(&alat3);

        let expected_alat3 = ZZ_AFFLAT.from_offset_and_linear_lattice(
            3,
            1,
            Matrix::from_rows(vec![
                vec![Integer::from(2)],
                vec![Integer::from(0)],
                vec![Integer::from(1)],
            ]),
            ZZ_LINLAT.from_span(
                3,
                1,
                vec![
                    Matrix::from_rows(vec![
                        vec![Integer::from(1)],
                        vec![Integer::from(-3)],
                        vec![Integer::from(0)],
                    ]),
                    Matrix::from_rows(vec![
                        vec![Integer::from(5)],
                        vec![Integer::from(-3)],
                        vec![Integer::from(0)],
                    ]),
                ],
            ),
        );

        assert!(ZZ_AFFLAT.eq(&alat3, &expected_alat3));
    }
}
