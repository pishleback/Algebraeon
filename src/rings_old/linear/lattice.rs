use std::borrow::Borrow;
use std::fmt::Display;

use super::matrix::*;
use super::super::numbers::nzq::*;
use super::super::ring::*;

//return a metamatrix whose rows are a basis for the joint row span of all the passed metamatricies
fn metamatrix_row_sum<Ring: BezoutDomain, MetaMatT: Borrow<Matrix<Ring>>>(
    cols: usize,
    metamats: Vec<MetaMatT>,
) -> Matrix<Ring> {
    for metamat in &metamats {
        assert_eq!(metamat.borrow().cols(), cols);
    }
    let joined_metamat = Matrix::join_rows(cols, metamats);
    let (h, _u, _u_det, pivs) = Matrix::row_hermite_algorithm(joined_metamat);
    h.submatrix((0..pivs.len()).collect(), (0..cols).collect()) //return the top non-zero and linearly independent rows from h
}

//return a metamatrix whose rows are a basis for the intersection of the row spans of the passed metamatricies
fn metamatrix_row_intersection<Ring: BezoutDomain, MetaMatT: Borrow<Matrix<Ring>>>(
    cols: usize,
    mut metamat1: MetaMatT,
    mut metamat2: MetaMatT,
) -> Matrix<Ring> {
    assert_eq!(metamat1.borrow().cols(), cols);
    assert_eq!(metamat2.borrow().cols(), cols);
    //metamats should have linearly independent rows
    debug_assert_eq!(metamat1.borrow().clone().rank(), metamat1.borrow().rows());
    debug_assert_eq!(metamat2.borrow().clone().rank(), metamat2.borrow().rows());
    if metamat1.borrow().rows() > metamat2.borrow().rows() {
        //optimize for when we take linear combinations of rows from metamat1 rather than metamat2 later
        (metamat1, metamat2) = (metamat2, metamat1);
    }
    let joined_metamat = Matrix::join_rows(
        cols,
        vec![metamat1.borrow(), &Matrix::neg(metamat2.borrow().clone())],
    );
    //the row kernel of joined_metamat tells us about which linear combinations of rows of metamat1 are equal to which linear combinations of rows of metamat2
    let row_ker = joined_metamat.row_kernel();
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
        (0..row_ker.rank())
            .map(|i| {
                Matrix::mul_refs(
                    &LinearLattice::basis_matrix(&row_ker, i)
                        .submatrix(vec![0], (0..metamat1.borrow().rows()).collect()),
                    metamat1.borrow(),
                )
                .unwrap()
            })
            .collect(),
    )
}

#[derive(Debug, Clone)]
pub struct LinearLattice<Ring: BezoutDomain> {
    //matrix whose rows are a basis of the linear lattice
    //NOT necessarily in row hermite normal form
    metamatrix: Matrix<Ring>,
    //each row represents a matrix of this shape
    rows: usize,
    cols: usize,
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

fn mats_to_rows<Ring: ComRing, MatT: Borrow<Matrix<Ring>>>(
    rows: usize,
    cols: usize,
    mats: Vec<MatT>,
) -> Matrix<Ring> {
    for mat in &mats {
        assert_eq!(mat.borrow().rows(), rows);
        assert_eq!(mat.borrow().cols(), cols);
    }
    Matrix::construct(mats.len(), rows * cols, |r, c| {
        let (mr, mc) = idx_to_rc(rows, cols, c);
        mats[r].borrow().at(mr, mc).unwrap().clone()
    })
}

impl<Ring: BezoutDomain> LinearLattice<Ring> {
    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn check_invariants(&self) -> Result<(), &'static str> {
        if self.rows * self.cols != self.metamatrix.cols() {
            return Err("the number of colnums of the meta_matrix should be rows*cols");
        }
        if self.metamatrix.clone().rank() != self.metamatrix.rows() {
            return Err("the rows of meta_matrix should be linearly independent");
        }
        Ok(())
    }

    pub fn from_span<MatT: Borrow<Matrix<Ring>>>(
        rows: usize,
        cols: usize,
        mats: Vec<MatT>,
    ) -> LinearLattice<Ring> {
        let spanning_meta_matrix = mats_to_rows(rows, cols, mats);
        let (h, _u, _u_det, pivs) = spanning_meta_matrix.row_hermite_algorithm();
        let metamatrix = h.submatrix((0..pivs.len()).collect(), (0..rows * cols).collect());
        let lattice = LinearLattice {
            metamatrix,
            rows,
            cols,
        };
        debug_assert!(lattice.check_invariants().is_ok());
        lattice
    }

    pub fn from_basis<MatT: Borrow<Matrix<Ring>>>(
        rows: usize,
        cols: usize,
        mats: Vec<MatT>,
    ) -> Self {
        let metamatrix = mats_to_rows(rows, cols, mats);
        let lattice = LinearLattice {
            metamatrix: metamatrix,
            rows,
            cols,
        };
        debug_assert!(lattice.check_invariants().is_ok());
        lattice
    }

    pub fn rank(&self) -> usize {
        self.metamatrix.rows()
    }

    fn basis_row(&self, basis_num: usize) -> Matrix<Ring> {
        self.metamatrix
            .submatrix(vec![basis_num], (0..self.metamatrix.cols()).collect())
    }

    pub fn basis_matrix(&self, r: usize) -> Matrix<Ring> {
        if self.rank() <= r {
            panic!();
        }
        Matrix::construct(self.rows, self.cols, |mr, mc| {
            self.metamatrix
                .at(r, rc_to_idx(self.rows, self.cols, mr, mc))
                .unwrap()
                .clone()
        })
    }

    pub fn basis_matrices(&self) -> Vec<Matrix<Ring>> {
        (0..self.rank()).map(|r| self.basis_matrix(r)).collect()
    }

    pub fn basis_matrix_element<'b>(&'b self, basis_num: usize, r: usize, c: usize) -> &'b Ring {
        self.metamatrix
            .at(basis_num, rc_to_idx(self.rows, self.cols, r, c))
            .unwrap()
    }

    fn contains_row<MatT: Borrow<Matrix<Ring>>>(&self, mat_as_row: MatT) -> bool {
        match self.metamatrix.row_solve(mat_as_row) {
            Some(_taps) => true,
            None => false,
        }
    }

    pub fn contains_point<MatT: Borrow<Matrix<Ring>>>(&self, mat: MatT) -> bool {
        self.contains_row(mats_to_rows(self.rows, self.cols, vec![mat]))
    }

    //is lat a subset of self?
    pub fn contains_sublattice<LatT: Borrow<LinearLattice<Ring>>>(&self, sublat: LatT) -> bool {
        assert_eq!(self.metamatrix.cols(), sublat.borrow().metamatrix.cols());
        for basis_num in 0..sublat.borrow().rank() {
            if !self.contains_row(sublat.borrow().basis_row(basis_num)) {
                return false;
            }
        }
        true
    }

    pub fn eq(lat1: &Self, lat2: &Self) -> bool {
        lat1.contains_sublattice(lat2) && lat2.contains_sublattice(lat1)
    }

    pub fn sum<LatT: Borrow<LinearLattice<Ring>>>(
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> Self {
        for lat in &lats {
            assert_eq!(lat.borrow().rows, rows);
            assert_eq!(lat.borrow().cols, cols);
        }
        LinearLattice {
            rows,
            cols,
            metamatrix: metamatrix_row_sum(
                rows * cols,
                lats.iter().map(|lat| &lat.borrow().metamatrix).collect(),
            ),
        }
    }

    pub fn sum_pair<LatT: Borrow<LinearLattice<Ring>>>(
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> Self {
        Self::sum(rows, cols, vec![lat1, lat2])
    }

    pub fn intersect<LatT: Borrow<LinearLattice<Ring>>>(
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> Self {
        if lats.len() == 0 {
            LinearLattice {
                rows,
                cols,
                metamatrix: Matrix::ident(rows * cols),
            }
        } else if lats.len() == 1 {
            lats[0].borrow().clone()
        } else {
            let mut int_lat = Self::intersect_pair(rows, cols, lats[0].borrow(), lats[1].borrow());
            for i in 2..lats.len() {
                int_lat = Self::intersect_pair(rows, cols, &int_lat, lats[i].borrow());
            }
            int_lat
        }
    }

    pub fn intersect_pair<LatT: Borrow<LinearLattice<Ring>>>(
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> Self {
        assert_eq!(lat1.borrow().rows, rows);
        assert_eq!(lat1.borrow().cols, cols);
        assert_eq!(lat2.borrow().rows, rows);
        assert_eq!(lat2.borrow().cols, cols);
        let intersection_lattice = LinearLattice {
            rows,
            cols,
            metamatrix: metamatrix_row_intersection(
                rows * cols,
                &lat1.borrow().metamatrix,
                &lat2.borrow().metamatrix,
            ),
        };
        debug_assert!(intersection_lattice.check_invariants().is_ok());
        intersection_lattice
    }

    pub fn as_hyperplane_intersection(&self) -> Vec<Self> {
        //extend the basis of lat to a full basis
        let mut extended_basis = vec![];
        let (rows, cols) = (self.rows, self.cols);
        let mut extended_lat = self.clone();
        for i in 0..rows {
            for j in 0..cols {
                let v = Matrix::construct(rows, cols, |r, c| {
                    if i == r && j == c {
                        Ring::one()
                    } else {
                        Ring::zero()
                    }
                });
                let new_extended_lat = Self::sum_pair(
                    rows,
                    cols,
                    &extended_lat,
                    &Self::from_basis(rows, cols, vec![v.clone()]),
                );

                debug_assert!(
                    new_extended_lat.rank() == extended_lat.rank()
                        || new_extended_lat.rank() == extended_lat.rank() + 1
                );

                if new_extended_lat.rank() == extended_lat.rank() + 1 {
                    extended_lat = new_extended_lat;
                    extended_basis.push(v);
                }
            }
        }
        debug_assert_eq!(extended_lat.rank(), rows * cols);

        //now there is one hyperplane for each subset of extended_basis which omits one element
        (0..extended_basis.len())
            .map(|i| {
                Self::sum_pair(
                    rows,
                    cols,
                    self,
                    &Self::from_basis(
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

impl<R: BezoutDomain + Display> LinearLattice<R> {
    pub fn pprint(&self) {
        println!("Start Linear Lattice");
        for r in 0..self.metamatrix.rows() {
            self.basis_matrix(r).pprint();
        }
        println!("End Linear Lattice");
    }
}

#[derive(Debug, Clone)]
pub enum AffineLatticeElements<Ring: BezoutDomain> {
    Empty(),
    NonEmpty {
        offset: Matrix<Ring>,        //offset.rows == 1 and offset.cols == self.cols
        linlat: LinearLattice<Ring>, //linlat.rows == self.rows and linlat.cols == self.cols
    },
}

#[derive(Debug, Clone)]
pub struct AffineLattice<Ring: BezoutDomain> {
    rows: usize,
    cols: usize,
    elems: AffineLatticeElements<Ring>,
}

impl<Ring: BezoutDomain> AffineLattice<Ring> {
    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn elems(&self) -> &AffineLatticeElements<Ring> {
        &self.elems
    }

    pub fn check_invariants(&self) -> Result<(), &'static str> {
        match &self.elems {
            AffineLatticeElements::Empty() => {}
            AffineLatticeElements::NonEmpty { offset, linlat } => match linlat.check_invariants() {
                Ok(()) => {
                    if offset.rows() != self.rows {
                        return Err("offset rows doesnt match self rows");
                    }
                    if offset.cols() != self.cols {
                        return Err("offset columns doesnt match self columns");
                    }
                    if linlat.rows() != self.rows {
                        return Err("linlat rows doesnt match self rows");
                    }
                    if linlat.cols() != self.cols {
                        return Err("linlat columns doesnt match self columns");
                    }
                }
                Err(msg) => {
                    return Err(msg);
                }
            },
        }
        Ok(())
    }

    pub fn rank(&self) -> Option<usize> {
        match &self.elems {
            AffineLatticeElements::Empty() => None,
            AffineLatticeElements::NonEmpty {
                offset: _offset,

                linlat,
            } => Some(linlat.rank()),
        }
    }

    pub fn empty(rows: usize, cols: usize) -> Self {
        AffineLattice {
            rows,
            cols,
            elems: AffineLatticeElements::Empty(),
        }
    }

    pub fn to_offset_and_linear_lattice(self) -> Option<(Matrix<Ring>, LinearLattice<Ring>)> {
        match self.elems {
            AffineLatticeElements::Empty() => None,
            AffineLatticeElements::NonEmpty { offset, linlat } => Some((offset, linlat)),
        }
    }

    pub fn from_offset_and_linear_lattice(
        rows: usize,
        cols: usize,
        offset: Matrix<Ring>,
        linlat: LinearLattice<Ring>,
    ) -> AffineLattice<Ring> {
        assert_eq!(offset.rows(), rows);
        assert_eq!(offset.cols(), cols);
        assert_eq!(linlat.rows(), rows);
        assert_eq!(linlat.cols(), cols);
        let afflat = AffineLattice {
            rows,
            cols,
            elems: AffineLatticeElements::NonEmpty { offset, linlat },
        };
        debug_assert!(afflat.check_invariants().is_ok());
        afflat
    }

    pub fn contains_point<MatT: Borrow<Matrix<Ring>>>(&self, mat: MatT) -> bool {
        match &self.elems {
            AffineLatticeElements::Empty() => false,
            AffineLatticeElements::NonEmpty { offset, linlat } => linlat.contains_point(
                Matrix::add_ref(Matrix::neg(offset.clone()), mat.borrow()).unwrap(),
            ),
        }
    }

    //is other a subset of self?
    pub fn contains_sublattice<LatT: Borrow<AffineLattice<Ring>>>(&self, other: LatT) -> bool {
        match &other.borrow().elems {
            AffineLatticeElements::Empty() => true,
            AffineLatticeElements::NonEmpty {
                offset: other_offset,
                linlat: other_linlat,
            } => match &self.elems {
                AffineLatticeElements::Empty() => false,
                AffineLatticeElements::NonEmpty {
                    offset: _self_offset,
                    linlat: _self_linlat,
                } => {
                    for bn in 0..other_linlat.rank() {
                        if !self.contains_point(
                            Matrix::add_ref(other_linlat.basis_matrix(bn), other_offset).unwrap(),
                        ) {
                            return false;
                        }
                    }
                    true
                }
            },
        }
    }

    pub fn eq(lat1: &AffineLattice<Ring>, lat2: &AffineLattice<Ring>) -> bool {
        lat1.contains_sublattice(lat2) && lat2.contains_sublattice(lat1)
    }

    pub fn sum<LatT: Borrow<AffineLattice<Ring>>>(
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> AffineLattice<Ring> {
        let mut sum_offset = Matrix::zero(rows, cols);
        let mut sum_linlats = vec![];
        for lat in &lats {
            assert_eq!(lat.borrow().rows(), rows);
            assert_eq!(lat.borrow().cols(), cols);
            match &lat.borrow().elems {
                AffineLatticeElements::Empty() => {
                    return Self::empty(rows, cols);
                }
                AffineLatticeElements::NonEmpty { offset, linlat } => {
                    Matrix::add_mut(&mut sum_offset, offset).unwrap();
                    sum_linlats.push(linlat);
                }
            }
        }
        AffineLattice {
            rows: rows,
            cols: cols,
            elems: AffineLatticeElements::NonEmpty {
                offset: sum_offset,
                linlat: LinearLattice::sum(rows, cols, sum_linlats),
            },
        }
    }

    pub fn sum_pair<LatT: Borrow<AffineLattice<Ring>>>(
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> AffineLattice<Ring> {
        Self::sum(rows, cols, vec![lat1, lat2])
    }

    pub fn intersect<LatT: Borrow<AffineLattice<Ring>>>(
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> AffineLattice<Ring> {
        if lats.len() == 0 {
            AffineLattice {
                rows,
                cols,
                elems: AffineLatticeElements::NonEmpty {
                    offset: Matrix::zero(rows, cols),
                    linlat: LinearLattice {
                        rows,
                        cols,
                        metamatrix: Matrix::ident(rows * cols),
                    },
                },
            }
        } else if lats.len() == 1 {
            lats[0].borrow().clone()
        } else {
            let mut int_lat = Self::intersect_pair(rows, cols, lats[0].borrow(), lats[1].borrow());
            for i in 2..lats.len() {
                int_lat = Self::intersect_pair(rows, cols, &int_lat, lats[i].borrow());
            }
            int_lat
        }
    }

    fn offset_linlat_to_metamat(
        rows: usize,
        cols: usize,
        offset: &Matrix<Ring>,
        linlat: &LinearLattice<Ring>,
    ) -> Matrix<Ring> {
        let mut metamat = Matrix::zero(1 + linlat.rank(), 1 + rows * cols);
        *metamat.at_mut(0, 0).unwrap() = Ring::one();
        for idx in 0..rows * cols {
            let (r, c) = idx_to_rc(rows, cols, idx);
            println!("rows={} cols={} r={} c={} idx={}", rows, cols, r, c, idx);
            *metamat.at_mut(0, 1 + idx).unwrap() = offset.at(r, c).unwrap().clone();
        }
        for bn in 0..linlat.rank() {
            for idx in 0..rows * cols {
                let (r, c) = idx_to_rc(rows, cols, idx);
                *metamat.at_mut(0 + 1 + bn, 1 + idx).unwrap() =
                    linlat.basis_matrix_element(bn, r, c).clone();
            }
        }
        metamat
    }

    pub fn intersect_pair<LatT: Borrow<AffineLattice<Ring>>>(
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> AffineLattice<Ring> {
        assert_eq!(lat1.borrow().rows(), rows);
        assert_eq!(lat1.borrow().cols(), cols);
        assert_eq!(lat2.borrow().rows(), rows);
        assert_eq!(lat2.borrow().cols(), cols);
        match &lat1.borrow().elems {
            AffineLatticeElements::Empty() => Self::empty(rows, cols),
            AffineLatticeElements::NonEmpty {
                offset: offset1,
                linlat: linlat1,
            } => match &lat2.borrow().elems {
                AffineLatticeElements::Empty() => Self::empty(rows, cols),
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

                    let metamat1 = Self::offset_linlat_to_metamat(rows, cols, offset1, linlat1);
                    let metamat2 = Self::offset_linlat_to_metamat(rows, cols, offset2, linlat2);
                    let int_metamat =
                        metamatrix_row_intersection(1 + rows * cols, metamat1, metamat2);

                    if int_metamat.rows() == 0 {
                        //the hyperlattice is just the origin, so the coresponding affine lattice - the intersection with the plane (1, *, ..., *) - is empty.
                        Self::empty(rows, cols)
                    } else {
                        let (int_metamat_h, _u, _u_det, pivs) = int_metamat.row_hermite_algorithm();
                        // int_metamat_h.pprint();
                        if Ring::is_unit(int_metamat_h.at(0, 0).unwrap().clone()) {
                            debug_assert_eq!(int_metamat_h.at(0, 0).unwrap(), &Ring::one());
                        }
                        if int_metamat_h.at(0, 0).unwrap() == &Ring::one() {
                            let mut int_offset = Matrix::zero(rows, cols);
                            for idx in 0..rows * cols {
                                let (r, c) = idx_to_rc(rows, cols, idx);
                                *int_offset.at_mut(r, c).unwrap() =
                                    int_metamat_h.at(0, 1 + idx).unwrap().clone();
                            }
                            let int_basis_mats: Vec<Matrix<Ring>> = (0..pivs.len() - 1)
                                .map(|bn| {
                                    debug_assert_eq!(
                                        int_metamat_h.at(1 + bn, 0).unwrap(),
                                        &Ring::zero()
                                    );
                                    let mut basis_mat = Matrix::zero(rows, cols);
                                    for idx in 0..rows * cols {
                                        let (r, c) = idx_to_rc(rows, cols, idx);
                                        *basis_mat.at_mut(r, c).unwrap() =
                                            int_metamat_h.at(1 + bn, 1 + idx).unwrap().clone();
                                    }
                                    basis_mat
                                })
                                .collect();
                            // int_offset.pprint();
                            // for basis_mat in &int_basis_mats {
                            //     basis_mat.pprint();
                            // }
                            Self::from_offset_and_linear_lattice(
                                rows,
                                cols,
                                int_offset,
                                LinearLattice::from_basis(rows, cols, int_basis_mats),
                            )
                        } else {
                            //the hyperlattice does not intersect the plane (1, *, ..., *) because int_metamat_h(0, 0) is not a unit
                            Self::empty(rows, cols)
                        }
                    }
                }
            },
        }
    }
}

impl<R: BezoutDomain + Display> AffineLattice<R> {
    pub fn pprint(&self) {
        println!("Start Affine Lattice");
        match &self.elems {
            AffineLatticeElements::Empty() => println!("Empty"),
            AffineLatticeElements::NonEmpty { offset, linlat } => {
                println!("Offset");
                offset.pprint();
                linlat.pprint();
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
        lattice.check_invariants().unwrap();

        let lattice = LinearLattice {
            metamatrix: Matrix::from_rows(vec![
                vec![Integer::from(0), Integer::from(3), Integer::from(0)],
                vec![Integer::from(2), Integer::from(0), Integer::from(1)],
            ]),
            rows: 2,
            cols: 2,
        };
        assert!(lattice.check_invariants().is_err());

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
        assert!(lattice.check_invariants().is_err());
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
        let lattice = LinearLattice::from_span(
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
            lattice.contains_point(&Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(3)],
                vec![Integer::from(0), Integer::from(1)],
            ]))
        );

        assert_eq!(
            false,
            lattice.contains_point(&Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(4)],
                vec![Integer::from(0), Integer::from(1)],
            ]))
        );

        assert_eq!(
            false,
            lattice.contains_point(&Matrix::from_rows(vec![
                vec![Integer::from(2), Integer::from(3)],
                vec![Integer::from(1), Integer::from(1)],
            ]))
        );

        assert!(!LinearLattice::eq(
            &LinearLattice::from_span(
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
            &LinearLattice::from_span(
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
            a.pprint();
            println!("b");
            b.pprint();
            println!("a & b");
            let int =
                LinearLattice::intersect_pair(3, 1, a.clone().col_span(), b.clone().col_span());
            int.pprint();
            println!("a + b");
            let sum = LinearLattice::sum_pair(3, 1, a.clone().col_span(), b.clone().col_span());
            sum.pprint();

            assert!(LinearLattice::eq(&int, &c.col_span()));
            assert!(LinearLattice::eq(&sum, &d.col_span()));
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
            a.pprint();
            println!("b");
            b.pprint();
            println!("a & b");
            let int =
                LinearLattice::intersect_pair(3, 1, a.clone().col_span(), b.clone().col_span());
            int.pprint();
            println!("a + b");
            let sum = LinearLattice::sum_pair(3, 1, a.clone().col_span(), b.clone().col_span());
            sum.pprint();

            assert!(LinearLattice::eq(&int, &c.col_span()));
            assert!(LinearLattice::eq(&sum, &d.col_span()));
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

            let int =
                LinearLattice::intersect(4, 1, vec![a.col_span(), b.col_span(), c.col_span()]);

            assert!(LinearLattice::eq(
                &int,
                &Matrix::from_rows(vec![
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
                ])
                .col_span()
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
            a.pprint();
            println!("b");
            b.pprint();
            println!("a & b");
            let int =
                LinearLattice::intersect_pair(3, 1, a.clone().col_span(), b.clone().col_span());
            int.pprint();
            println!("a + b");
            let sum = LinearLattice::sum_pair(3, 1, a.clone().col_span(), b.clone().col_span());
            sum.pprint();

            assert!(LinearLattice::eq(&int, &c.col_span()));
            assert!(LinearLattice::eq(&sum, &d.col_span()));
        }
    }

    #[test]
    fn affine_lattice_invariants() {
        let afflat = AffineLattice::<Integer> {
            rows: 2,
            cols: 2,
            elems: AffineLatticeElements::Empty(),
        };
        assert!(afflat.check_invariants().is_ok());

        let afflat = AffineLattice {
            rows: 2,
            cols: 2,
            elems: AffineLatticeElements::NonEmpty {
                linlat: LinearLattice::from_basis(
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
        assert!(afflat.check_invariants().is_ok());
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

        let alat1 = a1.col_solution_lattice(y1);
        let alat2 = a2.col_solution_lattice(y2);

        alat1.pprint();
        println!();
        alat2.pprint();

        let alat3 = AffineLattice::sum(3, 1, vec![alat1, alat2]);
        println!();
        alat3.pprint();

        let expected_alat3 = AffineLattice::from_offset_and_linear_lattice(
            3,
            1,
            Matrix::from_rows(vec![
                vec![Integer::from(2)],
                vec![Integer::from(0)],
                vec![Integer::from(1)],
            ]),
            LinearLattice::from_span(
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

        assert!(AffineLattice::eq(&alat3, &expected_alat3));
    }
}
