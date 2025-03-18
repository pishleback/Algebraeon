use std::borrow::Borrow;
use std::rc::Rc;

use super::super::structure::structure::*;
use super::matrix::*;
use algebraeon_sets::structure::*;

//return a metamatrix whose rows are a basis for the joint row span of all the passed metamatricies
fn metamatrix_row_sum<RS: BezoutDomainStructure, MetaMatT: Borrow<Matrix<RS::Set>>>(
    ring: Rc<RS>,
    cols: usize,
    metamats: Vec<MetaMatT>,
) -> Matrix<RS::Set> {
    let mat_struct = MatrixStructure::new(ring);
    for metamat in &metamats {
        assert_eq!(metamat.borrow().cols(), cols);
    }
    let joined_metamat = Matrix::join_rows(cols, metamats);
    let (h, _u, _u_det, pivs) = mat_struct.row_hermite_algorithm(joined_metamat);
    h.submatrix((0..pivs.len()).collect(), (0..cols).collect()) //return the top non-zero and linearly independent rows from h
}

//return a metamatrix whose rows are a basis for the intersection of the row spans of the passed metamatricies
fn metamatrix_row_intersection<RS: BezoutDomainStructure, MetaMatT: Borrow<Matrix<RS::Set>>>(
    ring: Rc<RS>,
    cols: usize,
    mut metamat1: MetaMatT,
    mut metamat2: MetaMatT,
) -> Matrix<RS::Set> {
    let mat_struct = MatrixStructure::new(ring.clone());
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
        (0..LinearLatticeStructure::new(ring.clone()).rank(&row_ker))
            .map(|i| {
                mat_struct
                    .mul(
                        &LinearLatticeStructure::new(ring.clone())
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
pub struct LinearLattice<Set: Clone> {
    //matrix whose rows are a basis of the linear lattice
    //NOT necessarily in row hermite normal form
    metamatrix: Matrix<Set>,
    //each row represents a matrix of this shape
    rows: usize,
    cols: usize,
}

impl<Set: Clone> LinearLattice<Set> {
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

fn mats_to_rows<Set: Clone, MatT: Borrow<Matrix<Set>>>(
    rows: usize,
    cols: usize,
    mats: Vec<MatT>,
) -> Matrix<Set> {
    for mat in &mats {
        assert_eq!(mat.borrow().rows(), rows);
        assert_eq!(mat.borrow().cols(), cols);
    }
    Matrix::construct(mats.len(), rows * cols, |r, c| {
        let (mr, mc) = idx_to_rc(rows, cols, c);
        mats[r].borrow().at(mr, mc).unwrap().clone()
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LinearLatticeStructure<RS: BezoutDomainStructure> {
    ring: Rc<RS>,
}

impl<RS: BezoutDomainStructure> Structure for LinearLatticeStructure<RS> {
    type Set = LinearLattice<RS::Set>;
}

impl<RS: BezoutDomainStructure> LinearLatticeStructure<RS> {
    pub fn new(ring: Rc<RS>) -> Self {
        Self { ring }
    }
}

impl<RS: BezoutDomainStructure> LinearLatticeStructure<RS> {
    pub fn check_invariants(&self, lat: &LinearLattice<RS::Set>) -> Result<(), &'static str> {
        if lat.rows * lat.cols != lat.metamatrix.cols() {
            return Err("the number of colnums of the meta_matrix should be rows*cols");
        }
        if MatrixStructure::new(self.ring.clone()).rank(lat.metamatrix.clone())
            != lat.metamatrix.rows()
        {
            return Err("the rows of meta_matrix should be linearly independent");
        }
        Ok(())
    }

    pub fn from_span<MatT: Borrow<Matrix<RS::Set>>>(
        &self,
        rows: usize,
        cols: usize,
        mats: Vec<MatT>,
    ) -> LinearLattice<RS::Set> {
        let spanning_meta_matrix = mats_to_rows(rows, cols, mats);
        let (h, _u, _u_det, pivs) =
            MatrixStructure::new(self.ring.clone()).row_hermite_algorithm(spanning_meta_matrix);
        let metamatrix = h.submatrix((0..pivs.len()).collect(), (0..rows * cols).collect());
        let lattice = LinearLattice {
            metamatrix,
            rows,
            cols,
        };
        debug_assert!(self.check_invariants(&lattice).is_ok());
        lattice
    }

    pub fn transpose(&self, lat: &LinearLattice<RS::Set>) -> LinearLattice<RS::Set> {
        LinearLattice {
            metamatrix: mats_to_rows(
                lat.cols,
                lat.rows,
                self.basis_matrices(lat)
                    .into_iter()
                    .map(|mat| mat.transpose())
                    .collect(),
            ),
            rows: lat.cols,
            cols: lat.rows,
        }
    }

    pub fn from_basis<MatT: Borrow<Matrix<RS::Set>>>(
        &self,
        rows: usize,
        cols: usize,
        mats: Vec<MatT>,
    ) -> LinearLattice<RS::Set> {
        let metamatrix = mats_to_rows(rows, cols, mats);
        let lattice = LinearLattice {
            metamatrix: metamatrix,
            rows,
            cols,
        };
        debug_assert!(self.check_invariants(&lattice).is_ok());
        lattice
    }

    pub fn zero(&self, rows: usize, cols: usize) -> LinearLattice<RS::Set> {
        self.from_basis::<Matrix<RS::Set>>(rows, cols, vec![])
    }

    pub fn rank(&self, lat: &LinearLattice<RS::Set>) -> usize {
        lat.metamatrix.rows()
    }

    pub fn take_nonzero_point(&self, lat: &LinearLattice<RS::Set>) -> Option<Matrix<RS::Set>> {
        if self.rank(&lat) == 0 {
            None
        } else {
            Some(self.basis_matrix(&lat, 0))
        }
    }

    fn basis_row(&self, lat: &LinearLattice<RS::Set>, basis_num: usize) -> Matrix<RS::Set> {
        lat.metamatrix
            .submatrix(vec![basis_num], (0..lat.metamatrix.cols()).collect())
    }

    pub fn basis_matrix(&self, lat: &LinearLattice<RS::Set>, r: usize) -> Matrix<RS::Set> {
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

    pub fn basis_matrices(&self, lat: &LinearLattice<RS::Set>) -> Vec<Matrix<RS::Set>> {
        (0..self.rank(lat))
            .map(|r| self.basis_matrix(lat, r))
            .collect()
    }

    pub fn basis_matrix_element<'b>(
        &self,
        lat: &'b LinearLattice<RS::Set>,
        basis_num: usize,
        r: usize,
        c: usize,
    ) -> &'b RS::Set {
        lat.metamatrix
            .at(basis_num, rc_to_idx(lat.rows, lat.cols, r, c))
            .unwrap()
    }

    fn contains_row(
        &self,
        lat: &LinearLattice<RS::Set>,
        mat_as_row: impl Borrow<Matrix<RS::Set>>,
    ) -> bool {
        match MatrixStructure::new(self.ring.clone()).row_solve(&lat.metamatrix, mat_as_row) {
            Some(_taps) => true,
            None => false,
        }
    }

    pub fn contains_point<MatT: Borrow<Matrix<RS::Set>>>(
        &self,
        lat: &LinearLattice<RS::Set>,
        mat: MatT,
    ) -> bool {
        self.contains_row(lat, mats_to_rows(lat.rows, lat.cols, vec![mat]))
    }

    //is lat a subset of self?
    pub fn contains_sublattice(
        &self,
        lat: &LinearLattice<RS::Set>,
        sublat: impl Borrow<LinearLattice<RS::Set>>,
    ) -> bool {
        assert_eq!(lat.metamatrix.cols(), sublat.borrow().metamatrix.cols());
        for basis_num in 0..self.rank(sublat.borrow()) {
            if !self.contains_row(lat, self.basis_row(sublat.borrow(), basis_num)) {
                return false;
            }
        }
        true
    }

    pub fn equal(&self, lat1: &LinearLattice<RS::Set>, lat2: &LinearLattice<RS::Set>) -> bool {
        self.contains_sublattice(lat1, lat2) && self.contains_sublattice(lat2, lat1)
    }

    //given a contained in b, find rank(b) - rank(a) basis vectors needed to extend a to b
    pub fn extension_basis(
        &self,
        lat_a: &LinearLattice<RS::Set>,
        lat_b: &LinearLattice<RS::Set>,
    ) -> Vec<Matrix<RS::Set>> {
        //https://math.stackexchange.com/questions/2554408/how-to-find-the-basis-of-a-quotient-space

        let rows = lat_a.rows();
        let cols = lat_b.cols();
        assert_eq!(rows, lat_b.rows());
        assert_eq!(cols, lat_b.cols());
        debug_assert!(self.contains_sublattice(lat_b, lat_a));

        let n = rows * cols;
        // form matrix of all vectors from [other self]
        // row reduce and get pivots - take cols from orig to form basies of the quotient space

        let mat_structure = MatrixStructure::new(self.ring.clone());
        let metamatrix = Matrix::join_rows(n, vec![&lat_a.metamatrix, &lat_b.metamatrix]);
        let (_h, _u, _u_det, pivs) = mat_structure.col_hermite_algorithm(metamatrix.clone());

        let mut extension_basis = vec![];
        for r in pivs {
            //dont take vectors which form a basis of lat_a
            if r >= self.rank(lat_a) {
                extension_basis.push(Matrix::construct(rows, cols, |mr, mc| {
                    metamatrix
                        .at(r, rc_to_idx(rows, cols, mr, mc))
                        .unwrap()
                        .clone()
                }));
            }
        }

        debug_assert_eq!(self.rank(&lat_a) + extension_basis.len(), self.rank(&lat_b));

        extension_basis
    }

    pub fn sum<LatT: Borrow<LinearLattice<RS::Set>>>(
        &self,
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> LinearLattice<RS::Set> {
        for lat in &lats {
            assert_eq!(lat.borrow().rows, rows);
            assert_eq!(lat.borrow().cols, cols);
        }
        LinearLattice {
            rows,
            cols,
            metamatrix: metamatrix_row_sum(
                self.ring.clone(),
                rows * cols,
                lats.iter().map(|lat| &lat.borrow().metamatrix).collect(),
            ),
        }
    }

    pub fn sum_pair<LatT: Borrow<LinearLattice<RS::Set>>>(
        &self,
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> LinearLattice<RS::Set> {
        self.sum(rows, cols, vec![lat1, lat2])
    }

    pub fn intersect<LatT: Borrow<LinearLattice<RS::Set>>>(
        &self,
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> LinearLattice<RS::Set> {
        if lats.len() == 0 {
            LinearLattice {
                rows,
                cols,
                metamatrix: MatrixStructure::new(self.ring.clone()).ident(rows * cols),
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

    pub fn intersect_pair(
        &self,
        rows: usize,
        cols: usize,
        lat1: impl Borrow<LinearLattice<RS::Set>>,
        lat2: impl Borrow<LinearLattice<RS::Set>>,
    ) -> LinearLattice<RS::Set> {
        assert_eq!(lat1.borrow().rows, rows);
        assert_eq!(lat1.borrow().cols, cols);
        assert_eq!(lat2.borrow().rows, rows);
        assert_eq!(lat2.borrow().cols, cols);
        let intersection_lattice = LinearLattice {
            rows,
            cols,
            metamatrix: metamatrix_row_intersection(
                self.ring.clone(),
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
        lat: &LinearLattice<RS::Set>,
    ) -> Vec<LinearLattice<RS::Set>> {
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

impl<RS: BezoutDomainStructure + ToStringStructure> LinearLatticeStructure<RS> {
    pub fn pprint(&self, lat: &LinearLattice<RS::Set>) {
        println!("Start Linear Lattice");
        for r in 0..lat.metamatrix.rows() {
            MatrixStructure::new(self.ring.clone()).pprint(&self.basis_matrix(lat, r));
        }
        println!("End Linear Lattice");
    }
}

impl<R: MetaType> MetaType for LinearLattice<R>
where
    R::Structure: BezoutDomainStructure,
{
    type Structure = LinearLatticeStructure<R::Structure>;

    fn structure() -> Rc<Self::Structure> {
        LinearLatticeStructure::new(R::structure()).into()
    }
}

impl<R: MetaType> LinearLattice<R>
where
    R::Structure: BezoutDomainStructure + ToStringStructure,
{
    pub fn pprint(&self) {
        Self::structure().pprint(self)
    }
}

impl<R: MetaType> PartialEq for LinearLattice<R>
where
    R::Structure: BezoutDomainStructure,
{
    fn eq(&self, other: &Self) -> bool {
        Self::structure().equal(self, other)
    }
}

impl<R: MetaType> LinearLattice<R>
where
    R::Structure: BezoutDomainStructure,
{
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        Self::structure().check_invariants(self)
    }

    pub fn from_span<MatT: Borrow<Matrix<R>>>(
        rows: usize,
        cols: usize,
        mats: Vec<MatT>,
    ) -> LinearLattice<R> {
        Self::structure().from_span(rows, cols, mats)
    }

    pub fn from_basis<MatT: Borrow<Matrix<R>>>(
        rows: usize,
        cols: usize,
        mats: Vec<MatT>,
    ) -> LinearLattice<R> {
        Self::structure().from_basis(rows, cols, mats)
    }

    pub fn rank(&self) -> usize {
        Self::structure().rank(self)
    }

    fn basis_row(&self, basis_num: usize) -> Matrix<R> {
        Self::structure().basis_row(self, basis_num)
    }

    pub fn basis_matrix(&self, r: usize) -> Matrix<R> {
        Self::structure().basis_matrix(self, r)
    }

    pub fn basis_matrices(&self) -> Vec<Matrix<R>> {
        Self::structure().basis_matrices(self)
    }

    pub fn basis_matrix_element<'b>(
        &self,
        lat: &'b LinearLattice<R>,
        basis_num: usize,
        r: usize,
        c: usize,
    ) -> &'b R {
        Self::structure().basis_matrix_element(lat, basis_num, r, c)
    }

    fn contains_row(&self, mat_as_row: impl Borrow<Matrix<R>>) -> bool {
        Self::structure().contains_row(self, mat_as_row)
    }

    pub fn contains_point<MatT: Borrow<Matrix<R>>>(&self, mat: MatT) -> bool {
        Self::structure().contains_point(self, mat)
    }

    //is lat a subset of self?
    pub fn contains_sublattice(&self, sublat: impl Borrow<LinearLattice<R>>) -> bool {
        Self::structure().contains_sublattice(self, sublat)
    }

    pub fn sum<LatT: Borrow<LinearLattice<R>>>(
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> LinearLattice<R> {
        Self::structure().sum(rows, cols, lats)
    }

    pub fn sum_pair<LatT: Borrow<LinearLattice<R>>>(
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> LinearLattice<R> {
        Self::structure().sum_pair(rows, cols, lat1, lat2)
    }

    pub fn intersect<LatT: Borrow<LinearLattice<R>>>(
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> LinearLattice<R> {
        Self::structure().intersect(rows, cols, lats)
    }

    pub fn intersect_pair(
        rows: usize,
        cols: usize,
        lat1: impl Borrow<LinearLattice<R>>,
        lat2: impl Borrow<LinearLattice<R>>,
    ) -> LinearLattice<R> {
        Self::structure().intersect_pair(rows, cols, lat1, lat2)
    }

    pub fn as_hyperplane_intersection(&self) -> Vec<LinearLattice<R>> {
        Self::structure().as_hyperplane_intersection(self)
    }
}

#[derive(Debug, Clone)]
pub enum AffineLatticeElements<Set: Clone> {
    Empty(),
    NonEmpty {
        offset: Matrix<Set>,        //offset.rows == 1 and offset.cols == self.cols
        linlat: LinearLattice<Set>, //linlat.rows == self.rows and linlat.cols == self.cols
    },
}

#[derive(Debug, Clone)]
pub struct AffineLattice<Set: Clone> {
    rows: usize,
    cols: usize,
    elems: AffineLatticeElements<Set>,
}

impl<Set: Clone> AffineLattice<Set> {
    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn elems(&self) -> &AffineLatticeElements<Set> {
        &self.elems
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AffineLatticeStructure<RS: BezoutDomainStructure> {
    ring: Rc<RS>,
}

impl<RS: BezoutDomainStructure> Structure for AffineLatticeStructure<RS> {
    type Set = AffineLattice<RS::Set>;
}

impl<RS: BezoutDomainStructure> AffineLatticeStructure<RS> {
    pub fn new(ring: Rc<RS>) -> Self {
        Self { ring }
    }
}

impl<RS: BezoutDomainStructure> AffineLatticeStructure<RS> {
    pub fn check_invariants(&self, lat: &AffineLattice<RS::Set>) -> Result<(), &'static str> {
        match &lat.elems {
            AffineLatticeElements::Empty() => {}
            AffineLatticeElements::NonEmpty { offset, linlat } => {
                match LinearLatticeStructure::new(self.ring.clone()).check_invariants(linlat) {
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

    pub fn rows(&self, lat: &AffineLattice<RS::Set>) -> usize {
        lat.rows
    }

    pub fn cols(&self, lat: &AffineLattice<RS::Set>) -> usize {
        lat.cols
    }

    pub fn rank(&self, lat: &AffineLattice<RS::Set>) -> Option<usize> {
        match &lat.elems {
            AffineLatticeElements::Empty() => None,
            AffineLatticeElements::NonEmpty {
                offset: _offset,

                linlat,
            } => Some(LinearLatticeStructure::new(self.ring.clone()).rank(linlat)),
        }
    }

    pub fn empty(&self, rows: usize, cols: usize) -> AffineLattice<RS::Set> {
        AffineLattice {
            rows,
            cols,
            elems: AffineLatticeElements::Empty(),
        }
    }

    pub fn to_offset_and_linear_lattice(
        &self,
        lat: AffineLattice<RS::Set>,
    ) -> Option<(Matrix<RS::Set>, LinearLattice<RS::Set>)> {
        match lat.elems {
            AffineLatticeElements::Empty() => None,
            AffineLatticeElements::NonEmpty { offset, linlat } => Some((offset, linlat)),
        }
    }

    pub fn take_point(&self, lat: AffineLattice<RS::Set>) -> Option<Matrix<RS::Set>> {
        match self.to_offset_and_linear_lattice(lat) {
            Some((offset, _linlat)) => Some(offset),
            None => None,
        }
    }

    pub fn from_offset_and_linear_lattice(
        &self,
        rows: usize,
        cols: usize,
        offset: Matrix<RS::Set>,
        linlat: LinearLattice<RS::Set>,
    ) -> AffineLattice<RS::Set> {
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

    pub fn contains_point(
        &self,
        lat: &AffineLattice<RS::Set>,
        mat: impl Borrow<Matrix<RS::Set>>,
    ) -> bool {
        let mat_struct = MatrixStructure::new(self.ring.clone());
        match &lat.elems {
            AffineLatticeElements::Empty() => false,
            AffineLatticeElements::NonEmpty { offset, linlat } => {
                LinearLatticeStructure::new(self.ring.clone()).contains_point(
                    linlat,
                    mat_struct
                        .add(&mat_struct.neg(offset.clone()), mat.borrow())
                        .unwrap(),
                )
            }
        }
    }

    //is other a subset of self?
    pub fn contains_sublattice(
        &self,
        lat: &AffineLattice<RS::Set>,
        other: impl Borrow<AffineLattice<RS::Set>>,
    ) -> bool {
        let mat_struct = MatrixStructure::new(self.ring.clone());
        let linlat_struct = LinearLatticeStructure::new(self.ring.clone());
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
                                .add(&linlat_struct.basis_matrix(other_linlat, bn), other_offset)
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

    pub fn equal(&self, lat1: &AffineLattice<RS::Set>, lat2: &AffineLattice<RS::Set>) -> bool {
        self.contains_sublattice(lat1, lat2) && self.contains_sublattice(lat2, lat1)
    }

    pub fn sum<LatT: Borrow<AffineLattice<RS::Set>>>(
        &self,
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> AffineLattice<RS::Set> {
        let mut sum_offset = MatrixStructure::new(self.ring.clone()).zero(rows, cols);
        let mut sum_linlats = vec![];
        for lat in &lats {
            assert_eq!(self.rows(lat.borrow()), rows);
            assert_eq!(self.cols(lat.borrow()), cols);
            match &lat.borrow().elems {
                AffineLatticeElements::Empty() => {
                    return self.empty(rows, cols);
                }
                AffineLatticeElements::NonEmpty { offset, linlat } => {
                    MatrixStructure::new(self.ring.clone())
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
                linlat: LinearLatticeStructure::new(self.ring.clone()).sum(rows, cols, sum_linlats),
            },
        }
    }

    pub fn sum_pair<LatT: Borrow<AffineLattice<RS::Set>>>(
        &self,
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> AffineLattice<RS::Set> {
        self.sum(rows, cols, vec![lat1, lat2])
    }

    pub fn intersect<LatT: Borrow<AffineLattice<RS::Set>>>(
        &self,
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> AffineLattice<RS::Set> {
        if lats.len() == 0 {
            AffineLattice {
                rows,
                cols,
                elems: AffineLatticeElements::NonEmpty {
                    offset: MatrixStructure::new(self.ring.clone()).zero(rows, cols),
                    linlat: LinearLattice {
                        rows,
                        cols,
                        metamatrix: MatrixStructure::new(self.ring.clone()).ident(rows * cols),
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
        offset: &Matrix<RS::Set>,
        linlat: &LinearLattice<RS::Set>,
    ) -> Matrix<RS::Set> {
        let mut metamat = MatrixStructure::new(self.ring.clone()).zero(
            1 + LinearLatticeStructure::new(self.ring.clone()).rank(linlat),
            1 + rows * cols,
        );
        *metamat.at_mut(0, 0).unwrap() = self.ring.one();
        for idx in 0..rows * cols {
            let (r, c) = idx_to_rc(rows, cols, idx);
            // println!("rows={} cols={} r={} c={} idx={}", rows, cols, r, c, idx);
            *metamat.at_mut(0, 1 + idx).unwrap() = offset.at(r, c).unwrap().clone();
        }
        for bn in 0..LinearLatticeStructure::new(self.ring.clone()).rank(linlat) {
            for idx in 0..rows * cols {
                let (r, c) = idx_to_rc(rows, cols, idx);
                *metamat.at_mut(0 + 1 + bn, 1 + idx).unwrap() =
                    LinearLatticeStructure::new(self.ring.clone())
                        .basis_matrix_element(linlat, bn, r, c)
                        .clone();
            }
        }
        metamat
    }

    pub fn intersect_pair(
        &self,
        rows: usize,
        cols: usize,
        lat1: impl Borrow<AffineLattice<RS::Set>>,
        lat2: impl Borrow<AffineLattice<RS::Set>>,
    ) -> AffineLattice<RS::Set> {
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
                    let int_metamat = metamatrix_row_intersection(
                        self.ring.clone(),
                        1 + rows * cols,
                        metamat1,
                        metamat2,
                    );

                    if int_metamat.rows() == 0 {
                        //the hyperlattice is just the origin, so the coresponding affine lattice - the intersection with the plane (1, *, ..., *) - is empty.
                        self.empty(rows, cols)
                    } else {
                        let (int_metamat_h, _u, _u_det, pivs) =
                            MatrixStructure::new(self.ring.clone())
                                .row_hermite_algorithm(int_metamat);
                        // MatrixStructure::new(self.ring).pprint(&int_metamat_h);
                        if self.ring.is_unit(int_metamat_h.at(0, 0).unwrap()) {
                            debug_assert!(self
                                .ring
                                .equal(int_metamat_h.at(0, 0).unwrap(), &self.ring.one()));
                        }
                        if self
                            .ring
                            .equal(int_metamat_h.at(0, 0).unwrap(), &self.ring.one())
                        {
                            let mut int_offset =
                                MatrixStructure::new(self.ring.clone()).zero(rows, cols);
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
                                        MatrixStructure::new(self.ring.clone()).zero(rows, cols);
                                    for idx in 0..rows * cols {
                                        let (r, c) = idx_to_rc(rows, cols, idx);
                                        *basis_mat.at_mut(r, c).unwrap() =
                                            int_metamat_h.at(1 + bn, 1 + idx).unwrap().clone();
                                    }
                                    basis_mat
                                })
                                .collect();
                            // MatrixStructure::new(self.ring).pprint(&int_offset);
                            // for basis_mat in &int_basis_mats {
                            // MatrixStructure::new(self.ring).pprint(basis_mat);
                            // }
                            self.from_offset_and_linear_lattice(
                                rows,
                                cols,
                                int_offset,
                                LinearLatticeStructure::new(self.ring.clone()).from_basis(
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

impl<RS: BezoutDomainStructure + ToStringStructure> AffineLatticeStructure<RS> {
    pub fn pprint(&self, lat: &AffineLattice<RS::Set>) {
        println!("Start Affine Lattice");
        match &lat.elems {
            AffineLatticeElements::Empty() => println!("Empty"),
            AffineLatticeElements::NonEmpty { offset, linlat } => {
                println!("Offset");
                MatrixStructure::new(self.ring.clone()).pprint(&offset);
                LinearLatticeStructure::new(self.ring.clone()).pprint(&linlat);
            }
        }
        println!("End Affine Lattice");
    }
}

impl<R: MetaType> MetaType for AffineLattice<R>
where
    R::Structure: BezoutDomainStructure,
{
    type Structure = AffineLatticeStructure<R::Structure>;

    fn structure() -> Rc<Self::Structure> {
        AffineLatticeStructure::new(R::structure()).into()
    }
}

impl<R: MetaType> AffineLattice<R>
where
    R::Structure: BezoutDomainStructure + ToStringStructure,
{
    pub fn pprint(&self) {
        Self::structure().pprint(self)
    }
}

impl<R: MetaType> PartialEq for AffineLattice<R>
where
    R::Structure: BezoutDomainStructure,
{
    fn eq(&self, other: &Self) -> bool {
        Self::structure().equal(self, other)
    }
}

impl<R: MetaType> AffineLattice<R>
where
    R::Structure: BezoutDomainStructure,
{
    pub fn check_invariants(&self) -> Result<(), &'static str> {
        Self::structure().check_invariants(self)
    }

    pub fn rank(&self) -> Option<usize> {
        Self::structure().rank(self)
    }

    pub fn empty(rows: usize, cols: usize) -> AffineLattice<R> {
        Self::structure().empty(rows, cols)
    }

    pub fn to_offset_and_linear_lattice(self) -> Option<(Matrix<R>, LinearLattice<R>)> {
        Self::structure().to_offset_and_linear_lattice(self)
    }

    pub fn from_offset_and_linear_lattice(
        rows: usize,
        cols: usize,
        offset: Matrix<R>,
        linlat: LinearLattice<R>,
    ) -> AffineLattice<R> {
        Self::structure().from_offset_and_linear_lattice(rows, cols, offset, linlat)
    }

    pub fn contains_point(&self, mat: impl Borrow<Matrix<R>>) -> bool {
        Self::structure().contains_point(self, mat)
    }

    //is other a subset of self?
    pub fn contains_sublattice(&self, other: impl Borrow<AffineLattice<R>>) -> bool {
        Self::structure().contains_sublattice(self, other)
    }

    pub fn sum<LatT: Borrow<AffineLattice<R>>>(
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> AffineLattice<R> {
        Self::structure().sum(rows, cols, lats)
    }

    pub fn sum_pair<LatT: Borrow<AffineLattice<R>>>(
        rows: usize,
        cols: usize,
        lat1: LatT,
        lat2: LatT,
    ) -> AffineLattice<R> {
        Self::structure().sum_pair(rows, cols, lat1, lat2)
    }

    pub fn intersect<LatT: Borrow<AffineLattice<R>>>(
        rows: usize,
        cols: usize,
        lats: Vec<LatT>,
    ) -> AffineLattice<R> {
        Self::structure().intersect(rows, cols, lats)
    }

    fn offset_linlat_to_metamat(
        rows: usize,
        cols: usize,
        offset: &Matrix<R>,
        linlat: &LinearLattice<R>,
    ) -> Matrix<R> {
        Self::structure().offset_linlat_to_metamat(rows, cols, offset, linlat)
    }

    pub fn intersect_pair(
        rows: usize,
        cols: usize,
        lat1: impl Borrow<Self>,
        lat2: impl Borrow<Self>,
    ) -> Self {
        Self::structure().intersect_pair(rows, cols, lat1, lat2)
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn linear_lattice_invariant() {
        let lattice = LinearLattice {
            metamatrix: Matrix::<Integer>::from_rows(vec![
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
            metamatrix: Matrix::<Integer>::from_rows(vec![
                vec![Integer::from(0), Integer::from(3), Integer::from(0)],
                vec![Integer::from(2), Integer::from(0), Integer::from(1)],
            ]),
            rows: 2,
            cols: 2,
        };
        assert!(lattice.check_invariants().is_err());

        let lattice = LinearLattice {
            metamatrix: Matrix::<Integer>::from_rows(vec![
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
        let lattice = LinearLattice::<Integer>::from_span(
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

        assert_ne!(
            LinearLattice::from_span(
                2,
                3,
                vec![
                    &Matrix::<Integer>::from_rows(vec![
                        vec![Integer::from(0), Integer::from(2), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                    &Matrix::from_rows(vec![
                        vec![Integer::from(0), Integer::from(4), Integer::from(0)],
                        vec![Integer::from(0), Integer::from(0), Integer::from(0)],
                    ]),
                ],
            ),
            LinearLattice::from_span(
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
        );
    }

    #[test]
    fn linear_lattice_sum_and_intersection() {
        {
            //standard basis sum and intersection
            let a = Matrix::<Integer>::from_rows(vec![
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

            assert_eq!(int, c.col_span());
            assert_eq!(sum, d.col_span());
        }

        {
            //sum and intersection as gcd and lcm
            let a = Matrix::<Integer>::from_rows(vec![
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

            assert_eq!(int, c.col_span());
            assert_eq!(sum, d.col_span());
        }

        {
            //triple intersection
            let a = Matrix::<Integer>::from_rows(vec![
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

            assert_eq!(
                int,
                Matrix::from_rows(vec![
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
            );
        }

        {
            //complex example
            let a = Matrix::<Integer>::from_rows(vec![
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

            assert_eq!(int, c.col_span());
            assert_eq!(sum, d.col_span());
        }
    }

    #[test]
    fn linear_lattice_extension_basis() {
        let a = Matrix::from_rows(vec![
            vec![Rational::from(1), Rational::from(0), Rational::from(0)],
            vec![Rational::from(1), Rational::from(0), Rational::from(0)],
            vec![Rational::from(-1), Rational::from(0), Rational::from(0)],
        ]);

        let b = Matrix::from_rows(vec![
            vec![Rational::from(1), Rational::from(1), Rational::from(0)],
            vec![Rational::from(1), Rational::from(1), Rational::from(0)],
            vec![Rational::from(1), Rational::from(-1), Rational::from(0)],
        ]);

        println!("a");
        a.pprint();
        println!("b");
        b.pprint();

        let ll_s = LinearLatticeStructure::new(Rational::structure());
        let ext = ll_s.extension_basis(&a.col_span(), &b.col_span());

        println!("ext");
        for v in ext {
            v.pprint();
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
                        Matrix::<Integer>::from_rows(vec![
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

    #[test]
    fn affine_lattice_sum_and_intersection() {
        let a1 = Matrix::<Integer>::from_rows(vec![
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

        assert_eq!(alat3, expected_alat3);
    }
}
