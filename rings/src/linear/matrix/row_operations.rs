use super::*;

#[derive(Debug)]
pub enum ElementaryOppType<RS: RingSignature> {
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

pub struct ElementaryOpp<RS: RingSignature> {
    ring: RS,
    transpose: bool, //false = row opp, true = column opp
    opp: ElementaryOppType<RS>,
}

impl<RS: IntegralDomainSignature> ElementaryOpp<RS> {
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
                let m = Matrix::construct(2, 2, |i, j| match (i, j) {
                    (0, 0) => a.clone(),
                    (1, 0) => b.clone(),
                    (0, 1) => c.clone(),
                    (1, 1) => d.clone(),
                    _ => unreachable!(),
                });
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

    pub fn new_row_opp(ring: RS, opp: ElementaryOppType<RS>) -> Self {
        Self {
            ring,
            transpose: false,
            opp,
        }
    }

    pub fn new_col_opp(ring: RS, opp: ElementaryOppType<RS>) -> Self {
        Self {
            ring,
            transpose: true,
            opp,
        }
    }

    pub fn det(&self) -> RS::Set {
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

    pub fn apply(&self, m: &mut Matrix<RS::Set>) {
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
