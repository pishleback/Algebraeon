use super::ring_of_integers::{
    RingOfIntegersWithIntegralBasisElement, RingOfIntegersWithIntegralBasisStructure,
};
use crate::{
    linear::{
        matrix::Matrix,
        subspace::{LinearLattice, LinearLatticeStructure},
    },
    structure::{IdealArithmeticStructure, IdealStructure},
};
use algebraeon_nzq::Integer;
use algebraeon_sets::structure::MetaType;

pub struct RingOfIntegersIdeal {
    // 1 column and n rows
    lattice: LinearLattice<Integer>,
}

impl RingOfIntegersWithIntegralBasisStructure {
    #[cfg(debug_assertions)]
    fn check_ideal(&self, ideal: &RingOfIntegersIdeal) {
        assert_eq!(ideal.lattice.rows(), self.degree());
        assert_eq!(ideal.lattice.cols(), 1);
        //TODO: check that it's actually an ideal
    }
}

impl RingOfIntegersIdeal {
    /// A basis of this ideal as a Z-module.
    pub fn integer_basis(&self) -> Vec<RingOfIntegersWithIntegralBasisElement> {
        self.lattice
            .basis_matrices()
            .into_iter()
            .map(|m| RingOfIntegersWithIntegralBasisElement::from_col(&m))
            .collect()
    }
}

impl IdealStructure for RingOfIntegersWithIntegralBasisStructure {
    type Ideal = RingOfIntegersIdeal;
}

impl IdealArithmeticStructure for RingOfIntegersWithIntegralBasisStructure {
    fn principal_ideal(&self, a: &Self::Set) -> Self::Ideal {
        todo!()
    }

    fn ideal_equal(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        LinearLatticeStructure::new(Integer::structure()).equal(&a.lattice, &b.lattice)
    }

    fn ideal_contains(&self, a: &Self::Ideal, b: &Self::Ideal) -> bool {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        LinearLatticeStructure::new(Integer::structure())
            .contains_sublattice(&a.lattice, &b.lattice)
    }

    fn ideal_intersection(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        todo!()
    }

    fn ideal_add(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        todo!()
    }

    fn ideal_mul(&self, a: &Self::Ideal, b: &Self::Ideal) -> Self::Ideal {
        #[cfg(debug_assertions)]
        {
            self.check_ideal(a);
            self.check_ideal(b);
        }
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn to_and_from_mat() {
    //     let m1 = Matrix::from_cols(vec![vec![1, 2], vec![3, 4]]);
    //     let i = RingOfIntegersIdeal::from_col_matrix(2, m1);
    //     let m2 = i.to_col_matrix();
    //     assert_eq!(m2.at(0, 0).unwrap(), &Integer::from(1));
    //     assert_eq!(m2.at(1, 0).unwrap(), &Integer::from(2));
    //     assert_eq!(m2.at(0, 1).unwrap(), &Integer::from(3));
    //     assert_eq!(m2.at(1, 1).unwrap(), &Integer::from(4));
    // }
}
