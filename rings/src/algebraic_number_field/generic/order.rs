use crate::{
    algebraic_number_field::{
        AlgebraicIntegerRingSignature, AlgebraicNumberFieldSignature, FullRankSublatticeWithBasis,
        FullRankSublatticeWithBasisSignature,
    },
    matrix::{Matrix, SymmetricMatrix},
    module::finitely_free_submodule::FinitelyFreeSubmodule,
    structure::{
        AdditiveGroupSignature, AdditiveMonoidEqSignature, AdditiveMonoidSignature,
        CharZeroRingSignature, CharacteristicSignature, DedekindDomainSignature,
        FinitelyFreeModuleSignature, IntegralDomainSignature, RingDivisionError, RingSignature,
        SemiModuleSignature, SemiRingSignature, SemiRingUnitsSignature,
    },
};
use algebraeon_nzq::{Integer, Natural};
use algebraeon_sets::structure::{
    BorrowedStructure, EqSignature, Function, InjectiveFunction, Morphism, SetSignature, Signature,
    ToStringSignature,
};
use std::marker::PhantomData;

pub type RingOfIntegersWithIntegralBasis<K, KB> = OrderWithBasis<K, KB, true>;

#[derive(Debug, Clone)]
pub struct OrderWithBasis<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
> {
    full_rank_z_sublattice: FullRankSublatticeWithBasis<K, KB>,
    products: SymmetricMatrix<Vec<Integer>>,
    one: Vec<Integer>,
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    OrderWithBasis<K, KB, MAXIMAL>
{
    fn check_is_order(&self) -> Result<(), String> {
        for b in self.full_rank_z_sublattice.basis() {
            if !self.full_rank_z_sublattice.anf().is_algebraic_integer(b) {
                return Err("Basis vectors must be algebraic integers for an order".to_string());
            }
        }
        Ok(())
    }

    fn new_impl(anf: KB, basis: Vec<K::Set>) -> Result<Self, String> {
        let n = anf.borrow().n();
        let products = SymmetricMatrix::construct_bottom_left(n, |r, c| {
            anf.borrow().mul(&basis[r], &basis[c])
        });
        let full_rank_z_sublattice = FullRankSublatticeWithBasis::new_unchecked(anf, basis);

        let products = products.map(|p| {
            full_rank_z_sublattice
                .outbound_order_to_anf_inclusion()
                .try_preimage(&p)
        });
        let products = products.unwrap_entries();
        if products.is_none() {
            return Err("An order must be closed under multiplication".to_string());
        }
        let products = products.unwrap();

        let one = full_rank_z_sublattice
            .outbound_order_to_anf_inclusion()
            .try_preimage(&full_rank_z_sublattice.anf().one());
        if one.is_none() {
            return Err("An order must contain 1".to_string());
        }
        let one = one.unwrap();
        let s = Self {
            full_rank_z_sublattice,
            products,
            one,
        };
        Ok(s)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> OrderWithBasis<K, KB, false> {
    pub fn new(anf: KB, basis: Vec<K::Set>) -> Result<Self, String> {
        let s = Self::new_impl(anf, basis)?;
        s.check_is_order()?;
        Ok(s)
    }

    pub fn new_unchecked(anf: KB, basis: Vec<K::Set>) -> Self {
        let s = Self::new_impl(anf, basis).unwrap();
        #[cfg(debug_assertions)]
        s.check_is_order().unwrap();
        s
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> OrderWithBasis<K, KB, true> {
    fn check_is_maximal(&self) -> Result<(), String> {
        let self_disc = self.discriminant();
        let anf_disc = FullRankSublatticeWithBasisSignature::<K>::anf(self).discriminant();
        if self_disc != anf_disc {
            debug_assert!(self_disc > anf_disc);
            return Err("Not maximal".to_string());
        }
        Ok(())
    }

    pub fn new_maximal(anf: KB, basis: Vec<K::Set>) -> Result<Self, String> {
        let s = Self::new_impl(anf, basis)?;
        s.check_is_order()?;
        s.check_is_maximal()?;
        Ok(s)
    }

    pub fn new_maximal_unchecked(anf: KB, basis: Vec<K::Set>) -> Self {
        let s = Self::new_impl(anf, basis).unwrap();
        #[cfg(debug_assertions)]
        s.check_is_order().unwrap();
        s.check_is_maximal().unwrap();
        s
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> PartialEq
    for OrderWithBasis<K, KB, MAXIMAL>
{
    fn eq(&self, other: &Self) -> bool {
        self.full_rank_z_sublattice == other.full_rank_z_sublattice
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> Eq
    for OrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> Signature
    for OrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> SetSignature
    for OrderWithBasis<K, KB, MAXIMAL>
{
    type Set = Vec<Integer>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        self.full_rank_z_sublattice.is_element(x)
    }
}

impl<
    K: AlgebraicNumberFieldSignature + ToStringSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
> ToStringSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn to_string(&self, elem: &Self::Set) -> String {
        self.full_rank_z_sublattice.to_string(elem)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> EqSignature
    for OrderWithBasis<K, KB, MAXIMAL>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.full_rank_z_sublattice.equal(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    AdditiveMonoidSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn zero(&self) -> Self::Set {
        self.full_rank_z_sublattice.zero()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.full_rank_z_sublattice.add(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    AdditiveGroupSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.full_rank_z_sublattice.neg(a)
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.full_rank_z_sublattice.sub(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    FullRankSublatticeWithBasisSignature<K> for OrderWithBasis<K, KB, MAXIMAL>
{
    fn anf(&self) -> &K {
        self.full_rank_z_sublattice.anf()
    }

    fn basis(&self) -> &Vec<<K>::Set> {
        self.full_rank_z_sublattice.basis()
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    SemiRingSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn one(&self) -> Self::Set {
        self.one.clone()
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let n = self.n();
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        let mut t = self.zero();
        for i in 0..n {
            for j in 0..n {
                self.add_mut(
                    &mut t,
                    &self
                        .free_lattice_restructure()
                        .scalar_mul(self.products.get(i, j).unwrap(), &(&a[i] * &b[j])),
                );
            }
        }
        debug_assert!(
            self.equal(
                &self
                    .outbound_order_to_anf_inclusion()
                    .try_preimage(&self.anf().mul(
                        &self.outbound_order_to_anf_inclusion().image(a),
                        &self.outbound_order_to_anf_inclusion().image(b)
                    ))
                    .unwrap(),
                &t
            )
        );
        t
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    OrderWithBasis<K, KB, MAXIMAL>
{
    pub fn quotient_sublattice(
        &self,
        a: &FinitelyFreeSubmodule<Integer>,
        b: &FinitelyFreeSubmodule<Integer>,
    ) -> FinitelyFreeSubmodule<Integer> {
        #[cfg(debug_assertions)]
        {
            self.free_lattice_restructure()
                .submodules()
                .is_element(a)
                .unwrap();
            self.free_lattice_restructure()
                .submodules()
                .is_element(b)
                .unwrap();
        }

        /*
        We want the ideal of all points x in the ring R such that xj belongs to I for all j in J.
        It's sufficient to check on a basis of points j in J.
        For each j in a basis for J we find the space of points x such that xj belongs to I and take their intersection.
         */

        let n = self.n();
        let module = self.free_lattice_restructure();
        let sublattices = module.submodules();

        sublattices.intersect_list(
            b.basis()
                .into_iter()
                .map(|bv| {
                    // column matrix for multiplication by jb wrt the integer basis for the ring
                    let bv_mulmat = Matrix::join_cols(
                        n,
                        (0..n)
                            .map(|c| {
                                Matrix::<Integer>::from_col(self.mul(&bv, &module.basis_element(c)))
                            })
                            .collect(),
                    );
                    bv_mulmat.col_preimage(a)
                })
                .collect(),
        )
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> RingSignature
    for OrderWithBasis<K, KB, MAXIMAL>
{
    fn is_reduced(&self) -> Result<bool, String> {
        Ok(true)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    CharacteristicSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    SemiRingUnitsSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, crate::structure::RingDivisionError> {
        if self.is_zero(a) {
            Err(RingDivisionError::DivideByZero)
        } else if let Some(a_inv) = self.outbound_order_to_anf_inclusion().try_preimage(
            &self
                .anf()
                .inv(&self.outbound_order_to_anf_inclusion().image(a))
                .unwrap(),
        ) {
            Ok(a_inv)
        } else {
            Err(RingDivisionError::NotDivisible)
        }
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    IntegralDomainSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn div(
        &self,
        a: &Self::Set,
        b: &Self::Set,
    ) -> Result<Self::Set, crate::structure::RingDivisionError> {
        match self.inv(b) {
            Ok(b_inv) => Ok(self.mul(a, &b_inv)),
            Err(err) => Err(err),
        }
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    CharZeroRingSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer> {
        self.anf()
            .try_to_int(&self.outbound_order_to_anf_inclusion().image(x))
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> DedekindDomainSignature
    for OrderWithBasis<K, KB, true>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> AlgebraicIntegerRingSignature<K>
    for OrderWithBasis<K, KB, true>
{
    fn anf(&self) -> &K {
        self.full_rank_z_sublattice.anf()
    }

    fn to_anf(&self, x: &Self::Set) -> K::Set {
        self.outbound_order_to_anf_inclusion().image(x)
    }

    fn try_from_anf(&self, y: &K::Set) -> Option<Self::Set> {
        self.outbound_order_to_anf_inclusion().try_preimage(y)
    }

    fn integral_basis(&self) -> Vec<Self::Set> {
        self.free_lattice_restructure().basis_vecs()
    }
}

mod anf_inclusion {
    use super::*;
    use crate::structure::RingHomomorphism;

    #[derive(Debug, Clone)]
    pub struct AlgebraicNumberFieldOrderWithBasisInclusion<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
    > {
        _k: PhantomData<K>,
        _kb: PhantomData<KB>,
        order: OB,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
    > AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        pub fn new(order: OB) -> Self {
            Self {
                _k: PhantomData,
                _kb: PhantomData,
                order,
            }
        }

        pub fn order(&self) -> &OrderWithBasis<K, KB, MAXIMAL> {
            self.order.borrow()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
    > Morphism<OrderWithBasis<K, KB, MAXIMAL>, K>
        for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        fn domain(&self) -> &OrderWithBasis<K, KB, MAXIMAL> {
            self.order()
        }

        fn range(&self) -> &K {
            self.order().anf()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
    > Function<OrderWithBasis<K, KB, MAXIMAL>, K>
        for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        fn image(&self, x: &Vec<Integer>) -> <K as SetSignature>::Set {
            self.order()
                .full_rank_z_sublattice
                .outbound_order_to_anf_inclusion()
                .image(x)
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
    > InjectiveFunction<OrderWithBasis<K, KB, MAXIMAL>, K>
        for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        fn try_preimage(&self, y: &<K as SetSignature>::Set) -> Option<Vec<Integer>> {
            self.order()
                .full_rank_z_sublattice
                .outbound_order_to_anf_inclusion()
                .try_preimage(y)
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<OrderWithBasis<K, KB, MAXIMAL>>,
    > RingHomomorphism<OrderWithBasis<K, KB, MAXIMAL>, K>
        for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        parsing::parse_rational_polynomial,
        polynomial::Polynomial,
        structure::{IdealsArithmeticSignature, IntoErgonomic, RingToIdealsSignature},
    };
    use algebraeon_nzq::Rational;

    #[test]
    fn ring_of_integer_arithmetic() {
        let x = Polynomial::<Rational>::var().into_ergonomic();

        // Take the integral basis (0 + x, 1/2 + 1/2x)
        let a = Polynomial::<Rational>::from_coeffs(vec![Rational::ZERO, Rational::ONE]);
        let b = Polynomial::<Rational>::from_coeffs(vec![Rational::ONE_HALF, Rational::ONE_HALF]);

        let anf = (x.pow(2) + 7)
            .into_verbose()
            .algebraic_number_field()
            .unwrap();
        let roi =
            RingOfIntegersWithIntegralBasis::new_maximal(anf, vec![a.clone(), b.clone()]).unwrap();

        {
            assert_eq!(
                roi.outbound_order_to_anf_inclusion()
                    .image(&vec![Integer::from(1), Integer::from(4)]),
                (2 + 3 * &x).into_verbose()
            );
        }

        {
            assert!(
                roi.outbound_order_to_anf_inclusion()
                    .try_preimage(&Polynomial::<Rational>::from_coeffs(vec![
                        Rational::ONE_HALF,
                        Rational::ONE,
                    ]))
                    .is_none()
            );

            let c = roi
                .outbound_order_to_anf_inclusion()
                .try_preimage(&Polynomial::<Rational>::from_coeffs(vec![
                    Rational::from(2),
                    Rational::from(3),
                ]))
                .unwrap();
            assert_eq!(c.len(), 2);
            assert_eq!(c[0], Integer::from(1));
            assert_eq!(c[1], Integer::from(4));
        }

        {
            // 0 = 0 * (0+x) + 0 * (1/2 + 1/2x)
            let zero = roi.zero();
            assert_eq!(zero.len(), 2);
            assert_eq!(zero[0], Integer::ZERO);
            assert_eq!(zero[1], Integer::ZERO);
        }

        {
            // 1 = -1 * (0+x) + 2 * (1/2 + 1/2x)
            let one = roi.one();
            assert_eq!(one.len(), 2);
            assert_eq!(one[0], Integer::from(-1));
            assert_eq!(one[1], Integer::from(2));
        }

        {
            let alpha = roi
                .outbound_order_to_anf_inclusion()
                .try_preimage(&(2 + 3 * &x).into_verbose())
                .unwrap();
            let beta = roi
                .outbound_order_to_anf_inclusion()
                .try_preimage(&(-1 + 2 * &x).into_verbose())
                .unwrap();

            {
                let gamma = roi
                    .outbound_order_to_anf_inclusion()
                    .try_preimage(&(1 + 5 * &x).into_verbose())
                    .unwrap();
                // (2 + 3x) + (-1 + 2x) = 1 + 5x
                assert_eq!(roi.add(&alpha, &beta), gamma);
            }

            {
                let gamma = roi
                    .outbound_order_to_anf_inclusion()
                    .try_preimage(&(-44 + &x).into_verbose())
                    .unwrap();
                // x^2 = -7 so
                // (2 + 3x) * (-1 + 2x) = -44 + x
                assert_eq!(roi.mul(&alpha, &beta), gamma);
            }

            {
                let gamma = roi
                    .outbound_order_to_anf_inclusion()
                    .try_preimage(&(-2 - 3 * &x).into_verbose())
                    .unwrap();
                // -(2 + 3x) = -2 - 3x
                assert_eq!(roi.neg(&alpha), gamma);
            }

            {
                assert_eq!(roi.inv(&roi.neg(&roi.one())).unwrap(), roi.neg(&roi.one()));
                assert_eq!(roi.inv(&roi.one()).unwrap(), roi.one());
                assert!(roi.inv(&alpha).is_err());
                assert!(roi.inv(&beta).is_err());
            }
        }

        println!("{:?}", roi);
    }

    #[test]
    fn order_ideals() {
        // Q[sqrt(-3)]
        let anf = parse_rational_polynomial("x^2+3", "x")
            .unwrap()
            .algebraic_number_field()
            .unwrap();

        // Z[sqrt(-3)]
        let order = anf
            .order(vec![
                parse_rational_polynomial("1", "x").unwrap(),
                parse_rational_polynomial("x", "x").unwrap(),
            ])
            .unwrap();

        let ideal6 = order
            .ideals()
            .principal_ideal(&order.from_int(Integer::from(6)));

        let ideal15 = order
            .ideals()
            .principal_ideal(&order.from_int(Integer::from(15)));

        assert!(
            order.ideals().equal(
                &order.ideals().add(&ideal6, &ideal15),
                &order
                    .ideals()
                    .principal_ideal(&order.from_int(Integer::from(3)))
            )
        );

        assert!(
            order.ideals().equal(
                &order.ideals().intersect(&ideal6, &ideal15),
                &order
                    .ideals()
                    .principal_ideal(&order.from_int(Integer::from(30)))
            )
        );

        assert!(
            order.ideals().equal(
                &order.ideals().quotient(&ideal15, &ideal6),
                &order
                    .ideals()
                    .principal_ideal(&order.from_int(Integer::from(5)))
            )
        );
    }
}
