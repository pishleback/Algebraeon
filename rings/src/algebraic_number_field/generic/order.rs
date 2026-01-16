use crate::{
    algebraic_number_field::{
        AlgebraicIntegerRingSignature, AlgebraicNumberFieldSignature,
        FullRankIntegerSubmoduleWithBasis, FullRankIntegerSubmoduleWithBasisSignature, OrderIdeal,
    },
    matrix::{Matrix, SymmetricMatrix},
    module::{
        finitely_free_module::RingToFinitelyFreeModuleSignature,
        finitely_free_submodule::FinitelyFreeSubmodule,
    },
    structure::{
        AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
        CancellativeAdditionSignature, CancellativeMultiplicationSignature, CharZeroRingSignature,
        CharacteristicSignature, CommutativeMultiplicationSignature, DedekindDomainSignature,
        FinitelyFreeModuleSignature, IntegralDomainSignature,
        LeftDistributiveMultiplicationOverAddition, MultiplicationSignature,
        MultiplicativeAbsorptionMonoidSignature, MultiplicativeIntegralMonoidSignature,
        MultiplicativeMonoidSignature, OneSignature, RightDistributiveMultiplicationOverAddition,
        RingSignature, RingToIdealsSignature, RinglikeSpecializationSignature, SemiModuleSignature,
        SemiRingSignature, TryNegateSignature, TryReciprocalSignature, ZeroEqSignature,
        ZeroSignature,
    },
};
use algebraeon_nzq::{Integer, Natural};
use algebraeon_sets::structure::{
    BorrowedStructure, EqSignature, Function, InjectiveFunction, MetaType, SetSignature, Signature,
    ToStringSignature,
};

pub type RingOfIntegersWithIntegralBasis<K, KB> = OrderWithBasis<K, KB, true>;

#[derive(Debug, Clone)]
pub struct OrderWithBasis<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
> {
    full_rank_z_integer_submodule: FullRankIntegerSubmoduleWithBasis<K, KB>,
    products: SymmetricMatrix<Vec<Integer>>,
    one: Vec<Integer>,
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    OrderWithBasis<K, KB, MAXIMAL>
{
    fn check_is_order(&self) -> Result<(), String> {
        for b in self.full_rank_z_integer_submodule.basis() {
            if !self
                .full_rank_z_integer_submodule
                .anf()
                .is_algebraic_integer(b)
            {
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
        let full_rank_z_integer_submodule =
            FullRankIntegerSubmoduleWithBasis::new_unchecked(anf, basis);

        let products = products.map(|p| {
            full_rank_z_integer_submodule
                .outbound_order_to_anf_inclusion()
                .try_preimage(&p)
        });
        let products = products.unwrap_entries();
        if products.is_none() {
            return Err("An order must be closed under multiplication".to_string());
        }
        let products = products.unwrap();

        let one = full_rank_z_integer_submodule
            .outbound_order_to_anf_inclusion()
            .try_preimage(&full_rank_z_integer_submodule.anf().one());
        if one.is_none() {
            return Err("An order must contain 1".to_string());
        }
        let one = one.unwrap();
        let s = Self {
            full_rank_z_integer_submodule,
            products,
            one,
        };
        Ok(s)
    }

    /// The conductor of an order is the ideal of the order consisting of all elements such that multiplying by any element of the ring of integers produces an element of the order.
    pub fn conductor(&self) -> OrderIdeal {
        self.anf()
            .ring_of_integers()
            .suborder_conductor_as_ideal_of_order(self)
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
        let anf_disc = FullRankIntegerSubmoduleWithBasisSignature::<K>::anf(self).discriminant();
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
        self.full_rank_z_integer_submodule == other.full_rank_z_integer_submodule
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
        self.full_rank_z_integer_submodule.is_element(x)
    }
}

impl<
    K: AlgebraicNumberFieldSignature + ToStringSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
> ToStringSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn to_string(&self, elem: &Self::Set) -> String {
        self.full_rank_z_integer_submodule.to_string(elem)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> EqSignature
    for OrderWithBasis<K, KB, MAXIMAL>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.full_rank_z_integer_submodule.equal(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    RinglikeSpecializationSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn try_ring_restructure(&self) -> Option<impl EqSignature<Set = Self::Set> + RingSignature> {
        Some(self.clone())
    }

    fn try_char_zero_ring_restructure(
        &self,
    ) -> Option<impl EqSignature<Set = Self::Set> + CharZeroRingSignature> {
        Some(self.clone())
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> ZeroSignature
    for OrderWithBasis<K, KB, MAXIMAL>
{
    fn zero(&self) -> Self::Set {
        self.full_rank_z_integer_submodule.zero()
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    AdditionSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.full_rank_z_integer_submodule.add(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    CancellativeAdditionSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn try_sub(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.sub(a, b))
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    TryNegateSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn try_neg(&self, a: &Self::Set) -> Option<Self::Set> {
        Some(self.neg(a))
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    AdditiveMonoidSignature for OrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    AdditiveGroupSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.full_rank_z_integer_submodule.neg(a)
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.full_rank_z_integer_submodule.sub(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    FullRankIntegerSubmoduleWithBasisSignature<K> for OrderWithBasis<K, KB, MAXIMAL>
{
    fn anf(&self) -> &K {
        self.full_rank_z_integer_submodule.anf()
    }

    fn basis(&self) -> &Vec<<K>::Set> {
        self.full_rank_z_integer_submodule.basis()
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> OneSignature
    for OrderWithBasis<K, KB, MAXIMAL>
{
    fn one(&self) -> Self::Set {
        self.one.clone()
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    MultiplicationSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        let n = self.n();
        debug_assert!(self.is_element(a).is_ok());
        debug_assert!(self.is_element(b).is_ok());
        let mut t = self.zero();
        #[allow(clippy::needless_range_loop)]
        for i in 0..n {
            for j in 0..n {
                self.add_mut(
                    &mut t,
                    &self
                        .free_integer_submodule_restructure()
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
    CommutativeMultiplicationSignature for OrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    MultiplicativeMonoidSignature for OrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    MultiplicativeAbsorptionMonoidSignature for OrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    LeftDistributiveMultiplicationOverAddition for OrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    RightDistributiveMultiplicationOverAddition for OrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    SemiRingSignature for OrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    OrderWithBasis<K, KB, MAXIMAL>
{
    pub fn quotient_integer_submodule(
        &self,
        a: &FinitelyFreeSubmodule<Integer>,
        b: &FinitelyFreeSubmodule<Integer>,
    ) -> FinitelyFreeSubmodule<Integer> {
        #[cfg(debug_assertions)]
        {
            self.free_integer_submodule_restructure()
                .submodules()
                .is_element(a)
                .unwrap();
            self.free_integer_submodule_restructure()
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
        let module = self.free_integer_submodule_restructure();
        let integer_submodules = module.submodules();

        integer_submodules.intersect_list(
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
    TryReciprocalSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn try_reciprocal(&self, a: &Self::Set) -> Option<Self::Set> {
        if self.is_zero(a) {
            None
        } else {
            self.outbound_order_to_anf_inclusion().try_preimage(
                &self
                    .anf()
                    .try_reciprocal(&self.outbound_order_to_anf_inclusion().image(a))
                    .unwrap(),
            )
        }
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    CancellativeMultiplicationSignature for OrderWithBasis<K, KB, MAXIMAL>
{
    fn try_divide(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        Some(self.mul(a, &self.try_reciprocal(b)?))
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    MultiplicativeIntegralMonoidSignature for OrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    IntegralDomainSignature for OrderWithBasis<K, KB, MAXIMAL>
{
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
        self.full_rank_z_integer_submodule.anf()
    }

    fn to_anf(&self, x: &Self::Set) -> K::Set {
        self.outbound_order_to_anf_inclusion().image(x)
    }

    fn try_from_anf(&self, y: &K::Set) -> Option<Self::Set> {
        self.outbound_order_to_anf_inclusion().try_preimage(y)
    }

    fn integral_basis(&self) -> Vec<Self::Set> {
        self.free_integer_submodule_restructure().basis_vecs()
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> OrderWithBasis<K, KB, true> {
    pub fn suborder_conductor_as_ideal_of_order<KOB: BorrowedStructure<K>, const MAXIMAL: bool>(
        &self,
        order: &OrderWithBasis<K, KOB, MAXIMAL>,
    ) -> OrderIdeal {
        let anf = order.anf();
        debug_assert_eq!(FullRankIntegerSubmoduleWithBasisSignature::anf(self), anf);
        debug_assert_eq!(AlgebraicIntegerRingSignature::anf(self), anf);

        let order_to_maximal_order = self
            .inbound_integer_submodule_inclusion_unchecked::<OrderWithBasis<K, KOB, MAXIMAL>, _>(
                order,
            );

        let ideal = order
            .ideals()
            .outbound_integer_submodules_inclusion()
            .try_preimage(
                &order_to_maximal_order
                    .integer_submodules_inclusion()
                    .try_preimage(
                        &self.quotient_integer_submodule(
                            &order_to_maximal_order.integer_submodules_inclusion().image(
                                &Integer::structure()
                                    .free_module(anf.n())
                                    .submodules()
                                    .full_submodule(),
                            ),
                            &Integer::structure()
                                .free_module(anf.n())
                                .submodules()
                                .full_submodule(),
                        ),
                    )
                    .unwrap(),
            )
            .unwrap();

        #[cfg(debug_assertions)]
        order.ideals().is_element(&ideal).unwrap();

        ideal
    }

    pub fn suborder_conductor_as_ideal_of_self<KOB: BorrowedStructure<K>, const MAXIMAL: bool>(
        &self,
        order: &OrderWithBasis<K, KOB, MAXIMAL>,
    ) -> OrderIdeal {
        let anf = order.anf();
        debug_assert_eq!(FullRankIntegerSubmoduleWithBasisSignature::anf(self), anf);
        debug_assert_eq!(AlgebraicIntegerRingSignature::anf(self), anf);

        let order_to_maximal_order = self
            .inbound_integer_submodule_inclusion_unchecked::<OrderWithBasis<K, KOB, MAXIMAL>, _>(
                order,
            );

        let ideal = self
            .ideals()
            .outbound_integer_submodules_inclusion()
            .try_preimage(
                &self.quotient_integer_submodule(
                    &order_to_maximal_order.integer_submodules_inclusion().image(
                        &Integer::structure()
                            .free_module(anf.n())
                            .submodules()
                            .full_submodule(),
                    ),
                    &Integer::structure()
                        .free_module(anf.n())
                        .submodules()
                        .full_submodule(),
                ),
            )
            .unwrap();

        #[cfg(debug_assertions)]
        self.ideals().is_element(&ideal).unwrap();

        ideal
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
                assert_eq!(
                    roi.try_reciprocal(&roi.neg(&roi.one())).unwrap(),
                    roi.neg(&roi.one())
                );
                assert_eq!(roi.try_reciprocal(&roi.one()).unwrap(), roi.one());
                assert!(roi.try_reciprocal(&alpha).is_none());
                assert!(roi.try_reciprocal(&beta).is_none());
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

    #[test]
    fn order_conductor() {
        // Q[sqrt(-3)]
        let anf = parse_rational_polynomial("x^2+3", "x")
            .unwrap()
            .algebraic_number_field()
            .unwrap();

        // Z[1/2 + 1/2 sqrt(-3)]
        let roi = anf.ring_of_integers();

        // Z[sqrt(-3)]
        let order = anf
            .order(vec![
                parse_rational_polynomial("1", "x").unwrap(),
                parse_rational_polynomial("x", "x").unwrap(),
            ])
            .unwrap();

        assert!(order.ideals().equal(
            &order.conductor(),
            &order.ideals().generated_ideal(vec![
                vec![Integer::from(1), Integer::from(1)],
                vec![Integer::from(0), Integer::from(2)]
            ])
        ));

        assert!(
            roi.ideals().equal(
                &roi.suborder_conductor_as_ideal_of_self(&order),
                &roi.ideals()
                    .principal_ideal(&roi.from_int(Integer::from(2)))
            )
        );
    }
}
