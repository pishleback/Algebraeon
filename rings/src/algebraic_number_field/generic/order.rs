use crate::{
    algebraic_number_field::{
        AlgebraicNumberFieldFullRankZSubmoduleWithBasis, AlgebraicNumberFieldSignature,
    },
    matrix::SymmetricMatrix,
    module::finitely_free_module::FinitelyFreeModuleStructure,
    structure::{
        AdditiveGroupSignature, AdditiveMonoidEqSignature, AdditiveMonoidSignature,
        CharZeroRingSignature, CharacteristicSignature, DedekindDomainSignature,
        FiniteDimensionalFieldExtension, FiniteRankFreeRingExtension, IntegralDomainSignature,
        MetaCharZeroRing, RingDivisionError, RingSignature, SemiModuleSignature, SemiRingSignature,
        SemiRingUnitsSignature,
    },
};
use algebraeon_nzq::{Integer, IntegerCanonicalStructure, Natural};
use algebraeon_sets::structure::{
    BorrowedStructure, EqSignature, Function, InjectiveFunction, Morphism, SetSignature, Signature,
};
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct AlgebraicNumberFieldOrderWithBasis<
    K: AlgebraicNumberFieldSignature,
    KB: BorrowedStructure<K>,
    const MAXIMAL: bool,
> {
    full_rank_z_submodule: AlgebraicNumberFieldFullRankZSubmoduleWithBasis<K, KB>,
    products: SymmetricMatrix<Vec<Integer>>,
    one: Vec<Integer>,
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    fn check_is_order(&self) -> Result<(), String> {
        for b in self.full_rank_z_submodule.basis() {
            if !self.full_rank_z_submodule.anf().is_algebraic_integer(b) {
                return Err("Basis vectors must be algebraic integers for an order".to_string());
            }
        }
        Ok(())
    }

    fn new_impl(anf: KB, basis: Vec<K::Set>) -> Result<Self, String> {
        fn make_products<K: AlgebraicNumberFieldSignature>(
            anf: &K,
            basis: &Vec<K::Set>,
        ) -> Result<SymmetricMatrix<Vec<Integer>>, String> {
            let n = anf.n();
            let mut s = SymmetricMatrix::filled(n, vec![]);
            for r in 0..n {
                for c in r..n {
                    let mut x = vec![];
                    for y in anf
                        .inbound_finite_dimensional_rational_extension()
                        .to_vec(&anf.mul(&basis[r], &basis[c]))
                    {
                        if let Some(y) = y.try_to_int() {
                            x.push(y)
                        } else {
                            return Err(format!("{} is not an integer", y));
                        }
                    }
                    s.set(r, c, x).unwrap();
                }
            }
            Ok(s)
        }
        let products = make_products(anf.borrow(), &basis)
            .map_err(|err| format!("Not an order: Not closed under multiplication: {}", err))
            .unwrap();
        let anf_one = anf.borrow().one();
        let abelian_group =
            AlgebraicNumberFieldFullRankZSubmoduleWithBasis::new_unchecked(anf, basis);
        let one = abelian_group
            .outbound_anf_inclusion()
            .try_preimage(&anf_one);
        if one.is_none() {
            return Err("An order must contain 1".to_string());
        }
        let one = one.unwrap();
        let s = Self {
            full_rank_z_submodule: abelian_group,
            products,
            one,
        };
        Ok(s)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>>
    AlgebraicNumberFieldOrderWithBasis<K, KB, false>
{
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

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>>
    AlgebraicNumberFieldOrderWithBasis<K, KB, true>
{
    fn check_is_maximal(&self) -> Result<(), String> {
        let self_disc = self.discriminant();
        let anf_disc = self.anf().discriminant();
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

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    pub fn anf(&self) -> &K {
        self.full_rank_z_submodule.anf()
    }

    pub fn basis(&self) -> &Vec<K::Set> {
        self.full_rank_z_submodule.basis()
    }

    pub fn n(&self) -> usize {
        self.full_rank_z_submodule.n()
    }

    pub fn full_rank_z_submodule_restructure(
        &self,
    ) -> &AlgebraicNumberFieldFullRankZSubmoduleWithBasis<K, KB> {
        &self.full_rank_z_submodule
    }

    pub fn free_z_module_restructure(
        &self,
    ) -> FinitelyFreeModuleStructure<IntegerCanonicalStructure, IntegerCanonicalStructure> {
        self.full_rank_z_submodule.free_z_module_restructure()
    }

    pub fn discriminant(&self) -> Integer {
        self.anf()
            .inbound_finite_dimensional_rational_extension()
            .discriminant(self.full_rank_z_submodule.basis())
            .try_to_int()
            .unwrap()
    }

    pub fn into_outbound_anf_inclusion(
        self,
    ) -> anf_inclusion::AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, Self> {
        anf_inclusion::AlgebraicNumberFieldOrderWithBasisInclusion::new(self)
    }

    pub fn outbound_anf_inclusion<'a>(
        &'a self,
    ) -> anf_inclusion::AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, &'a Self> {
        anf_inclusion::AlgebraicNumberFieldOrderWithBasisInclusion::new(self)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> PartialEq
    for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    fn eq(&self, other: &Self) -> bool {
        self.full_rank_z_submodule == other.full_rank_z_submodule
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> Eq
    for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> Signature
    for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> SetSignature
    for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    type Set = Vec<Integer>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        self.full_rank_z_submodule.is_element(x)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> EqSignature
    for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        self.full_rank_z_submodule_restructure().equal(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    AdditiveMonoidSignature for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    fn zero(&self) -> Self::Set {
        self.full_rank_z_submodule_restructure().zero()
    }

    fn add(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.full_rank_z_submodule_restructure().add(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    AdditiveGroupSignature for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    fn neg(&self, a: &Self::Set) -> Self::Set {
        self.full_rank_z_submodule_restructure().neg(a)
    }

    fn sub(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        self.full_rank_z_submodule_restructure().sub(a, b)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    SemiRingSignature for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
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
                        .free_z_module_restructure()
                        .scalar_mul(self.products.get(i, j).unwrap(), &(&a[i] * &b[j])),
                );
            }
        }
        debug_assert!(
            self.equal(
                &self
                    .outbound_anf_inclusion()
                    .try_preimage(&self.anf().mul(
                        &self.outbound_anf_inclusion().image(a),
                        &self.outbound_anf_inclusion().image(b)
                    ))
                    .unwrap(),
                &t
            )
        );
        t
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool> RingSignature
    for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    fn is_reduced(&self) -> Result<bool, String> {
        Ok(true)
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    CharacteristicSignature for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    fn characteristic(&self) -> Natural {
        Natural::ZERO
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    SemiRingUnitsSignature for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    fn inv(&self, a: &Self::Set) -> Result<Self::Set, crate::structure::RingDivisionError> {
        if self.is_zero(a) {
            Err(RingDivisionError::DivideByZero)
        } else if let Some(a_inv) = self.outbound_anf_inclusion().try_preimage(
            &self
                .anf()
                .inv(&self.outbound_anf_inclusion().image(a))
                .unwrap(),
        ) {
            Ok(a_inv)
        } else {
            Err(RingDivisionError::NotDivisible)
        }
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>, const MAXIMAL: bool>
    IntegralDomainSignature for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
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
    CharZeroRingSignature for AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>
{
    fn try_to_int(&self, x: &Self::Set) -> Option<Integer> {
        self.anf()
            .try_to_int(&self.outbound_anf_inclusion().image(x))
    }
}

impl<K: AlgebraicNumberFieldSignature, KB: BorrowedStructure<K>> DedekindDomainSignature
    for AlgebraicNumberFieldOrderWithBasis<K, KB, true>
{
}

mod anf_inclusion {
    use super::*;

    #[derive(Debug, Clone)]
    pub struct AlgebraicNumberFieldOrderWithBasisInclusion<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>>,
    > {
        _k: PhantomData<K>,
        _kb: PhantomData<KB>,
        order: OB,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>>,
    > AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        pub fn new(order: OB) -> Self {
            Self {
                _k: PhantomData,
                _kb: PhantomData,
                order,
            }
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>>,
    > Morphism<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>, K>
        for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        fn domain(&self) -> &AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL> {
            self.order.borrow()
        }

        fn range(&self) -> &K {
            self.order.borrow().anf()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>>,
    > Function<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>, K>
        for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        fn image(&self, x: &Vec<Integer>) -> <K as SetSignature>::Set {
            self.order
                .borrow()
                .full_rank_z_submodule_restructure()
                .outbound_anf_inclusion()
                .image(x)
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        KB: BorrowedStructure<K>,
        const MAXIMAL: bool,
        OB: BorrowedStructure<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>>,
    > InjectiveFunction<AlgebraicNumberFieldOrderWithBasis<K, KB, MAXIMAL>, K>
        for AlgebraicNumberFieldOrderWithBasisInclusion<K, KB, MAXIMAL, OB>
    {
        fn try_preimage(&self, y: &<K as SetSignature>::Set) -> Option<Vec<Integer>> {
            self.order
                .borrow()
                .full_rank_z_submodule_restructure()
                .outbound_anf_inclusion()
                .try_preimage(y)
        }
    }
}
