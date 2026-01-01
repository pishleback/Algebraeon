use std::marker::PhantomData;

use crate::{
    algebraic_number_field::{OrderIdeal, OrderWithBasis, RingOfIntegersExtension},
    integer::ideal::IntegerIdealsStructure,
    module::finitely_free_module::{
        FinitelyFreeModuleStructure, RingToFinitelyFreeModuleSignature,
    },
    structure::{
        AdditiveGroupSignature, CharZeroFieldSignature, CharZeroRingSignature,
        DedekindDomainIdealsSignature, DedekindDomainSignature, FiniteDimensionalFieldExtension,
        FiniteRankFreeRingExtension, MetaGreatestCommonDivisor, RingHomomorphism,
        RingToIdealsSignature,
    },
};
use algebraeon_nzq::{
    Integer, IntegerCanonicalStructure, Rational, RationalCanonicalStructure, traits::Fraction,
};
use algebraeon_sets::structure::{
    BorrowedStructure, FiniteSetSignature, Function, InjectiveFunction, MetaType, Morphism,
    SetSignature,
};

/// An algebraic number field is a field of characteristic zero such that
/// the inclusion of its rational subfield is finite dimensional
pub trait AlgebraicNumberFieldSignature: CharZeroFieldSignature {
    type Basis: FiniteSetSignature;
    type RationalInclusion<B: BorrowedStructure<Self>>: FiniteDimensionalFieldExtension<RationalCanonicalStructure, Self>;

    fn inbound_finite_dimensional_rational_extension<'a>(
        &'a self,
    ) -> Self::RationalInclusion<&'a Self>;
    fn into_inbound_finite_dimensional_rational_extension(self) -> Self::RationalInclusion<Self>;

    /// The dimension of this algebraic number field as a vector space over the rational numbers
    fn n(&self) -> usize {
        self.inbound_finite_dimensional_rational_extension()
            .degree()
    }

    /// An element which generates this algebraic number field when adjoined to the rational numbers
    /// Such an element always exists by the primitive element theorem
    fn generator(&self) -> Self::Set;

    /// Determine whether an element is integral over the integers i.e. is it a root of a monic integer polynomial
    fn is_algebraic_integer(&self, a: &Self::Set) -> bool;

    /// The discriminant of this algebraic number field i.e. the discriminant of its ring of integers
    /// Implementations should not compute this by constructing the ring of integers, as the constructor for a maximal OrderWithBasis calls this function to validate its input
    fn discriminant(&self) -> Integer;

    /// A list of self.n() elements which generate the ring of integers as a Z-module
    fn integral_basis(&self) -> Vec<Self::Set>;

    fn ring_of_integers<'a>(&'a self) -> OrderWithBasis<Self, &'a Self, true> {
        OrderWithBasis::new_maximal_unchecked(self, self.integral_basis())
    }
    fn into_ring_of_integers(self) -> OrderWithBasis<Self, Self, true> {
        let basis = self.integral_basis();
        OrderWithBasis::new_maximal_unchecked(self, basis)
    }

    fn order<'a>(
        &'a self,
        basis: Vec<Self::Set>,
    ) -> Result<OrderWithBasis<Self, &'a Self, false>, String> {
        OrderWithBasis::new(self, basis)
    }
    fn into_order(
        self,
        basis: Vec<Self::Set>,
    ) -> Result<OrderWithBasis<Self, Self, false>, String> {
        OrderWithBasis::new(self, basis)
    }

    /// The LCM of the denominators of the coefficients of the minimal polynomial of a.
    ///
    /// It may well be >1 even when the element a is an algebraic integer.
    fn min_poly_denominator_lcm(&self, a: &Self::Set) -> Integer {
        Integer::lcm_list(
            self.inbound_finite_dimensional_rational_extension()
                .min_poly(a)
                .coeffs()
                .into_iter()
                .map(|c| Integer::from(c.denominator()))
                .collect(),
        )
    }

    /// A scalar multiple of $a$ which is an algebraic integer.
    ///
    /// It need not return $a$ itself when $a$ is already an algebraic integer.
    fn integral_multiple(&self, a: &Self::Set) -> Self::Set {
        let m = self.min_poly_denominator_lcm(a);
        let b = self.mul(&self.try_from_rat(&Rational::from(m)).unwrap(), a);
        debug_assert!(self.is_algebraic_integer(&b));
        b
    }
}

pub trait AlgebraicIntegerRingSignature<K: AlgebraicNumberFieldSignature>:
    DedekindDomainSignature + CharZeroRingSignature
{
    fn n(&self) -> usize {
        self.anf().n()
    }

    fn anf(&self) -> &K;

    /// A list of self.n() elements which generate this ring as a Z-module
    fn integral_basis(&self) -> Vec<Self::Set>;

    fn to_anf(&self, x: &Self::Set) -> K::Set;

    fn try_from_anf(&self, y: &K::Set) -> Option<Self::Set>;

    fn order<'a>(&'a self) -> OrderWithBasis<K, &'a K, true> {
        OrderWithBasis::new_maximal_unchecked(
            self.anf(),
            self.integral_basis()
                .into_iter()
                .map(|v| self.outbound_roi_to_anf_inclusion().image(&v))
                .collect(),
        )
    }

    fn into_outbound_roi_to_anf_inclusion(
        self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<K, Self, Self> {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(self)
    }

    fn outbound_roi_to_anf_inclusion<'a>(
        &'a self,
    ) -> RingOfIntegersToAlgebraicNumberFieldInclusion<K, Self, &'a Self> {
        RingOfIntegersToAlgebraicNumberFieldInclusion::from_ring_of_integers(self)
    }

    /// The conductor of an order is an ideal of the order.
    ///
    /// It consists of all elements of the order such that multiplying by any element of the ring of integers produces an element of the order.
    ///
    /// Equivalently, it consists of all elements of the ring of integers such that multiplying by any element of the ring of integers produces an element of the order.
    ///
    /// The returned value represents an element of `order.ideals()` _not_ of `self.ideals()`.
    fn conductor<KB: BorrowedStructure<K>, const MAXIMAL: bool>(
        &self,
        order: &OrderWithBasis<K, KB, MAXIMAL>,
    ) -> OrderIdeal {
        todo!()
        // order
        //     .ideals()
        //     .outbound_sublattices_inclusion()
        //     .try_preimage(&order.quotient_sublattice(a, b))
        //     .unwrap()
    }
}

#[derive(Debug, Clone)]
pub struct RingOfIntegersToAlgebraicNumberFieldInclusion<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> {
    _roi: PhantomData<R>,
    roi: RB,
    _anf: PhantomData<K>,
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
{
    pub fn from_ring_of_integers(roi: RB) -> Self {
        Self {
            _roi: PhantomData,
            _anf: PhantomData,
            roi,
        }
    }

    pub fn roi(&self) -> &R {
        self.roi.borrow()
    }

    pub fn anf(&self) -> &K {
        self.roi().anf()
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> Morphism<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
{
    fn domain(&self) -> &R {
        self.roi()
    }

    fn range(&self) -> &K {
        self.anf()
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> Function<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
{
    fn image(&self, x: &<R as SetSignature>::Set) -> <K as SetSignature>::Set {
        self.roi().to_anf(x)
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> InjectiveFunction<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
{
    fn try_preimage(&self, y: &<K as SetSignature>::Set) -> Option<<R as SetSignature>::Set> {
        self.roi().try_from_anf(y)
    }
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> RingHomomorphism<R, K> for RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
{
}

impl<
    K: AlgebraicNumberFieldSignature,
    R: AlgebraicIntegerRingSignature<K>,
    RB: BorrowedStructure<R>,
> RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, RB>
{
    pub fn zq_extension<'a>(
        &'a self,
    ) -> RingOfIntegersExtension<
        K,
        R,
        &'a R,
        RingOfIntegersToAlgebraicNumberFieldInclusion<K, R, &'a R>,
        IntegerIdealsStructure<IntegerCanonicalStructure>,
        <R as RingToIdealsSignature>::Ideals<&'a R>,
    >
    where
        R: RingToIdealsSignature,
        <R as RingToIdealsSignature>::Ideals<&'a R>: DedekindDomainIdealsSignature<R, &'a R>,
    {
        let ideals_z = Integer::structure().into_ideals();
        let ideals_r = self.domain().ideals();
        RingOfIntegersExtension::new_integer_extension(
            RingOfIntegersToAlgebraicNumberFieldInclusion {
                _roi: PhantomData,
                _anf: PhantomData,
                roi: self.domain(),
            },
            ideals_z,
            ideals_r,
        )
    }
}

// pub trait AlgebraicIntegerRingInAlgebraicNumberFieldSignature:
//     RingHomomorphism<Self::RingOfIntegers, Self::AlgebraicNumberField>
//     + InjectiveFunction<Self::RingOfIntegers, Self::AlgebraicNumberField>
// {
//     type AlgebraicNumberField: AlgebraicNumberFieldSignature;
//     type RingOfIntegers: AlgebraicIntegerRingSignature<Self::AlgebraicNumberField>;

//     fn roi(&self) -> &Self::RingOfIntegers {
//         self.domain()
//     }

//     fn anf(&self) -> &Self::AlgebraicNumberField {
//         self.range()
//     }

//     fn discriminant(&self) -> Integer;

//     fn roi_to_anf(
//         &self,
//         x: &<Self::RingOfIntegers as SetSignature>::Set,
//     ) -> <Self::AlgebraicNumberField as SetSignature>::Set;

//     fn try_anf_to_roi(
//         &self,
//         y: &<Self::AlgebraicNumberField as SetSignature>::Set,
//     ) -> Option<<Self::RingOfIntegers as SetSignature>::Set>;
// }

mod anf_inclusion {
    use crate::{matrix::Matrix, structure::MetaCharZeroRing};

    use super::*;

    #[derive(Debug, Clone)]
    pub struct AlgebraicNumberFieldFullRankSublatticeWithBasisInclusion<
        K: AlgebraicNumberFieldSignature,
        Lat: FullRankSublatticeWithBasisSignature<K>,
        LatB: BorrowedStructure<Lat>,
    > {
        _k: PhantomData<K>,
        _lat: PhantomData<Lat>,
        sublattice: LatB,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        Lat: FullRankSublatticeWithBasisSignature<K>,
        LatB: BorrowedStructure<Lat>,
    > AlgebraicNumberFieldFullRankSublatticeWithBasisInclusion<K, Lat, LatB>
    {
        pub fn new(sublattice: LatB) -> Self {
            Self {
                _k: PhantomData,
                _lat: PhantomData,
                sublattice,
            }
        }

        pub fn sublattice(&self) -> &Lat {
            self.sublattice.borrow()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        Lat: FullRankSublatticeWithBasisSignature<K>,
        LatB: BorrowedStructure<Lat>,
    > Morphism<Lat, K> for AlgebraicNumberFieldFullRankSublatticeWithBasisInclusion<K, Lat, LatB>
    {
        fn domain(&self) -> &Lat {
            self.sublattice()
        }

        fn range(&self) -> &K {
            self.sublattice().anf()
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        Lat: FullRankSublatticeWithBasisSignature<K>,
        LatB: BorrowedStructure<Lat>,
    > Function<Lat, K> for AlgebraicNumberFieldFullRankSublatticeWithBasisInclusion<K, Lat, LatB>
    {
        fn image(&self, x: &Vec<Integer>) -> <K as SetSignature>::Set {
            debug_assert!(self.sublattice().is_element(x).is_ok());
            let k = self.sublattice().anf();
            let n = k.n();
            debug_assert_eq!(n, x.len());
            k.sum(
                (0..n)
                    .map(|i| k.mul(&k.from_int(&x[i]), self.sublattice().basis_vector(i)))
                    .collect(),
            )
        }
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        Lat: FullRankSublatticeWithBasisSignature<K>,
        LatB: BorrowedStructure<Lat>,
    > InjectiveFunction<Lat, K>
        for AlgebraicNumberFieldFullRankSublatticeWithBasisInclusion<K, Lat, LatB>
    {
        fn try_preimage(&self, y: &<K as SetSignature>::Set) -> Option<Vec<Integer>> {
            let k = self.sublattice().anf();
            let n = k.n();
            debug_assert!(k.is_element(y).is_ok());
            let mat = Matrix::join_cols(
                n,
                (0..n)
                    .map(|i| {
                        k.inbound_finite_dimensional_rational_extension()
                            .to_col(self.sublattice().basis_vector(i))
                    })
                    .collect(),
            );
            let y = k.inbound_finite_dimensional_rational_extension().to_vec(y);
            let x_rat = mat.col_solve(&y)?;
            let mut x_int = Vec::with_capacity(n);
            for c_rat in x_rat {
                x_int.push(c_rat.try_to_int()?);
            }
            Some(x_int)
        }
    }
}

mod sublattice_inclusion {
    use super::*;
    use algebraeon_sets::structure::BorrowedStructure;
    use std::marker::PhantomData;

    pub struct SublatticeInclusion<
        K: AlgebraicNumberFieldSignature,
        Lat: FullRankSublatticeWithBasisSignature<K>,
        LatB: BorrowedStructure<Lat>,
        Sublat: FullRankSublatticeWithBasisSignature<K>,
        SublatB: BorrowedStructure<Sublat>,
    > {
        _k: PhantomData<K>,
        _lat: PhantomData<Lat>,
        _sublat: PhantomData<Sublat>,
        lattice: LatB,
        sublattice: SublatB,
    }

    impl<
        K: AlgebraicNumberFieldSignature,
        Lat: FullRankSublatticeWithBasisSignature<K>,
        LatB: BorrowedStructure<Lat>,
        Sublat: FullRankSublatticeWithBasisSignature<K>,
        SublatB: BorrowedStructure<Sublat>,
    > SublatticeInclusion<K, Lat, LatB, Sublat, SublatB>
    {
        pub fn lattice(&self) -> &Lat {
            self.lattice.borrow()
        }

        pub fn sublattice(&self) -> &Sublat {
            self.sublattice.borrow()
        }

        fn check(&self) -> Result<(), String> {
            if self.lattice().contains_sublattice(self.sublattice()) {
                Ok(())
            } else {
                Err("Lattice does not contain sublattice".to_string())
            }
        }

        fn new_impl(lattice: LatB, sublattice: SublatB) -> Self {
            Self {
                _k: PhantomData,
                _lat: PhantomData,
                _sublat: PhantomData,
                lattice,
                sublattice,
            }
        }

        pub fn new(lattice: LatB, sublattice: SublatB) -> Option<Self> {
            let s = Self::new_impl(lattice, sublattice);
            if s.check().is_ok() { Some(s) } else { None }
        }

        pub fn new_unchecked(lattice: LatB, sublattice: SublatB) -> Self {
            let s = Self::new_impl(lattice, sublattice);
            #[cfg(debug_assertions)]
            s.check().unwrap();
            s
        }
    }
}

pub trait FullRankSublatticeWithBasisSignature<K: AlgebraicNumberFieldSignature>:
    AdditiveGroupSignature<Set = Vec<Integer>>
{
    fn anf(&self) -> &K;

    fn basis(&self) -> &Vec<K::Set>;

    fn basis_vector(&self, i: usize) -> &K::Set {
        debug_assert!(i < self.n());
        self.basis().get(i).unwrap()
    }

    fn n(&self) -> usize {
        debug_assert_eq!(self.anf().n(), self.basis().len());
        self.anf().n()
    }

    fn free_lattice_restructure(
        &self,
    ) -> FinitelyFreeModuleStructure<IntegerCanonicalStructure, &'static IntegerCanonicalStructure>
    {
        Integer::structure_ref().free_module(self.n())
    }

    fn contains_element(&self, p: &K::Set) -> bool {
        self.outbound_order_to_anf_inclusion()
            .try_preimage(p)
            .is_some()
    }

    fn contains_sublattice<Sublat: FullRankSublatticeWithBasisSignature<K>>(
        &self,
        sublattice: &Sublat,
    ) -> bool {
        sublattice
            .basis()
            .iter()
            .all(|sublat_basis_vector| self.contains_element(sublat_basis_vector))
    }

    fn discriminant(&self) -> Rational {
        self.anf()
            .inbound_finite_dimensional_rational_extension()
            .discriminant(self.basis())
    }

    fn into_outbound_order_to_anf_inclusion(
        self,
    ) -> anf_inclusion::AlgebraicNumberFieldFullRankSublatticeWithBasisInclusion<K, Self, Self>
    {
        anf_inclusion::AlgebraicNumberFieldFullRankSublatticeWithBasisInclusion::new(self)
    }

    fn outbound_order_to_anf_inclusion<'a>(
        &'a self,
    ) -> anf_inclusion::AlgebraicNumberFieldFullRankSublatticeWithBasisInclusion<K, Self, &'a Self>
    {
        anf_inclusion::AlgebraicNumberFieldFullRankSublatticeWithBasisInclusion::new(self)
    }

    fn into_inbound_sublattice_inclusion<
        Sublat: FullRankSublatticeWithBasisSignature<K>,
        SublatB: BorrowedStructure<Sublat>,
    >(
        self,
        sublattice: SublatB,
    ) -> Option<sublattice_inclusion::SublatticeInclusion<K, Self, Self, Sublat, SublatB>> {
        sublattice_inclusion::SublatticeInclusion::new(self, sublattice)
    }

    fn inbound_sublattice_inclusion<
        'a,
        Sublat: FullRankSublatticeWithBasisSignature<K>,
        SublatB: BorrowedStructure<Sublat>,
    >(
        &'a self,
        sublattice: SublatB,
    ) -> Option<sublattice_inclusion::SublatticeInclusion<K, Self, &'a Self, Sublat, SublatB>> {
        sublattice_inclusion::SublatticeInclusion::new(self, sublattice)
    }

    fn into_inbound_sublattice_inclusion_unchecked<
        Sublat: FullRankSublatticeWithBasisSignature<K>,
        SublatB: BorrowedStructure<Sublat>,
    >(
        self,
        sublattice: SublatB,
    ) -> sublattice_inclusion::SublatticeInclusion<K, Self, Self, Sublat, SublatB> {
        sublattice_inclusion::SublatticeInclusion::new_unchecked(self, sublattice)
    }

    fn inbound_sublattice_inclusion_unchecked<
        'a,
        Sublat: FullRankSublatticeWithBasisSignature<K>,
        SublatB: BorrowedStructure<Sublat>,
    >(
        &'a self,
        sublattice: SublatB,
    ) -> sublattice_inclusion::SublatticeInclusion<K, Self, &'a Self, Sublat, SublatB> {
        sublattice_inclusion::SublatticeInclusion::new_unchecked(self, sublattice)
    }
}
