use crate::structure::{
    AdditiveMonoidEqSignature, FavoriteAssociateSignature, FieldSignature,
    MultiplicativeGroupSignature, MultiplicativeIntegralMonoidSignature,
    MultiplicativeMonoidSignature, MultiplicativeMonoidUnitsSignature,
    MultiplicativeMonoidWithZeroSignature, SemiRingSignature, SetWithZeroSignature,
};
use algebraeon_nzq::{Natural, NaturalCanonicalStructure};
use algebraeon_sets::structure::{
    BorrowedStructure, EqSignature, MetaType, OrdSignature, SetSignature, Signature,
};
use std::fmt::Debug;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct NonZeroFactored<ObjectSet: Debug + Clone, ExponentSet: Debug + Clone> {
    pub unit: ObjectSet,
    pub powers: Vec<(ObjectSet, ExponentSet)>,
}

#[derive(Debug, Clone)]
pub enum Factored<ObjectSet: Debug + Clone, ExponentSet: Debug + Clone> {
    Zero,
    NonZero(NonZeroFactored<ObjectSet, ExponentSet>),
}

impl<ObjectSet: Debug + Clone, ExponentSet: Debug + Clone> Factored<ObjectSet, ExponentSet> {
    pub fn powers<'a>(&'a self) -> Option<&'a Vec<(ObjectSet, ExponentSet)>> {
        match self {
            Factored::Zero => None,
            Factored::NonZero(a) => Some(&a.powers),
        }
    }

    pub fn into_powers(self) -> Option<Vec<(ObjectSet, ExponentSet)>> {
        match self {
            Factored::Zero => None,
            Factored::NonZero(a) => Some(a.powers),
        }
    }

    pub fn unit_and_powers<'a>(
        &'a self,
    ) -> Option<(&'a ObjectSet, &'a Vec<(ObjectSet, ExponentSet)>)> {
        match self {
            Factored::Zero => None,
            Factored::NonZero(a) => Some((&a.unit, &a.powers)),
        }
    }

    pub fn into_unit_and_powers(self) -> Option<(ObjectSet, Vec<(ObjectSet, ExponentSet)>)> {
        match self {
            Factored::Zero => None,
            Factored::NonZero(a) => Some((a.unit, a.powers)),
        }
    }

    pub fn distinct_irreducibles<'a>(&'a self) -> Option<Vec<&'a ObjectSet>> {
        match self {
            Factored::Zero => None,
            Factored::NonZero(a) => Some(a.powers.iter().map(|(p, _)| p).collect()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FactoringStructure<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: SemiRingSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> {
    _objects: PhantomData<Object>,
    objects: ObjectB,
    _exponents: PhantomData<Exponent>,
    exponents: ExponentB,
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: SemiRingSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    pub fn new(objects: ObjectB, exponents: ExponentB) -> Self {
        Self {
            _objects: PhantomData,
            objects,
            _exponents: PhantomData,
            exponents,
        }
    }

    pub fn objects(&self) -> &Object {
        self.objects.borrow()
    }

    pub fn exponents(&self) -> &Exponent {
        self.exponents.borrow()
    }
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: SemiRingSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> Signature for FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: SemiRingSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> SetSignature for FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    type Set = Factored<Object::Set, Exponent::Set>;

    fn is_element(&self, x: &Self::Set) -> Result<(), String> {
        match x {
            Factored::Zero => Ok(()),
            Factored::NonZero(x) => {
                self.objects().units().is_element(&x.unit)?;
                for (prime, exponent) in &x.powers {
                    if !self.objects().is_fav_assoc(prime) {
                        return Err("Factors should be their favorite associate".to_string());
                    }
                    if self.objects().try_is_irreducible(prime) == Some(false) {
                        return Err("Failed to be irreducible".to_string());
                    }
                    self.exponents().is_element(exponent)?;
                    if self.exponents().is_zero(exponent) {
                        return Err("Exponent is zero".to_string());
                    }
                }

                for i in 0..x.powers.len() {
                    for j in 0..i {
                        let (pi, _ki) = &x.powers[i];
                        let (pj, _kj) = &x.powers[j];
                        if self.objects().equal(pi, pj) {
                            return Err("primes must be distinct".to_string());
                        }
                    }
                }

                Ok(())
            }
        }
    }
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: SemiRingSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> EqSignature for FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        #[cfg(debug_assertions)]
        {
            self.is_element(a).unwrap();
            self.is_element(b).unwrap();
        }
        todo!()
    }
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: SemiRingSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> SetWithZeroSignature for FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    fn zero(&self) -> Self::Set {
        Factored::Zero
    }
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: SemiRingSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> MultiplicativeMonoidSignature for FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    fn one(&self) -> Self::Set {
        Factored::NonZero(NonZeroFactored {
            unit: self.objects().one(),
            powers: vec![],
        })
    }

    fn mul(&self, a: &Self::Set, b: &Self::Set) -> Self::Set {
        #[cfg(debug_assertions)]
        {
            self.is_element(a).unwrap();
            self.is_element(b).unwrap();
        }
        let mut s = a.clone();
        self.mul_mut(&mut s, b);
        s
    }
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: SemiRingSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> MultiplicativeMonoidUnitsSignature for FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    fn try_inv(&self, a: &Self::Set) -> Option<Self::Set> {
        #[cfg(debug_assertions)]
        self.is_element(a).unwrap();
        match a {
            Factored::Zero => None,
            Factored::NonZero(a) => Some(Factored::NonZero(NonZeroFactored {
                unit: self.objects().units().inv(&a.unit),
                powers: a
                    .powers
                    .iter()
                    .map(|(p, k)| self.exponents().try_neg(k).map(|neg_k| (p.clone(), neg_k)))
                    .collect::<Option<Vec<_>>>()?,
            })),
        }
    }
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: SemiRingSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> MultiplicativeIntegralMonoidSignature
    for FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    fn try_div(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        let mut s = Some(a.clone());
        self.try_div_mut_impl(&mut s, b);
        s
    }
}

// These methods are all completely unchecked, to be used by things which construct factorizations.
impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: SemiRingSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    pub(crate) fn new_irreducible_impl(
        &self,
        p: Object::Set,
    ) -> Factored<Object::Set, Exponent::Set> {
        Factored::NonZero(NonZeroFactored {
            unit: self.objects().one(),
            powers: vec![(p, self.exponents().one())],
        })
    }

    pub(crate) fn new_unit_impl(&self, unit: Object::Set) -> Factored<Object::Set, Exponent::Set> {
        Factored::NonZero(NonZeroFactored {
            unit,
            powers: vec![],
        })
    }

    pub(crate) fn new_unit_and_powers_impl(
        &self,
        unit: Object::Set,
        powers: Vec<(Object::Set, Exponent::Set)>,
    ) -> Factored<Object::Set, Exponent::Set> {
        Factored::NonZero(NonZeroFactored { unit, powers })
    }

    pub(crate) fn mul_mut_impl(
        &self,
        a: &mut Factored<Object::Set, Exponent::Set>,
        b: &Factored<Object::Set, Exponent::Set>,
    ) {
        match b {
            Factored::Zero => {
                *a = Factored::Zero;
            }
            Factored::NonZero(b) => match a {
                Factored::Zero => {}
                Factored::NonZero(a) => {
                    a.unit = self.objects().units().mul(&a.unit, &b.unit);
                    for (b_p, b_k) in &b.powers {
                        debug_assert!(self.objects().is_fav_assoc(b_p));
                        'A_LOOP: {
                            for (a_p, a_k) in &mut a.powers {
                                debug_assert!(self.objects().is_fav_assoc(a_p));
                                if self.objects().equal(a_p, b_p) {
                                    *a_k = self.exponents().add(a_k, b_k);
                                    break 'A_LOOP;
                                }
                            }
                            a.powers.push((b_p.clone(), b_k.clone()));
                        }
                    }
                    a.powers.retain(|(_, k)| !self.exponents().is_zero(k));
                }
            },
        }
    }

    pub(crate) fn try_div_mut_impl(
        &self,
        a: &mut Option<Factored<Object::Set, Exponent::Set>>,
        b: &Factored<Object::Set, Exponent::Set>,
    ) {
        match b {
            Factored::Zero => {
                *a = None;
            }
            Factored::NonZero(b) => match a {
                None => {}
                Some(Factored::Zero) => {}
                Some(Factored::NonZero(a_nz)) => {
                    a_nz.unit = self
                        .objects()
                        .units()
                        .mul(&a_nz.unit, &self.objects().units().inv(&b.unit));
                    for (b_p, b_k) in &b.powers {
                        debug_assert!(self.objects().is_fav_assoc(b_p));
                        'A_LOOP: {
                            for (a_p, a_k) in &mut a_nz.powers {
                                debug_assert!(self.objects().is_fav_assoc(a_p));
                                if self.objects().equal(a_p, b_p) {
                                    if let Some(a_k_minus_b_k) = self.exponents().try_sub(a_k, b_k)
                                    {
                                        *a_k = a_k_minus_b_k;
                                    } else {
                                        *a = None;
                                        return;
                                    }
                                    break 'A_LOOP;
                                }
                            }
                            if let Some(b_k_neg) = self.exponents().try_neg(b_k) {
                                a_nz.powers.push((b_p.clone(), b_k_neg));
                            } else {
                                *a = None;
                                return;
                            }
                        }
                    }
                    a_nz.powers.retain(|(_, k)| !self.exponents().is_zero(k));
                }
            },
        }
    }

    pub(crate) fn is_irreducible_impl(&self, a: &Factored<Object::Set, Exponent::Set>) -> bool {
        match a {
            Factored::Zero => {}
            Factored::NonZero(a) => {
                if a.powers.len() == 1 {
                    let (_, k) = &a.powers[0];
                    if self.exponents().equal(k, &self.exponents().one()) {
                        return true;
                    }
                }
            }
        }
        false
    }
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: SemiRingSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    pub fn new_unit_and_powers_unchecked(
        &self,
        unit: Object::Set,
        powers: Vec<(Object::Set, Exponent::Set)>,
    ) -> Factored<Object::Set, Exponent::Set> {
        let f = self.new_unit_and_powers_impl(unit, powers);
        #[cfg(debug_assertions)]
        self.is_element(&f).unwrap();
        f
    }

    pub fn expand_squarefree(&self, a: &Factored<Object::Set, Exponent::Set>) -> Object::Set {
        #[cfg(debug_assertions)]
        self.is_element(a).unwrap();
        match a {
            Factored::Zero => self.objects().zero(),
            Factored::NonZero(a) => {
                let mut prod = a.unit.clone();
                for (p, _) in &a.powers {
                    self.objects().mul_mut(&mut prod, p);
                }
                prod
            }
        }
    }
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    ExponentB: BorrowedStructure<NaturalCanonicalStructure>,
> FactoringStructure<Object, ObjectB, NaturalCanonicalStructure, ExponentB>
{
    pub fn expand(&self, a: &Factored<Object::Set, Natural>) -> Object::Set {
        #[cfg(debug_assertions)]
        self.is_element(a).unwrap();
        match a {
            Factored::Zero => self.objects().zero(),
            Factored::NonZero(a) => {
                let mut prod = a.unit.clone();
                for (p, k) in &a.powers {
                    self.objects()
                        .mul_mut(&mut prod, &self.objects().nat_pow(p, k));
                }
                prod
            }
        }
    }

    pub fn is_irreducible(&self, a: &Factored<Object::Set, Natural>) -> bool {
        #[cfg(debug_assertions)]
        self.is_element(a).unwrap();
        self.is_irreducible_impl(a)
    }

    /// Return an iterator over all divisors of a factorization
    pub fn divisors<'a>(
        &'a self,
        a: &'a Factored<Object::Set, Natural>,
    ) -> Option<Box<dyn Iterator<Item = Object::Set> + 'a>> {
        match a {
            Factored::Zero => None,
            Factored::NonZero(a) => {
                let factors = &a.powers;
                if factors.is_empty() {
                    Some(Box::new(vec![self.objects().one()].into_iter()))
                } else {
                    let mut factor_powers = vec![];
                    for (p, k) in factors {
                        let j = factor_powers.len();
                        factor_powers.push(vec![]);
                        let mut p_pow = self.objects().one();
                        let mut i = Natural::from(0u8);
                        while &i <= k {
                            factor_powers[j].push(p_pow.clone());
                            p_pow = self.objects().mul(&p_pow, p);
                            i += Natural::from(1u8);
                        }
                    }

                    #[allow(clippy::redundant_closure_for_method_calls)]
                    Some(Box::new(
                        itertools::Itertools::multi_cartesian_product(
                            factor_powers.into_iter().map(|p_pows| p_pows.into_iter()),
                        )
                        .map(move |prime_power_factors| {
                            self.objects()
                                .product(prime_power_factors.iter().collect())
                                .clone()
                        }),
                    ))
                }
            }
        }
    }

    /// The number of divisors of a factorization
    pub fn count_divisors(&self, a: &Factored<Object::Set, Natural>) -> Option<Natural> {
        #[cfg(debug_assertions)]
        self.is_element(a).unwrap();
        match a {
            Factored::Zero => None,
            Factored::NonZero(a) => {
                let factors = &a.powers;
                let mut count = Natural::from(1u8);
                for (_p, k) in factors {
                    count *= k + Natural::ONE;
                }
                Some(count)
            }
        }
    }
}

pub trait UniqueFactorizationMonoidSignature:
    MultiplicativeMonoidWithZeroSignature + FavoriteAssociateSignature + EqSignature
{
    type FactoredExponent: SemiRingSignature + OrdSignature;

    fn factorization_exponents<'a>(&'a self) -> &'a Self::FactoredExponent;
    fn into_factorization_exponents(self) -> Self::FactoredExponent;

    fn factorizations<'a>(
        &'a self,
    ) -> FactoringStructure<Self, &'a Self, Self::FactoredExponent, &'a Self::FactoredExponent>
    {
        FactoringStructure::new(self, self.factorization_exponents())
    }
    fn into_factorizations(
        self,
    ) -> FactoringStructure<Self, Self, Self::FactoredExponent, Self::FactoredExponent> {
        FactoringStructure::new(self.clone(), self.into_factorization_exponents())
    }

    fn try_is_irreducible(&self, a: &Self::Set) -> Option<bool>;
}
pub trait MetaUniqueFactorizationMonoid: MetaType
where
    Self::Signature: UniqueFactorizationMonoidSignature,
{
    fn try_is_irreducible(&self) -> Option<bool> {
        Self::structure().try_is_irreducible(self)
    }
}
impl<R: MetaType> MetaUniqueFactorizationMonoid for R where
    Self::Signature: UniqueFactorizationMonoidSignature
{
}

pub trait FactoringMonoidSignature: UniqueFactorizationMonoidSignature {
    fn is_irreducible(&self, a: &Self::Set) -> bool {
        self.factorizations()
            .is_irreducible_impl(&self.factor_unchecked(a))
    }
    fn factor_unchecked(
        &self,
        a: &Self::Set,
    ) -> Factored<Self::Set, <Self::FactoredExponent as SetSignature>::Set>;
    fn factor(
        &self,
        a: &Self::Set,
    ) -> Factored<Self::Set, <Self::FactoredExponent as SetSignature>::Set> {
        let f = self.factor_unchecked(a);
        #[cfg(debug_assertions)]
        self.factorizations().is_element(&f).unwrap();
        f
    }
}
pub trait MetaFactoringMonoid: MetaType
where
    Self::Signature: FactoringMonoidSignature,
{
    fn is_irreducible(&self) -> bool {
        Self::structure().is_irreducible(self)
    }
    fn factor_unchecked(
        &self,
    ) -> Factored<
        Self,
        <<Self::Signature as UniqueFactorizationMonoidSignature>::FactoredExponent as SetSignature>::Set,
    >{
        Self::structure().factor_unchecked(self)
    }
    fn factor(
        &self,
    ) -> Factored<
        Self,
        <<Self::Signature as UniqueFactorizationMonoidSignature>::FactoredExponent as SetSignature>::Set,
    >{
        Self::structure().factor(self)
    }
}
impl<R: MetaType> MetaFactoringMonoid for R where Self::Signature: FactoringMonoidSignature {}

impl<FS: FieldSignature> UniqueFactorizationMonoidSignature for FS {
    type FactoredExponent = NaturalCanonicalStructure;

    fn factorization_exponents<'a>(&'a self) -> &'a Self::FactoredExponent {
        Natural::structure_ref()
    }

    fn into_factorization_exponents(self) -> Self::FactoredExponent {
        Natural::structure()
    }

    fn try_is_irreducible(&self, a: &Self::Set) -> Option<bool> {
        Some(false)
    }
}

impl<FS: FieldSignature> FactoringMonoidSignature for FS {
    fn factor_unchecked(
        &self,
        a: &Self::Set,
    ) -> Factored<Self::Set, <Self::FactoredExponent as SetSignature>::Set> {
        if self.is_zero(a) {
            Factored::Zero
        } else {
            Factored::NonZero(NonZeroFactored {
                unit: a.clone(),
                powers: vec![],
            })
        }
    }
}
