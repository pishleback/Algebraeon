use crate::structure::{
    AdditiveMonoidSignature, FavoriteAssociateSignature, MultiplicativeGroupSignature,
    MultiplicativeIntegralMonoidSignature, MultiplicativeMonoidSignature,
    MultiplicativeMonoidUnitsSignature, MultiplicativeMonoidWithZeroSignature,
    SetWithZeroSignature,
};
use algebraeon_sets::structure::{
    BorrowedStructure, EqSignature, OrdSignature, SetSignature, Signature,
};
use std::fmt::Debug;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct NonZeroFactored<ObjectSet: Debug + Clone, ExponentSet: Debug + Clone> {
    unit: ObjectSet,
    powers: Vec<(ObjectSet, ExponentSet)>,
}

#[derive(Debug, Clone)]
pub enum Factored<ObjectSet: Debug + Clone, ExponentSet: Debug + Clone> {
    Zero,
    NonZero(NonZeroFactored<ObjectSet, ExponentSet>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FactoringStructure<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: AdditiveMonoidSignature + OrdSignature,
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
    Exponent: AdditiveMonoidSignature + OrdSignature,
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
    Exponent: AdditiveMonoidSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> Signature for FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: AdditiveMonoidSignature + OrdSignature,
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
                    todo!("check prime is prime");
                    self.exponents().is_element(exponent)?;
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
    Exponent: AdditiveMonoidSignature + OrdSignature,
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
    Exponent: AdditiveMonoidSignature + OrdSignature,
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
    Exponent: AdditiveMonoidSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    fn mul_mut(
        &self,
        a: &mut Factored<Object::Set, Exponent::Set>,
        b: &Factored<Object::Set, Exponent::Set>,
    ) {
        #[cfg(debug_assertions)]
        {
            self.is_element(a).unwrap();
            self.is_element(b).unwrap();
        }
        match b {
            Factored::Zero => {
                *a = Factored::Zero;
            }
            Factored::NonZero(b) => match a {
                Factored::Zero => {}
                Factored::NonZero(a) => {
                    a.unit = self.objects().mul(&a.unit, &b.unit);
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
                }
            },
        }
    }
}

impl<
    Object: UniqueFactorizationMonoidSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: AdditiveMonoidSignature + OrdSignature,
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
    Exponent: AdditiveMonoidSignature + OrdSignature,
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
    Exponent: AdditiveMonoidSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> MultiplicativeIntegralMonoidSignature
    for FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    fn try_div(&self, a: &Self::Set, b: &Self::Set) -> Option<Self::Set> {
        todo!()
    }
}

pub trait UniqueFactorizationMonoidSignature:
    MultiplicativeMonoidWithZeroSignature + FavoriteAssociateSignature + EqSignature
{
    type FactoredExponent: AdditiveMonoidSignature + OrdSignature;

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

    fn is_irreducible(&self, a: &Self::Set) -> bool;
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
