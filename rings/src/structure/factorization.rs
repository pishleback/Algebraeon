use crate::structure::{
    AdditiveMonoidSignature, FavoriteAssociateSignature, MultiplicativeMonoidSignature,
    MultiplicativeMonoidWithZeroSignature, SetWithZeroSignature,
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
    Object: MultiplicativeMonoidWithZeroSignature + FavoriteAssociateSignature + EqSignature,
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
    Object: MultiplicativeMonoidWithZeroSignature + FavoriteAssociateSignature + EqSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: AdditiveMonoidSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    pub fn objects(&self) -> &Object {
        self.objects.borrow()
    }

    pub fn exponents(&self) -> &Exponent {
        self.exponents.borrow()
    }
}

impl<
    Object: MultiplicativeMonoidWithZeroSignature + FavoriteAssociateSignature + EqSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: AdditiveMonoidSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> Signature for FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
}

impl<
    Object: MultiplicativeMonoidWithZeroSignature + FavoriteAssociateSignature + EqSignature,
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
    Object: MultiplicativeMonoidWithZeroSignature + FavoriteAssociateSignature + EqSignature,
    ObjectB: BorrowedStructure<Object>,
    Exponent: AdditiveMonoidSignature + OrdSignature,
    ExponentB: BorrowedStructure<Exponent>,
> EqSignature for FactoringStructure<Object, ObjectB, Exponent, ExponentB>
{
    fn equal(&self, a: &Self::Set, b: &Self::Set) -> bool {
        todo!()
    }
}

impl<
    Object: MultiplicativeMonoidWithZeroSignature + FavoriteAssociateSignature + EqSignature,
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
    Object: MultiplicativeMonoidWithZeroSignature + FavoriteAssociateSignature + EqSignature,
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
    Object: MultiplicativeMonoidWithZeroSignature + FavoriteAssociateSignature + EqSignature,
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
        match (a, b) {
            (Factored::NonZero(a), Factored::NonZero(b)) => {
                todo!();
                todo!()
            }
            _ => Factored::Zero,
        }
    }
}
