use algebraeon_structures::*;
use itertools::Itertools;

/// The set of all endofunctions on a finite set X: functions X → X
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteSetEndofunctions<X: FiniteSetSignature + EqSignature> {
    set: X,
}

impl<X: FiniteSetSignature + EqSignature> FiniteSetEndofunctions<X> {
    pub fn new(set: X) -> Self {
        Self { set }
    }
}

impl<X: FiniteSetSignature + EqSignature> Signature for FiniteSetEndofunctions<X> {}

impl<X: FiniteSetSignature + EqSignature> SetSignature for FiniteSetEndofunctions<X> {
    type Elem = Vec<X::Elem>;

    fn validate_element(&self, f: &Self::Elem) -> Result<(), String> {
        if Natural::from(f.len()) != self.set.nat_size() {
            return Err("Function must have one value per element in the domain.".to_string());
        }
        for y in f {
            self.set.validate_element(y)?;
        }
        Ok(())
    }
}

impl<X: FiniteSetSignature + EqSignature> CountableSetSignature for FiniteSetEndofunctions<X> {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        let n: usize = self.set.nat_size().try_into().unwrap_or(usize::MAX);
        (0..n)
            .map(|_| self.set.list_all_elements())
            .multi_cartesian_product()
    }
}

impl<X: FiniteSetSignature + EqSignature> FiniteSetSignature for FiniteSetEndofunctions<X> {}
