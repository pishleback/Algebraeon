use crate::structure::{
    AdditionSignature, AdditiveGroupSignature, AdditiveMonoidSignature,
    CancellativeAdditionSignature, CancellativeMultiplicationSignature,
    CommutativeMultiplicationSignature, LeftDistributiveMultiplicationOverAddition,
    MultiplicationSignature, MultiplicativeAbsorptionMonoidSignature,
    MultiplicativeMonoidSignature, OneSignature, RightDistributiveMultiplicationOverAddition,
    RingSignature, RinglikeSpecializationSignature, SemiRingSignature, TryNegateSignature,
    ZeroSignature,
};
use algebraeon_sets::structure::{EqSignature, SetSignature, Signature, ToStringSignature};

fn xgcd(mut x: u64, mut y: u64) -> (u64, i64, i64) {
    let mut pa = 1;
    let mut a = 0;
    let mut pb = 0;
    let mut b = 1;
    while y != 0 {
        let (q, r) = (x / y, x % y);
        let new_a = pa - q as i64 * a;
        (a, pa) = (new_a, a);
        let new_b = pb - q as i64 * b;
        (b, pb) = (new_b, b);
        (x, y) = (y, r);
    }
    (x, pa, pb)
}

// Return the inverse of x modulo n if it exists
fn inv_mod_n(x: u64, n: u64) -> Option<u64> {
    let (g, a, _) = xgcd(x, n);
    if g == 1 {
        Some(a.rem_euclid(n as i64).try_into().unwrap())
    } else {
        None
    }
}

/// The ring of integers modulo odd n
/// where elements x are represented in Montgomery form as a 0 <= (x*r % n) < n
/// where r is the smallest power of 2 such that r > n.
#[derive(Debug, Clone)]
pub struct MontgomeryModuloOddStructure {
    n: u64,
    r: u64,
    r_pow: u32,
    inv_cache: Option<Vec<u64>>, // len = n if not None
}

impl MontgomeryModuloOddStructure {
    pub fn new_without_inv_cache_unchecked(n: u64) -> Self {
        debug_assert_eq!(n % 2, 1);
        // Need a bound such that all operations remain in u64
        debug_assert!(n < (1u64 << 31));
        let r = n.next_power_of_two();
        let r_pow = r.ilog2();
        debug_assert_eq!(1 << r_pow, r);
        debug_assert_eq!(r.count_ones(), 1);
        debug_assert!(n < r);
        Self {
            n,
            r,
            r_pow,
            inv_cache: None,
        }
    }
}

impl PartialEq for MontgomeryModuloOddStructure {
    fn eq(&self, other: &Self) -> bool {
        self.n == other.n
    }
}

impl Eq for MontgomeryModuloOddStructure {}

impl Signature for MontgomeryModuloOddStructure {}

impl SetSignature for MontgomeryModuloOddStructure {
    type Elem = u64; // 0 <= x < n

    fn validate_element(&self, x: &Self::Elem) -> Result<(), String> {
        if x >= &self.n {
            return Err(format!(
                "{x} is out of the valid range 0 <= {x} < {}",
                self.n
            ));
        }
        Ok(())
    }
}

impl ToStringSignature for MontgomeryModuloOddStructure {
    fn to_string(&self, elem: &Self::Elem) -> String {
        todo!()
    }
}

impl EqSignature for MontgomeryModuloOddStructure {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        #[cfg(debug_assertions)]
        {
            self.validate_element(a).unwrap();
            self.validate_element(b).unwrap();
        }
        a == b
    }
}

impl RinglikeSpecializationSignature for MontgomeryModuloOddStructure {}

impl ZeroSignature for MontgomeryModuloOddStructure {
    fn zero(&self) -> Self::Elem {
        0
    }
}

impl AdditionSignature for MontgomeryModuloOddStructure {
    fn add(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        #[cfg(debug_assertions)]
        {
            self.validate_element(a).unwrap();
            self.validate_element(b).unwrap();
        }
        let s = a + b;
        let s = if s >= self.n { s - self.n } else { s };
        #[cfg(debug_assertions)]
        self.validate_element(&s).unwrap();
        s
    }
}

impl TryNegateSignature for MontgomeryModuloOddStructure {
    fn try_neg(&self, a: &Self::Elem) -> Option<Self::Elem> {
        Some(self.neg(a))
    }
}

impl CancellativeAdditionSignature for MontgomeryModuloOddStructure {
    fn try_sub(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        Some(self.add(a, &self.neg(b)))
    }
}

impl AdditiveMonoidSignature for MontgomeryModuloOddStructure {}

impl AdditiveGroupSignature for MontgomeryModuloOddStructure {
    fn neg(&self, a: &Self::Elem) -> Self::Elem {
        #[cfg(debug_assertions)]
        self.validate_element(a).unwrap();
        let b = if *a == 0 { *a } else { self.n - a };
        #[cfg(debug_assertions)]
        self.validate_element(&b).unwrap();
        b
    }
}

impl OneSignature for MontgomeryModuloOddStructure {
    fn one(&self) -> Self::Elem {
        1
    }
}

impl MultiplicationSignature for MontgomeryModuloOddStructure {
    fn mul(&self, a: &Self::Elem, b: &Self::Elem) -> Self::Elem {
        // Use the REDC algorithm
        // https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
        // This algorithm is the whole reason why Montgomery form is useful.
        todo!()
    }
}

impl MultiplicativeMonoidSignature for MontgomeryModuloOddStructure {}

impl CommutativeMultiplicationSignature for MontgomeryModuloOddStructure {}

impl CancellativeMultiplicationSignature for MontgomeryModuloOddStructure {
    fn try_divide(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Self::Elem> {
        todo!()
    }
}

impl RightDistributiveMultiplicationOverAddition for MontgomeryModuloOddStructure {}

impl LeftDistributiveMultiplicationOverAddition for MontgomeryModuloOddStructure {}

impl MultiplicativeAbsorptionMonoidSignature for MontgomeryModuloOddStructure {}

impl SemiRingSignature for MontgomeryModuloOddStructure {}

impl RingSignature for MontgomeryModuloOddStructure {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn additive() {
        let ring = MontgomeryModuloOddStructure::new_without_inv_cache_unchecked(7);
        assert!(ring.equal(&ring.zero(), &0));
        assert!(ring.equal(&ring.add(&1, &2), &3));
        assert!(ring.equal(&ring.add(&3, &4), &0));
        assert!(ring.equal(&ring.add(&5, &6), &4));
        assert!(ring.equal(&ring.neg(&0), &0));
        assert!(ring.equal(&ring.neg(&4), &3));
    }

    #[test]
    fn test() {
        let ring = MontgomeryModuloOddStructure::new_without_inv_cache_unchecked(13);

        println!("{:?}", ring);

        // println!("{:?}", ring.to_string(&3));
    }
}
