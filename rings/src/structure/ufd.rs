use crate::structure::{
    AdditiveMonoidEqSignature, Factored, MultiplicativeMonoidSignature, SemiRingSignature,
    UniqueFactorizationMonoidSignature,
};
use algebraeon_sets::structure::*;
use std::fmt::Debug;

#[derive(Debug)]
pub enum FindFactorResult<Element> {
    Irreducible,
    Composite(Element, Element),
}

pub fn factorize_by_find_factor<
    Exponent: SemiRingSignature + OrdSignature,
    RS: UniqueFactorizationMonoidSignature<FactoredExponent = Exponent> + AdditiveMonoidEqSignature,
>(
    ring: &RS,
    elem: RS::Set,
    partial_factor: &impl Fn(RS::Set) -> FindFactorResult<RS::Set>,
) -> Factored<RS::Set, Exponent::Set> {
    debug_assert!(!ring.is_zero(&elem));
    if ring.is_unit(&elem) {
        ring.factorizations().new_unit_impl(elem)
    } else {
        debug_assert!(!ring.is_unit(&elem));
        match partial_factor(elem.clone()) {
            FindFactorResult::Composite(g, h) => ring.factorizations().mul(
                &factorize_by_find_factor(ring, g, partial_factor),
                &factorize_by_find_factor(ring, h, partial_factor),
            ),
            FindFactorResult::Irreducible => {
                //f is irreducible
                ring.factorizations().new_irreducible_impl(&elem)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::FactoringMonoidSignature;
    use algebraeon_nzq::{Integer, Natural};
    use algebraeon_sets::structure::MetaType;

    #[test]
    fn factorization_invariants() {
        let f = Integer::structure()
            .factorizations()
            .new_unit_and_powers_unchecked(
                Integer::from(-1),
                vec![
                    (Integer::from(2), Natural::from(2u8)),
                    (Integer::from(3), Natural::from(1u8)),
                ],
            );
        Integer::structure()
            .factorizations()
            .is_element(&f)
            .unwrap();

        let f = Integer::structure()
            .factorizations()
            .new_unit_and_powers_unchecked(Integer::from(1), vec![]);
        Integer::structure()
            .factorizations()
            .is_element(&f)
            .unwrap();

        let f = Integer::structure()
            .factorizations()
            .new_unit_and_powers_unchecked(
                Integer::from(-1),
                vec![
                    (Integer::from(2), Natural::from(2u8)),
                    (Integer::from(3), Natural::from(1u8)),
                    (Integer::from(5), Natural::from(0u8)),
                ],
            );
        assert!(
            Integer::structure()
                .factorizations()
                .is_element(&f)
                .is_err(),
            "can't have a power of zero"
        );

        let f = Integer::structure()
            .factorizations()
            .new_unit_and_powers_unchecked(
                Integer::from(3),
                vec![(Integer::from(2), Natural::from(2u8))],
            );
        assert!(
            Integer::structure()
                .factorizations()
                .is_element(&f)
                .is_err(),
            "unit should be a unit"
        );

        let f = Integer::structure()
            .factorizations()
            .new_unit_and_powers_unchecked(
                Integer::from(1),
                vec![
                    (Integer::from(0), Natural::from(1u8)),
                    (Integer::from(3), Natural::from(1u8)),
                ],
            );
        assert!(
            Integer::structure()
                .factorizations()
                .is_element(&f)
                .is_err(),
            "prime factors must not be zero"
        );

        let f = Integer::structure()
            .factorizations()
            .new_unit_and_powers_unchecked(
                Integer::from(-1),
                vec![
                    (Integer::from(4), Natural::from(1u8)),
                    (Integer::from(3), Natural::from(1u8)),
                ],
            );
        assert!(
            Integer::structure()
                .factorizations()
                .is_element(&f)
                .is_err(),
            "prime factors must be prime"
        );

        let f = Integer::structure()
            .factorizations()
            .new_unit_and_powers_unchecked(
                Integer::from(-1),
                vec![
                    (Integer::from(-2), Natural::from(2u8)),
                    (Integer::from(3), Natural::from(1u8)),
                ],
            );
        assert!(
            Integer::structure()
                .factorizations()
                .is_element(&f)
                .is_err(),
            "prime factors must be favoriate associate"
        );
    }

    #[test]
    fn test_count_divisors() {
        for a in 1..25 {
            println!("a = {}", a);
            let b = Integer::from(a);
            println!("b = {}", b);
            let fs = Integer::structure().factor(&b);
            println!("fs = {:?}", fs);
            assert_eq!(
                Integer::structure()
                    .factorizations()
                    .count_divisors(&fs)
                    .unwrap(),
                Natural::from(
                    Integer::structure()
                        .factorizations()
                        .divisors(&fs)
                        .unwrap()
                        .collect::<Vec<Integer>>()
                        .len()
                )
            );
        }
    }
}
