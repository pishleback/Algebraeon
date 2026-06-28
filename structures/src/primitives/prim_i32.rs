use crate::*;
use std::cmp::Ordering;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimitiveI32CanonicalStructure {}

impl MetaType for i32 {
    type Signature = PrimitiveI32CanonicalStructure;

    fn structure() -> Self::Signature {
        PrimitiveI32CanonicalStructure {}
    }
}

impl Signature for PrimitiveI32CanonicalStructure {}

impl SetSignature for PrimitiveI32CanonicalStructure {
    type Elem = i32;

    fn validate_element(&self, _x: &i32) -> Result<(), String> {
        Ok(())
    }
}

impl EqSignature for PrimitiveI32CanonicalStructure {
    fn equal(&self, a: &Self::Elem, b: &Self::Elem) -> bool {
        a == b
    }
}

impl PartialOrdSignature for PrimitiveI32CanonicalStructure {
    fn partial_cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Option<Ordering> {
        Some(self.cmp(a, b))
    }
}

impl OrdSignature for PrimitiveI32CanonicalStructure {
    fn cmp(&self, a: &Self::Elem, b: &Self::Elem) -> Ordering {
        a.cmp(b)
    }
}

impl CountableSetSignature for PrimitiveI32CanonicalStructure {
    fn generate_all_elements(&self) -> impl Iterator<Item = Self::Elem> {
        i32::MIN..=i32::MAX
    }
}

impl FiniteSetSignature for PrimitiveI32CanonicalStructure {
    fn size(&self) -> Natural {
        Natural::from(u32::MAX) + Natural::ONE
    }
}

impl EnumeratedOrdFiniteSetSignature for PrimitiveI32CanonicalStructure {
    fn list_all_elements_ordered(&self) -> Vec<Self::Elem> {
        self.list_all_elements()
    }

    fn element_to_enumeration(&self, elem: &i32) -> Natural {
        Natural::from(elem.wrapping_sub(i32::MIN) as u32)
    }

    fn enumeration_to_element(&self, num: &Natural) -> Option<Self::Elem> {
        if *num >= self.size() {
            None
        } else {
            Some(((u32::try_from(num).unwrap()) as i32).wrapping_add(i32::MIN))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MetaEnumeratedOrdFiniteSetSignature;

    #[test]
    fn test_enumeration() {
        assert_eq!(i32::MIN.element_to_enumeration(), Natural::from(0u32));
        assert_eq!((i32::MIN + 1).element_to_enumeration(), Natural::from(1u32));
        assert_eq!(
            (i32::MAX - 1).element_to_enumeration(),
            i32::structure().size() - Natural::from(2usize)
        );
        assert_eq!(
            (i32::MAX).element_to_enumeration(),
            i32::structure().size() - Natural::from(1usize)
        );

        assert_eq!(
            Some(i32::MIN),
            i32::structure().enumeration_to_element(&Natural::from(0u32))
        );
        assert_eq!(
            Some(i32::MIN + 1),
            i32::structure().enumeration_to_element(&Natural::from(1u32))
        );
        assert_eq!(
            Some(i32::MAX - 1),
            i32::structure()
                .enumeration_to_element(&(i32::structure().size() - Natural::from(2usize)))
        );
        assert_eq!(
            Some(i32::MAX),
            i32::structure()
                .enumeration_to_element(&(i32::structure().size() - Natural::from(1usize)))
        );
        assert_eq!(
            None,
            i32::structure().enumeration_to_element(&i32::structure().size())
        );
        assert_eq!(
            None,
            i32::structure()
                .enumeration_to_element(&(i32::structure().size() + Natural::from(1usize)))
        );
    }
}
