use super::group::Group;

#[derive(Debug, Clone)]
struct UniversalPermutation {
    forward: Vec<usize>,
    backward: Vec<usize>,
}

impl PartialEq for UniversalPermutation {
    fn eq(&self, other: &Self) -> bool {
        self.forward == other.forward && self.backward == other.backward
    }
}

impl Eq for UniversalPermutation {}

impl Group for UniversalPermutation {
    fn identity() -> Self {
        todo!()
    }

    fn inverse(self) -> Self {
        todo!()
    }

    fn compose_refs(a: &Self, b: &Self) -> Self {
        todo!()
    }

    fn compose_mut(&mut self, other: &Self) {
        *self = Self::compose_refs(self, other);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_whatever() {}
}
