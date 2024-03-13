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

    fn compose(a: Self, b: Self) -> Self {
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_whatever() {}
}
