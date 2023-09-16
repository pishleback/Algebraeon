use super::set::*;
use std::collections::HashSet;

#[derive(Clone, Debug)]
pub struct Function<DomainT: SetT, RangeT: SetT> {
    domain: DomainT,
    range: RangeT,
    func: Vec<usize>,
}

impl<DomainT: SetT, RangeT: SetT> Function<DomainT, RangeT> {
    pub fn check_state(&self) -> Result<(), &'static str> {
        if self.func.len() != self.domain.borrow().size() {
            return Err("function len and domain size dont match");
        }

        for x in self.domain.borrow().elems() {
            if !self.range.borrow().contains(self.func[x]) {
                return Err("function range does not contain image of element");
            }
        }

        Ok(())
    }

    pub fn apply(&self, x: usize) -> usize {
        self.func[x]
    }

    pub fn new_unchecked(domain: DomainT, range: RangeT, func: Vec<usize>) -> Self {
        Self {
            domain,
            range,
            func,
        }
    }

    pub fn domain(&self) -> DomainT {
        self.domain.clone()
    }

    pub fn range(&self) -> RangeT {
        self.range.clone()
    }

    pub fn is_injective(&self) -> bool {
        let mut img = HashSet::new();
        for x in self.domain.borrow().elems() {
            let y = self.apply(x);
            if img.contains(&y) {
                return false;
            } else {
                img.insert(y);
            }
        }
        return true;
    }

    pub fn is_surjective(&self) -> bool {
        let mut img = HashSet::new();
        for x in self.domain.borrow().elems() {
            img.insert(self.apply(x));
        }
        println!("img {:?} {}", img, self.range.borrow().size());
        return img.len() == self.range.borrow().size();
    }
}

#[derive(Clone, Debug)]
pub struct Bijection<DomainT: SetT, RangeT: SetT> {
    domain: DomainT,
    range: RangeT,
    func: Vec<usize>,
    inv: Vec<usize>,
}

impl<DomainT: SetT, RangeT: SetT> Bijection<DomainT, RangeT> {
    pub fn check_state(&self) -> Result<(), &'static str> {
        if self.domain.borrow().size() != self.range.borrow().size() {
            return Err("bijection sets have different sizes");
        }

        match self.func().check_state() {
            Err(msg) => {
                return Err(msg);
            }
            Ok(()) => {}
        }

        match self.inverse().check_state() {
            Err(msg) => {
                return Err(msg);
            }
            Ok(()) => {}
        }

        for x in self.domain.borrow().elems() {
            if x != self.inv[self.func[x]] {
                return Err("bijection inverse failed");
            }
        }

        for y in self.range.borrow().elems() {
            if y != self.func[self.inv[y]] {
                return Err("bijection inverse failed");
            }
        }

        Ok(())
    }

    pub fn func(&self) -> Function<DomainT, RangeT> {
        Function {
            domain: self.domain.clone(),
            range: self.range.clone(),
            func: self.func.clone(),
        }
    }

    pub fn inverse(&self) -> Function<RangeT, DomainT> {
        Function {
            domain: self.range.clone(),
            range: self.domain.clone(),
            func: self.inv.clone(),
        }
    }

    pub fn apply(&self, x: usize) -> usize {
        self.func[x]
    }

    pub fn apply_inverse(&self, y: usize) -> usize {
        self.inv[y]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_function_state() {
        let a = Set::new(5);
        let b = Set::new(7);

        let f = Function {
            domain: &a,
            range: &b,
            func: vec![0, 4, 6, 1, 2],
        };
        assert!(f.check_state().is_ok());

        let f = Function {
            domain: &a,
            range: &b,
            func: vec![0, 4, 7, 1, 2],
        };
        assert!(f.check_state().is_err());

        let f = Function {
            domain: &a,
            range: &b,
            func: vec![0, 4, 1, 2],
        };
        assert!(f.check_state().is_err());
    }

    #[test]
    pub fn test_bijection_state() {
        let a = Set::new(5);
        let b = Set::new(5);

        let f = Bijection {
            domain: &a,
            range: &b,
            func: vec![0, 2, 3, 4, 1],
            inv: vec![0, 4, 1, 2, 3],
        };
        assert!(f.check_state().is_ok());

        let f = Bijection {
            domain: &a,
            range: &b,
            func: vec![0, 1, 4, 2, 3],
            inv: vec![0, 1, 4, 2, 3],
        };
        assert!(f.check_state().is_err());

        let a = Set::new(3);
        let b = Set::new(5);
        let f = Bijection {
            domain: &a,
            range: &b,
            func: vec![0, 1, 2],
            inv: vec![0, 1, 2, 0, 1],
        };
        assert!(f.check_state().is_err());
    }
}
