#![allow(dead_code)]

use std::collections::BTreeSet;
use std::rc::Rc;

use super::function::*;
use super::set::*;

#[derive(Debug)]
pub struct Subset<'a> {
    set: &'a Set,
    elems: BTreeSet<usize>,
}

impl<'a> Subset<'a> {
    pub fn check_state(&self) -> Result<(), &'static str> {
        for x in &self.elems {
            if !self.set.contains(*x) {
                return Err("subset contains invalid set element");
            }
        }
        Ok(())
    }

    pub fn natural_inclusion(&self) -> Function<Rc<Set>, &'a Set> {
        Function::new_unchecked(
            Rc::new(Set::new(self.elems.len())),
            self.set,
            self.elems.clone().into_iter().collect(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn test_subset_state() {
        let set = Set::new(0);
        //empty subset of empty set
        let subset = Subset {
            set: &set,
            elems: vec![].into_iter().collect(),
        };
        assert!(subset.check_state().is_ok());
        //bad subset of empty set
        let subset = Subset {
            set: &set,
            elems: vec![0].into_iter().collect(),
        };
        assert!(subset.check_state().is_err());

        let set = Set::new(7);
        //typical
        let subset = Subset {
            set: &set,
            elems: vec![1, 2, 3].into_iter().collect(),
        };
        assert!(subset.check_state().is_ok());
        //empty
        let subset = Subset {
            set: &set,
            elems: vec![].into_iter().collect(),
        };
        assert!(subset.check_state().is_ok());
        //full
        let subset = Subset {
            set: &set,
            elems: vec![0, 1, 2, 3, 4, 5, 6].into_iter().collect(),
        };
        assert!(subset.check_state().is_ok());
        //partially invalid
        let subset = Subset {
            set: &set,
            elems: vec![5, 6, 7, 8].into_iter().collect(),
        };
        assert!(subset.check_state().is_err());
        //fully invalid
        let subset = Subset {
            set: &set,
            elems: vec![7, 8, 9].into_iter().collect(),
        };
        assert!(subset.check_state().is_err());
    }

    #[test]
    pub fn test_subset_natural_inclusion() {
        let set = Set::new(7);
        let subset = Subset {
            set: &set,
            elems: vec![2, 3, 4].into_iter().collect(),
        };
        let f = subset.natural_inclusion();
        f.check_state().unwrap();
        assert_eq!(f.domain().size(), 3);
        assert_eq!(f.range().size(), 7);
        assert!(std::ptr::eq(&set, f.range()));
        assert!(f.is_injective());
    }
}
