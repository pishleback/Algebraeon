use orthoclase_groups::examples::symmetric::*;
use orthoclase_rings::{number::modulo::Modulo, ring_structure::cannonical::*};


pub type TernaryField = Modulo<3>;
pub const ZERO: TernaryField = TernaryField::new(0);
pub const PLUS: TernaryField = TernaryField::new(1);
pub const MINUS: TernaryField = TernaryField::new(2);

struct ExtendedTernaryGolayCode {
    basis: [[usize; 12]; 6],
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_something() {
        let a = Permutation::new([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]).unwrap();

        println!("hiya {:?}", a);
    }
}
