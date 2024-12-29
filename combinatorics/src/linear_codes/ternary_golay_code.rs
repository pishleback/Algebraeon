use algebraeon_rings::number::modulo::Modulo;

pub type TernaryField = Modulo<3>;
pub const ZERO: TernaryField = TernaryField::new(0);
pub const PLUS: TernaryField = TernaryField::new(1);
pub const MINUS: TernaryField = TernaryField::new(2);

struct ExtendedTernaryGolayCode {
    basis: [[usize; 12]; 6],
}

#[cfg(test)]
mod tests {
    use algebraeon_groups::examples::symmetric::Permutation;

    #[test]
    fn test_something() {
        let a = Permutation::new([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]).unwrap();

        println!("hiya {:?}", a);
    }
}
