use orthoclase_groups::examples::symmetric::Permutation;
use orthoclase_rings::number::modulo::Modulo;

pub fn modular_permutation<const N: usize>(
    f: impl Fn(Modulo<N>) -> Modulo<N>,
) -> Result<Permutation<N>, &'static str> {
    let mut perm = [0; N];
    for i in 0..N {
        perm[i] = f(i.into()).into()
    }
    Permutation::new(perm)
}

#[cfg(test)]
mod tests {
    use orthoclase_rings::structure::StructuredType;

    use super::*;

    #[test]
    pub fn test() {
        println!(
            "{:?}",
            modular_permutation::<23>(|i| (i.into_ring() * 7).into_set())
                .unwrap()
                .disjoint_cycles()
        );
    }
}
