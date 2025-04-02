pub fn compositions_sized(n: usize, x: usize) -> impl Iterator<Item = Vec<usize>> {
    let c: Box<dyn Iterator<Item = Vec<usize>>> = if n == 0 && x == 0 {
        Box::new(vec![vec![]].into_iter())
    } else if n == 0 || x == 0 {
        Box::new(vec![].into_iter())
    } else {
        Box::new(super::subsets(n - 1, x - 1).into_iter().map(move |mut s| {
            s.push(n - 1);
            let mut c = vec![];
            c.push(s[0] + 1);
            for i in 1..x {
                c.push(s[i] - s[i - 1]);
            }
            c
        }))
    };
    c
}

pub fn compositions_sized_zero(n: usize, x: usize) -> impl Iterator<Item = Vec<usize>> {
    compositions_sized(n + x, x).map(|c| c.into_iter().map(|i| i - 1).collect())
}

pub fn compositions(n: usize) -> impl Iterator<Item = Vec<usize>> {
    (0..(n + 1))
        .map(move |x| compositions_sized(n, x))
        .flatten()
}

#[cfg(test)]
mod tests {
    use crate::combinatorics::compositions_sized;

    #[test]
    fn number_composition_count() {
        assert_eq!(compositions_sized(0, 0).collect::<Vec<_>>().len(), 1);
        assert_eq!(compositions_sized(0, 1).collect::<Vec<_>>().len(), 0);
        assert_eq!(compositions_sized(0, 2).collect::<Vec<_>>().len(), 0);
        assert_eq!(compositions_sized(1, 0).collect::<Vec<_>>().len(), 0);
        assert_eq!(compositions_sized(1, 1).collect::<Vec<_>>().len(), 1);
        assert_eq!(compositions_sized(1, 2).collect::<Vec<_>>().len(), 0);
        assert_eq!(compositions_sized(2, 0).collect::<Vec<_>>().len(), 0);
        assert_eq!(compositions_sized(2, 1).collect::<Vec<_>>().len(), 1);
        assert_eq!(compositions_sized(2, 2).collect::<Vec<_>>().len(), 1);
        assert_eq!(compositions_sized(7, 3).collect::<Vec<_>>().len(), 15);
    }
}
