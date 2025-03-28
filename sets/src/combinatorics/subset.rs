/// Returns all size k subsets of {0, 1, ..., n-1}.
/// ```
/// use algebraeon_sets::combinatorics::subsets;
/// assert_eq!(subsets(5, 3).collect::<Vec<_>>(), vec![
///     vec![0, 1, 2],
///     vec![0, 1, 3],
///     vec![0, 1, 4],
///     vec![0, 2, 3],
///     vec![0, 2, 4],
///     vec![0, 3, 4],
///     vec![1, 2, 3],
///     vec![1, 2, 4],
///     vec![1, 3, 4],
///     vec![2, 3, 4],
/// ]);
/// ```
pub fn subsets(n: usize, k: usize) -> impl Iterator<Item = Vec<usize>> {
    if k == 0 {
        vec![vec![]].into_iter()
    } else if k > n {
        vec![].into_iter()
    } else {
        //a size k subsets of n is the same as
        //some element of n and a size k-1 subset of everything after that element
        let mut ss = vec![];
        for first in 0..n {
            for rest in subsets(n - first - 1, k - 1) {
                let mut s = vec![first];
                for r in rest {
                    s.push(r + first + 1);
                }
                ss.push(s);
            }
        }
        ss.into_iter()
    }
}

/// Returns all size k subsets of items.
pub fn subsets_vec<'a, T: 'a + Clone>(
    items: Vec<T>,
    k: usize,
) -> impl 'a + Iterator<Item = Vec<T>> {
    subsets(items.len(), k)
        .map(move |subset| subset.into_iter().map(|idx| items[idx].clone()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    pub fn run() {
        println!("{:?}", subsets(5, 3).collect::<Vec<_>>());
        assert_eq!(subsets(5, 3).collect::<Vec<_>>().len(), 10)
    }
}
