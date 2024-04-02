//all size k subsets of n
pub fn subsets(n: usize, k: usize) -> Vec<Vec<usize>> {
    if k == 0 {
        vec![vec![]]
    } else if k > n {
        return vec![];
    } else {
        //a size k subsets of n is the same as
        //some element of n and a size k-1 subset of everything after that element
        let mut ss = vec![];
        for first in 0..n {
            for mut rest in subsets(n - first - 1, k - 1) {
                let mut s = vec![first];
                for r in rest {
                    s.push(r + first + 1);
                }
                ss.push(s);
            }
        }
        ss
    }
}

pub fn subsets_vec<T: Clone>(items: Vec<T>, k: usize) -> Vec<Vec<T>> {
    subsets(items.len(), k)
        .into_iter()
        .map(|subset| subset.into_iter().map(|idx| items[idx].clone()).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // pub fn run() {
    //     println!("{:?}", subsets(5, 3));
    // }
}
