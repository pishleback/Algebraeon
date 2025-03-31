struct PartitionIterator<P: Fn(usize) -> bool + Copy> {
    n: usize,
    x: usize,
    predicate: P,
    first: usize,
    min: usize,
    rest: Option<Box<PartitionIterator<P>>>,
}

impl<P: Fn(usize) -> bool + Copy> PartitionIterator<P> {
    pub fn new(n: usize, x: usize, predicate: P) -> Self {
        Self {
            n,
            x,
            predicate,
            first: 1,
            min: 1,
            rest: None,
        }
    }
}

impl<P: Fn(usize) -> bool + Copy> Iterator for PartitionIterator<P> {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.n < self.x {
            None //optimization when there are no partitions
        } else if self.n == 0 && self.x == 0 {
            //there is one partition of 0 into 0 pieces, namely []
            if self.first == 1 {
                self.first += 1;
                Some(vec![])
            } else {
                None
            }
        } else if self.n == 0 || self.x == 0 {
            //if n=0 or x=0 but not both, then there are no partitions
            return None;
        } else if self.x == 1 {
            //the only partition of n into 1 piece is [n]
            if self.first <= self.n {
                self.first = self.n + 1;
                match (self.predicate)(self.n) {
                    true => Some(vec![self.n]),
                    false => None,
                }
            } else {
                None
            }
        } else {
            if self.first > self.n {
                return None;
            }
            //looking for more than 2 pieces
            if self.rest.is_none() {
                self.rest = Some(Box::new(PartitionIterator {
                    n: self.n - self.first,
                    x: self.x - 1,
                    predicate: self.predicate,
                    min: self.min + 1,
                    first: self.first,
                    rest: None,
                }));
            }

            if !(self.predicate)(self.first) {
                self.first += 1;
                self.rest = None;
                return self.next();
            }

            match self.rest.as_mut().unwrap().as_mut().next().as_mut() {
                Some(rest_part) => {
                    //yield [self.first, ...rest_part...]
                    let mut part = vec![self.first];
                    part.append(rest_part);
                    Some(part)
                }
                None => {
                    //exhausted all partitions of the form [self.first, ...]
                    self.first += 1;
                    self.rest = None;
                    self.next()
                }
            }
        }
    }
}

/// Returns all partitions of n into exactly x parts such that each part in the partition satisfies the predicate and all parts are non-zero.
/// ```
/// use algebraeon_sets::combinatorics::predicated_partitions_sized_zero;
/// assert_eq!(predicated_partitions_sized_zero(6, 2, |k| k % 2 == 1).collect::<Vec<_>>(), vec![
///     vec![1, 5],
///     vec![3, 3],
/// ]);
/// ```
/// ```
/// use algebraeon_sets::combinatorics::predicated_partitions_sized_zero;
/// assert_eq!(predicated_partitions_sized_zero(6, 3, |k| k <= 3).collect::<Vec<_>>(), vec![
///     vec![0, 3, 3],
///     vec![1, 2, 3],
///     vec![2, 2, 2],
/// ]);
/// ```
pub fn predicated_partitions_sized<P: Fn(usize) -> bool + Copy>(
    n: usize,
    x: usize,
    predicate: P,
) -> impl Iterator<Item = Vec<usize>> {
    PartitionIterator::new(n, x, predicate)
}

/// Returns all partitions of n into exactly x parts such that each part in the partition satisfies the predicate where parts may be zero.
/// ```
/// use algebraeon_sets::combinatorics::predicated_partitions_sized_zero;
/// assert_eq!(predicated_partitions_sized_zero(6, 2, |k| k % 2 == 1).collect::<Vec<_>>(), vec![
///     vec![1, 5],
///     vec![3, 3],
/// ]);
/// ```
/// ```
/// use algebraeon_sets::combinatorics::predicated_partitions_sized_zero;
/// assert_eq!(predicated_partitions_sized_zero(6, 3, |k| k <= 3).collect::<Vec<_>>(), vec![
///     vec![0, 3, 3],
///     vec![1, 2, 3],
///     vec![2, 2, 2],
/// ]);
/// ```
pub fn predicated_partitions_sized_zero<P: Fn(usize) -> bool + Copy>(
    n: usize,
    x: usize,
    predicate: P,
) -> impl Iterator<Item = Vec<usize>> {
    PartitionIterator::new(n + x, x, move |k| predicate(k - 1))
        .map(|part| part.into_iter().map(|k| k - 1).collect::<Vec<usize>>())
}

/// Returns all partitions of n such that each part in the partition satisfies the predicate.
/// ```
/// use algebraeon_sets::combinatorics::predicated_partitions;
/// assert_eq!(predicated_partitions(6, |k| k % 2 == 1).collect::<Vec<_>>(), vec![
///     vec![1, 5],
///     vec![3, 3],
///     vec![1, 1, 1, 3],
///     vec![1, 1, 1, 1, 1, 1],
/// ]);
/// ```
/// ```
/// use algebraeon_sets::combinatorics::predicated_partitions;
/// assert_eq!(predicated_partitions(4, |k| k <= 2).collect::<Vec<_>>(), vec![
///     vec![2, 2],
///     vec![1, 1, 2],
///     vec![1, 1, 1, 1],
/// ]);
/// ```
pub fn predicated_partitions<P: Fn(usize) -> bool + Copy>(
    n: usize,
    predicate: P,
) -> impl Iterator<Item = Vec<usize>> {
    (1..n + 1).flat_map(move |x| predicated_partitions_sized::<P>(n, x, predicate))
}

/// Returns all partitions of n into exactly x parts where all parts are non-zero.
/// ```
/// use algebraeon_sets::combinatorics::partitions_sized;
/// assert_eq!(partitions_sized(8, 3).collect::<Vec<_>>(), vec![
///     vec![1, 1, 6],
///     vec![1, 2, 5],
///     vec![1, 3, 4],
///     vec![2, 2, 4],
///     vec![2, 3, 3],
/// ]);
/// ```
pub fn partitions_sized(n: usize, x: usize) -> impl Iterator<Item = Vec<usize>> {
    predicated_partitions_sized(n, x, |_| true)
}

/// Returns all partitions of n into exactly x parts where parts may be zero.
/// ```
/// use algebraeon_sets::combinatorics::partitions_sized_zero;
/// assert_eq!(partitions_sized_zero(4, 3).collect::<Vec<_>>(), vec![
///     vec![0, 0, 4],
///     vec![0, 1, 3],
///     vec![0, 2, 2],
///     vec![1, 1, 2],
/// ]);
/// ```
pub fn partitions_sized_zero(n: usize, x: usize) -> impl Iterator<Item = Vec<usize>> {
    predicated_partitions_sized_zero(n, x, |_| true)
}

/// Returns all partitions of n.
/// ```
/// use algebraeon_sets::combinatorics::partitions;
/// assert_eq!(partitions(4).collect::<Vec<_>>(), vec![
///     vec![4],
///     vec![1, 3],
///     vec![2, 2],
///     vec![1, 1, 2],
///     vec![1, 1, 1, 1]
/// ]);
/// ```
pub fn partitions(n: usize) -> impl Iterator<Item = Vec<usize>> {
    predicated_partitions(n, |_| true)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partitions() {
        let parts = predicated_partitions(12, |k| k % 2 == 1);
        println!("start");
        for part in parts {
            println!("{:?}", part);
        }
        println!("end");

        let parts = predicated_partitions(12, |_k| true);
        println!("start");
        for part in parts {
            println!("{:?}", part);
        }
        println!("end");

        let parts = partitions_sized_zero(5, 3);
        println!("start");
        for part in parts {
            println!("{:?}", part);
        }
        println!("end");

        let parts = partitions(6);
        println!("start");
        for part in parts {
            println!("{:?}", part);
        }
        println!("end");
    }
}
