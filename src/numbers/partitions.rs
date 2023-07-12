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

impl<P: Fn(usize) -> bool + Copy> Iterator
    for PartitionIterator<P>
{
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


pub fn partitions_sized<P: Fn(usize) -> bool + Copy>(
    n: usize,
    x: usize,
    predicate: P,
) -> impl Iterator<Item = Vec<usize>> {
    PartitionIterator::new(n, x, predicate)
}

pub fn partitions_sized_zero<P: Fn(usize) -> bool + Copy>(
    n: usize,
    x: usize,
    predicate: P,
) -> impl Iterator<Item = Vec<usize>> {
    PartitionIterator::new(n + x, x, move |k| predicate(k - 1))
        .map(|part| part.into_iter().map(|k| k - 1).collect::<Vec<usize>>())
}

pub fn partitions<P: Fn(usize) -> bool + Copy>(
    n: usize,
    predicate: P,
) -> impl Iterator<Item = Vec<usize>> {
    (1..n + 1).flat_map(move |x| partitions_sized::<P>(n, x, predicate))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partitions() {
        let parts = partitions(12, |k| k % 2 == 1);
        println!("start");
        for part in parts {
            println!("{:?}", part);
        }
        println!("end");
        let parts = partitions(12, |k| true);
        println!("start");
        for part in parts {
            println!("{:?}", part);
        }
        println!("end");
    }
}
