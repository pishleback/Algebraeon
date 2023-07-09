pub struct PartitionIterator {
    n: usize,
    x: usize,
    first: usize,
    min: usize,
    rest: Option<Box<PartitionIterator>>,
}

impl PartitionIterator {
    pub fn new(n: usize, x: usize) -> Self {
        Self {
            n,
            x,
            min: 1,
            first: 1,
            rest: None,
        }
    }
}

impl Iterator for PartitionIterator {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.n == 0 && self.x == 0 {
            if self.first == 0 {
                self.first += 1;
                Some(vec![])
            } else {
                None
            }
        } else if self.n == 0 || self.x == 0 {
            return None;
        } else if self.x == 1 {
            if self.first <= self.n {
                self.first = self.n + 1;
                Some(vec![self.n])
            } else {
                None
            }
        } else {
            if self.rest.is_none() {
                if self.first > self.n {
                    return None;
                } else {
                    self.rest = Some(Box::new(PartitionIterator {
                        n: self.n - self.first,
                        x: self.x - 1,
                        min: self.min + 1,
                        first: self.first,
                        rest: None,
                    }));
                }
            }
            match self.rest.as_mut().unwrap().as_mut().next().as_mut() {
                Some(rest_part) => {
                    let mut part = vec![self.first];
                    part.append(rest_part);
                    Some(part)
                }
                None => {
                    self.rest = None;
                    self.first += 1;
                    self.next()
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partitions() {
        let parts = PartitionIterator::new(10, 3);
        println!("start");
        for part in parts {
            println!("{:?}", part);
        }
        println!("end");
    }
}
