use algebraeon_nzq::Integer;

pub fn stirling_number1(n: usize, k: usize) -> Integer {
    if k > n {
        return Integer::ZERO;
    }
    if n == 0 {
        return if k == 0 { Integer::ONE } else { Integer::ZERO };
    }
    if k == 0 {
        return Integer::ZERO;
    }
    if k == n {
        return Integer::ONE;
    }

    let mut prev = vec![Integer::ZERO; k + 1];
    prev[0] = Integer::ONE;

    for i in 1..=n {
        let mut curr = vec![Integer::ZERO; k + 1];
        let limit = usize::min(i, k);
        for j in 1..=limit {
            curr[j] = prev[j - 1].clone() - &prev[j] * Integer::from(i - 1);
        }
        prev = curr;
    }

    prev[k].clone()
}

pub fn stirling_number2(n: usize, k: usize) -> Integer {
    if k > n {
        return Integer::ZERO;
    }
    if n == 0 {
        return if k == 0 { Integer::ONE } else { Integer::ZERO };
    }
    if k == 0 {
        return Integer::ZERO;
    }
    if k == n {
        return Integer::ONE;
    }

    let mut prev = vec![Integer::ZERO; k + 1];
    prev[0] = Integer::ONE;

    for i in 1..=n {
        let mut curr = vec![Integer::ZERO; k + 1];
        let limit = usize::min(i, k);
        for j in 1..=limit {
            curr[j] = prev[j - 1].clone() + &prev[j] * Integer::from(j);
        }
        prev = curr;
    }

    prev[k].clone()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stirling_first_kind_small_values() {
        assert_eq!(stirling_number1(0, 0), Integer::ONE);
        assert_eq!(stirling_number1(5, 5), Integer::ONE);
        assert_eq!(stirling_number1(5, 1), Integer::from(24));
        assert_eq!(stirling_number1(5, 2), Integer::from(-50));
        assert_eq!(stirling_number1(6, 3), Integer::from(-225));
    }

    #[test]
    fn stirling_second_kind_small_values() {
        assert_eq!(stirling_number2(0, 0), Integer::ONE);
        assert_eq!(stirling_number2(5, 5), Integer::ONE);
        assert_eq!(stirling_number2(5, 3), Integer::from(25));
        assert_eq!(stirling_number2(6, 2), Integer::from(31));
        assert_eq!(stirling_number2(8, 2), Integer::from(127));
    }

    #[test]
    fn stirling_invalid_indices_are_zero() {
        assert_eq!(stirling_number1(4, 5), Integer::ZERO);
        assert_eq!(stirling_number2(4, 5), Integer::ZERO);
        assert_eq!(stirling_number1(3, 0), Integer::ZERO);
        assert_eq!(stirling_number2(3, 0), Integer::ZERO);
    }
}
