use crate::{Integer, Natural, traits::Abs};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    InvalidInput,
}

pub fn stirling_number1_signed(n: usize, k: usize) -> Result<Integer, Error> {
    if k > n {
        return Err(Error::InvalidInput);
    }

    if n == 0 {
        return if k == 0 {
            Ok(Integer::ONE)
        } else {
            Ok(Integer::ZERO)
        };
    }
    if k == 0 {
        return Ok(Integer::ZERO);
    }
    if k == n {
        return Ok(Integer::ONE);
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

    Ok(prev[k].clone())
}

pub fn stirling_number1_unsigned(n: usize, k: usize) -> Result<Natural, Error> {
    Ok(stirling_number1_signed(n, k)?.abs())
}

pub fn stirling_number2(n: usize, k: usize) -> Result<Natural, Error> {
    if k > n {
        return Err(Error::InvalidInput);
    }
    if n == 0 {
        return if k == 0 {
            Ok(Natural::ONE)
        } else {
            Ok(Natural::ZERO)
        };
    }
    if k == 0 {
        return Ok(Natural::ZERO);
    }
    if k == n {
        return Ok(Natural::ONE);
    }

    let mut prev = vec![Natural::ZERO; k + 1];
    prev[0] = Natural::ONE;

    for i in 1..=n {
        let mut curr = vec![Natural::ZERO; k + 1];
        let limit = usize::min(i, k);
        for j in 1..=limit {
            curr[j] = prev[j - 1].clone() + &prev[j] * Natural::from(j);
        }
        prev = curr;
    }

    Ok(prev[k].clone())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stirling_first_kind_small_values() {
        assert_eq!(stirling_number1_signed(0, 0), Ok(Integer::ONE));

        assert_eq!(stirling_number1_signed(1, 0), Ok(Integer::ZERO));
        assert_eq!(stirling_number1_signed(1, 1), Ok(Integer::ONE));

        assert_eq!(stirling_number1_signed(2, 0), Ok(Integer::ZERO));
        assert_eq!(stirling_number1_signed(2, 1), Ok(Integer::from(-1)));
        assert_eq!(stirling_number1_signed(2, 2), Ok(Integer::ONE));

        assert_eq!(stirling_number1_signed(3, 0), Ok(Integer::ZERO));
        assert_eq!(stirling_number1_signed(3, 1), Ok(Integer::from(2)));
        assert_eq!(stirling_number1_signed(3, 2), Ok(Integer::from(-3)));
        assert_eq!(stirling_number1_signed(3, 3), Ok(Integer::ONE));

        assert_eq!(stirling_number1_signed(5, 5), Ok(Integer::ONE));
        assert_eq!(stirling_number1_signed(5, 1), Ok(Integer::from(24)));
        assert_eq!(stirling_number1_signed(5, 2), Ok(Integer::from(-50)));
        assert_eq!(stirling_number1_signed(6, 3), Ok(Integer::from(-225)));
    }

    #[test]
    fn stirling_second_kind_small_values() {
        assert_eq!(stirling_number2(0, 0), Ok(Natural::ONE));

        assert_eq!(stirling_number2(1, 0), Ok(Natural::ZERO));
        assert_eq!(stirling_number2(1, 1), Ok(Natural::ONE));

        assert_eq!(stirling_number2(2, 0), Ok(Natural::ZERO));
        assert_eq!(stirling_number2(2, 1), Ok(Natural::ONE));
        assert_eq!(stirling_number2(2, 2), Ok(Natural::ONE));

        assert_eq!(stirling_number2(3, 0), Ok(Natural::ZERO));
        assert_eq!(stirling_number2(3, 1), Ok(Natural::ONE));
        assert_eq!(stirling_number2(3, 2), Ok(Natural::from(3u64)));
        assert_eq!(stirling_number2(3, 3), Ok(Natural::ONE));

        assert_eq!(stirling_number2(5, 5), Ok(Natural::ONE));
        assert_eq!(stirling_number2(5, 3), Ok(Natural::from(25u64)));
        assert_eq!(stirling_number2(6, 2), Ok(Natural::from(31u64)));
        assert_eq!(stirling_number2(8, 2), Ok(Natural::from(127u64)));
    }

    #[test]
    fn stirling_invalid_indices_are_zero() {
        assert_eq!(stirling_number1_signed(4, 5), Err(Error::InvalidInput));
        assert_eq!(stirling_number2(4, 5), Err(Error::InvalidInput));
    }
}
