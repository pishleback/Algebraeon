pub trait DivMod<T> {
    type DivOutput;
    type ModOutput;
    fn div_mod(self, other: T) -> (Self::DivOutput, Self::ModOutput);
}

pub trait AbsDiff<T> {
    type Output;
    fn abs_diff(self, rhs: T) -> Self::Output;
}
