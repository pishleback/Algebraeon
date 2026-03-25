pub struct ModuloOddStructure {
    n: u64,
}

impl ModuloOddStructure {
    pub fn new(n: u64) -> Self {
        Self { n }
    }
}

#[cfg(test)]
mod tests {}
