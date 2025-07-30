#[derive(Debug, Clone)]
pub struct Colour {
    pub rgb: [f32; 3],
}

impl Colour {
    pub fn black() -> Self {
        Self {
            rgb: [0.0, 0.0, 0.0],
        }
    }
    pub fn white() -> Self {
        Self {
            rgb: [1.0, 1.0, 1.0],
        }
    }
    pub fn red() -> Self {
        Self {
            rgb: [1.0, 0.0, 0.0],
        }
    }
    pub fn green() -> Self {
        Self {
            rgb: [0.0, 1.0, 0.0],
        }
    }
    pub fn blue() -> Self {
        Self {
            rgb: [0.0, 0.0, 1.0],
        }
    }
    pub fn cyan() -> Self {
        Self {
            rgb: [0.0, 1.0, 1.0],
        }
    }
    pub fn magenta() -> Self {
        Self {
            rgb: [1.0, 0.0, 1.0],
        }
    }
    pub fn yellow() -> Self {
        Self {
            rgb: [1.0, 1.0, 0.0],
        }
    }
}
