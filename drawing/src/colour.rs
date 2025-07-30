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
    pub fn orange() -> Self {
        Self {
            rgb: [1.0, 0.5, 0.0],
        }
    }
    pub fn lime() -> Self {
        Self {
            rgb: [1.0, 0.5, 0.0],
        }
    }
    pub fn pink() -> Self {
        Self {
            rgb: [1.0, 0.0, 0.5],
        }
    }
    pub fn purple() -> Self {
        Self {
            rgb: [0.5, 0.0, 1.0],
        }
    }
    pub fn aquamarine() -> Self {
        Self {
            rgb: [0.0, 1.0, 0.5],
        }
    }
    pub fn azure() -> Self {
        Self {
            rgb: [0.0, 0.5, 1.0],
        }
    }
    pub fn grey() -> Self {
        Self {
            rgb: [0.5, 0.5, 0.5],
        }
    }
    pub fn lighten(self) -> Self {
        let t = 0.5;
        Self {
            rgb: [
                1.0 - t * (1.0 - self.rgb[0]),
                1.0 - t * (1.0 - self.rgb[1]),
                1.0 - t * (1.0 - self.rgb[2]),
            ],
        }
    }
    pub fn darken(self) -> Self {
        let t = 0.5;
        Self {
            rgb: [t * self.rgb[0], t * self.rgb[1], t * self.rgb[2]],
        }
    }
}
