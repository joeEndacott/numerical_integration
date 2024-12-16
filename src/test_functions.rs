/// This function takes a real number, x, as an input and returns x^3
pub fn cubic(x: f64) -> f64 {
    x * x * x
}

/// This function takes a real number, x, as an input and returns x^4
pub fn quartic(x: f64) -> f64 {
    x * x * x * x
}

/// This function takes a real number, x, as an input, and returns exp(- x^2)
pub fn gaussian(x: f64) -> f64 {
    f64::exp(-x * x)
}
