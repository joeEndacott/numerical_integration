pub mod numerical_integration_algorithms;
pub mod test_functions;

use numerical_integration_algorithms::fixed_step_simpson as numeric_integral;
use test_functions::{cubic, gaussian, quartic};

fn main() {
    println!(
        "Integral of x^3 from 0 to 1: {}",
        numeric_integral(0.0, 1.0, 100, cubic)
    );
    println!(
        "Integral of x^4 from 0 to 1: {}",
        numeric_integral(0.0, 1.0, 100, quartic)
    );
    println!(
        "Integral of e^(-x^2) from -10 to 10: {}",
        numeric_integral(-10.0, 10.0, 100, gaussian)
    );
}
