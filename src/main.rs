pub mod numerical_integration_algorithms;
pub mod test_functions;
pub mod variable_step_size_numerical_integration_algorithm;

use test_functions::{cubic, gaussian, quartic};
use variable_step_size_numerical_integration_algorithm::variable_step_simpson as numeric_integral;

fn main() {
    println!(
        "Integral of x^3 from 0 to 1: {}",
        numeric_integral(cubic, 0.0, 1.0, 0.1, 0.000001)
    );
    println!(
        "Integral of x^4 from 0 to 1: {}",
        numeric_integral(quartic, 0.0, 1.0, 0.1, 0.000001)
    );
    println!(
        "Integral of e^(-x^2) from -10 to 10: {}",
        numeric_integral(gaussian, -10.0, 10.0, 1.0, 0.0000001)
    );
}
