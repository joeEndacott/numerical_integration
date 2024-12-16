/// # Riemann sum numerical integration algorithm.
///
/// Details:
/// Author -> Joe Endacott (https://github.com/joeEndacott).
/// Date of last modification -> 16/12/2024.
///
/// Description:
/// fixed_step_riemann approximates the definite integral of some function,
/// func,
/// between lower_limit and upper_limit.
/// The domain is split up into num_steps equal size steps, and func is
/// approximated as being flat over the range of each step.
/// The integral over each step is computed, and all of these integrals are
/// summed to approximate the total integral.
///
/// Inputs:
/// lower_limit -> lower limit of integration.
/// upper_limit -> upper limit of integration.
/// num_steps -> number of steps that the domain is divided up into (to get a
/// better approximation, make num_steps bigger).
/// func -> function to be integrated.
///
/// Output: integral -> approximate value of definite integral.
///
pub fn fixed_step_riemann<F>(lower_limit: f64, upper_limit: f64, num_steps: u32, func: F) -> f64
where
    F: Fn(f64) -> f64,
{
    // Error handling - num_steps must be greater than 0, and upper_limit must
    // be greater than lower limit.
    if num_steps == 0 || lower_limit >= upper_limit {
        panic!(
            "Invalid inputs: num_steps must be > 0 and lower_limit 
        < upper_limit."
        );
    }

    let mut integral = 0.0;
    let step_size = (upper_limit - lower_limit) / (num_steps as f64);

    for n in 0..num_steps {
        // Adding the contribution to the integral from each step.
        integral += func(lower_limit + (n as f64) * step_size) * step_size;
    }

    integral
}

/// # Simpson's rule numerical integration algorithm.
///
/// Details:
/// Author -> Joe Endacott (https://github.com/joeEndacott).
/// Date of last modification -> 16/12/2024.
///
/// Description:
/// fixed_step_simpson approximates the definite integral of some function,
/// func, between lower_limit and upper_limit.
/// The domain is split up into 2 * num_steps equal size steps, and func is
/// approximated as being quadratic over the range of each step.
/// The integral over each step is computed, and all of these integrals are
/// summed to approximate the total integral.
///
/// Inputs:
/// lower_limit -> lower limit of integration.
/// upper_limit -> upper limit of integration.
/// num_steps -> half the number of steps that the domain is divided up into
/// (to get a better approximation, make num_steps bigger).
/// The Simpson algorithm requires an even number of steps, hence why the number
/// of steps used is twice num_steps.
/// func -> function to be integrated.
///
/// Output:
/// integral -> approximate value of definite integral.
///
pub fn fixed_step_simpson<F>(lower_limit: f64, upper_limit: f64, num_steps: u32, func: F) -> f64
where
    F: Fn(f64) -> f64,
{
    // Error handling - num_steps must be greater than 0 and upper_limit must be
    // greater than lower limit.
    if num_steps == 0 || lower_limit >= upper_limit {
        panic!(
            "Invalid inputs: num_steps must be > 0 and 
        lower_limit < upper_limit."
        );
    }

    let mut integral = 0.0;
    let step_size = (upper_limit - lower_limit) / (2.0 * num_steps as f64);

    // Adding the contribution to the integral from the lower and upper limits.
    integral += step_size * (func(lower_limit) + func(upper_limit)) / 3.0;

    for n in 1..(2 * num_steps) {
        // The coefficients for odd n and even n are different.
        let coefficient = if n % 2 == 0 { 2.0 } else { 4.0 };

        // Adding the contribution to the integral from each step.
        integral += coefficient * step_size * func(lower_limit + (n as f64) * step_size) / 3.0;
    }

    integral
}
