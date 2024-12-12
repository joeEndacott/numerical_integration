fn main() {}

fn flat(_x: f64) -> f64 {
    // This function takes a real number, xå, as an input and returns 1
    1.0
}

fn linear(x: f64) -> f64 {
    // This function takes a real number, x, as an input and returns x
    x
}

fn quadratic(x: f64) -> f64 {
    // This function takes a real number, x, as an input and returns x^2
    x * x
}

fn cubic(x: f64) -> f64 {
    // This function takes a real number, x, as an input and returns x^3
    x * x * x
}

fn quartic(x: f64) -> f64 {
    // This function takes a real number, x, as an input and returns x^4
    x * x * x * x
}

fn gaussian(x: f64) -> f64 {
    // This function takes a real number, x, as an input, and returns exp(- x^2)
    f64::exp(-x * x)
}

fn riemann_sum<F>(lower_limit: f64, upper_limit: f64, num_steps: u32, func: F) -> f64
where
    F: Fn(f64) -> f64,
{
    // riemann_sum approximates the definite integral of some function, func, between lower_limit and upper_limit
    // The domain is split up into num_steps equal size steps, and func is approximated as being flat over the range of each step
    // The integral over each step is computed, and all of these integrals are summed to approximate the total integral

    // Inputs:
    // lower_limit -> lower limit of integration
    // upper_limit -> upper limit of integration
    // num_steps -> number of steps that the domain is divided up into (to get a better approximation, make num_steps bigger)
    // func -> function to be integrated

    // Output: approximate value of definite integral

    // Error handling - num_steps must be > 0, and upper_limit must be > lower limit
    if num_steps == 0 || lower_limit >= upper_limit {
        panic!("Invalid inputs: num_steps must be > 0 and lower_limit < upper_limit.");
    }

    let mut integral = 0.0;
    let step_size = (upper_limit - lower_limit) / (num_steps as f64);

    for n in 0..num_steps {
        // Adding the contribution to the integral from each step
        integral += func(lower_limit + (n as f64) * step_size) * step_size;
    }

    integral
}

fn fixed_step_simpson<F>(lower_limit: f64, upper_limit: f64, num_steps: u32, func: F) -> f64
where
    F: Fn(f64) -> f64,
{
    // fixed_step_simpson approximates the definite integral of some function, func, between lower_limit and upper_limit
    // The domain is split up into 2 * num_steps equal size steps, and func is approximated as being quadratic over the range of each step
    // The integral over each step is computed, and all of these integrals are summed to approximate the total integral

    // Inputs:
    // lower_limit -> lower limit of integration
    // upper_limit -> upper limit of integration
    // num_steps -> half the number of steps that the domain is divided up into (to get a better approximation, make num_steps bigger)
    // The Simpson algorithm requires an even number of steps, hence why the number of steps used is twice num_steps
    // func -> function to be integrated

    // Output: approximate value of definite integral

    // Error handling - num_steps must be > 0, and upper_limit must be > lower limit
    if num_steps == 0 || lower_limit >= upper_limit {
        panic!("Invalid inputs: num_steps must be > 0 and lower_limit < upper_limit.");
    }

    let mut integral = 0.0;
    let step_size = (upper_limit - lower_limit) / (2.0 * num_steps as f64);

    // Adding the contribution to the integral from the lower and upper limits
    integral += step_size * (func(lower_limit) + func(upper_limit)) / 3.0;

    for n in 1..(2 * num_steps) {
        // The coefficients for odd n and even n are different
        let coefficient = if n % 2 == 0 { 2.0 } else { 4.0 };

        // Adding the contribution to the integral from each step
        integral += coefficient * step_size * func(lower_limit + (n as f64) * step_size) / 3.0;
    }

    integral
}

fn variable_step_simpson<F>(
    lower_limit: f64,
    upper_limit: f64,
    step_size: f64,
    max_err: f64,
    func: F,
) -> f64
where
    F: Fn(f64) -> f64,
{
    // Error handling - init_step_size must be > 0, and upper_limit must be > lower limit
    if step_size <= 0.0 || lower_limit >= upper_limit {
        panic!("Invalid inputs! Ensure num_steps > 0, lower_limit < upper_limit.");
    }

    let mut integral = 0.0;
    let mut current_position = lower_limit;
    let mut current_segment_size = 2.0 * step_size;

    while current_position < upper_limit {
        let (segment_integral, segment_size) = adaptive_simpson_step(
            current_position,
            upper_limit,
            step_size,
            current_segment_size,
            max_err,
            &func,
        );
        integral += segment_integral;
        current_segment_size = segment_size;
        current_position += segment_size;
    }

    integral
}

fn adaptive_simpson_step<F>(
    start_position: f64,
    upper_limit: f64,
    step_size: f64,
    init_segment_size: f64,
    max_err: f64,
    func: F,
) -> (f64, f64)
where
    F: Fn(f64) -> f64,
{
    // Error handling for inputs to function
    if start_position >= upper_limit
        || step_size <= 0.0
        || init_segment_size <= 0.0
        || max_err <= 0.0
    {
        panic!("Invalid inputs! Ensure start_position < upper_limit, step_size > 0, init_segment_size > 0, max_err > 0.");
    }

    // Handles edge case if start_position is less than 2 steps away from upper_limit
    if start_position + 2.0 * step_size > upper_limit {
        let segment_integral =
            0.5 * (func(start_position) + func(upper_limit)) * (upper_limit - start_position); // Applies the trapezium rule to calculate the integral
        let segment_size = upper_limit - start_position;
        return (segment_integral, segment_size);
    }

    // Calculate the function values at the start_position, and the two adjacent points to the right.
    let f_0 = func(start_position);
    let f_1 = func(start_position + step_size);
    let f_2 = func(start_position + 2.0 * step_size);

    // Quadratic interpolation coefficients
    let a = (f_0 - 2.0 * f_1 + f_2) / 2.0 * step_size * step_size;
    let b = (-3.0 * f_0 + 4.0 * f_1 - f_2) / 2.0 * step_size;
    let c = f_0;

    let mut segment_size = init_segment_size.min(upper_limit - start_position); // If init_segment_size is too big, shrinks the integration segment

    // f_temp stores the output when func is evaluated at the rightmost point of the integration segment
    let mut f_temp = func(start_position + segment_size);

    // error stores the absolute value of the difference between the value of the true function and the quadratic interpolation
    let mut error = (f_temp - a * segment_size * segment_size + b * segment_size + c).abs();

    // If error is smaller/bigger than max_err, step_direction is positive/negative, and so the integration segment increases/decreases in size
    let step_direction = if error < max_err {
        step_size
    } else {
        -step_size
    };

    // Increases/decreases the size of the integration segment until the error is saturated, or the segment grows too small or too big
    while (error < max_err) == (step_direction > 0.0)
        && segment_size >= 2.0 * step_size
        && start_position + segment_size < upper_limit
    {
        segment_size += step_direction;
        f_temp = func(start_position + segment_size);
        error = (f_temp - a * segment_size * segment_size + b * segment_size + c).abs();
    }

    // Corrects for segment_size overshooting when step_direction > 0
    segment_size += (step_direction + step_size) / 2.0;

    // Approximates the function as a quadratic and integrates over the integration segment
    let segment_integral = (a * segment_size * segment_size * segment_size / 3.0)
        + (b * segment_size * segment_size / 2.0)
        + c * segment_size;

    (segment_integral, segment_size)
}
