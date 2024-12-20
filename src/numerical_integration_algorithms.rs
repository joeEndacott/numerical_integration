/// # Riemann sum numerical integration algorithm.
///
/// Description: this function approximates the definite integral of a
/// real-valued function of a real variable, func, between lower_limit and
/// upper_limit.
/// The domain is split up into num_steps equal size steps, and func is
/// approximated as being flat over the range of each step.
/// The integral over each step is computed, and all of these integrals are
/// summed to approximate the total integral.
///
/// Example use case: todo: add example use case.
///
pub fn fixed_step_riemann<F>(
    func: F,
    lower_limit: f64,
    upper_limit: f64,
    num_steps: u32,
) -> f64
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
/// Description: this function approximates the definite integral of a
/// real-valued function of a real variable, func, between lower_limit and
/// upper_limit.
/// The domain is split up into 2 * num_steps equal size steps, and func is
/// approximated as being quadratic over the range of each step.
/// The integral over each step is computed, and all of these integrals are
/// summed to approximate the total integral.
///
/// Example use case: todo: add example use case.
///
pub fn fixed_step_simpson<F>(
    func: F,
    lower_limit: f64,
    upper_limit: f64,
    num_steps: u32,
) -> f64
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
        integral += coefficient
            * step_size
            * func(lower_limit + (n as f64) * step_size)
            / 3.0;
    }

    integral
}

/// # Simpson's rule algorithm with variable step size
///
/// Description: this function approximates the definite integral of a
/// real-valued function of a real variable, func, between lower_limit and
/// upper_limit.
/// The domain is split up into "segments" of different sizes, and func is approximated as being quadratic over the range of each segment.
/// The integral over each segment is computed, and all of these integrals are
/// summed to approximate the total integral.
///
/// Example use case: todo: add example use case.
///
pub fn variable_step_simpson<F>(
    func: F,
    lower_limit: f64,
    upper_limit: f64,
    min_segment_size: f64,
    max_err: f64,
) -> f64
where
    F: Fn(f64) -> f64,
{
    // Error handling - min_step_size must be greater than 0, and upper_limit must be greater than lower limit.
    if min_segment_size <= 0.0 || lower_limit >= upper_limit {
        panic!(
            "Invalid inputs! Ensure num_steps > 0, lower_limit < upper_limit."
        );
    }

    let mut integral = 0.0;

    // current_position is the current position in the domain.
    let mut current_position = lower_limit;

    // float_err is used to compensate for the error when adding floating point numbers.
    // let float_err: f64 = f64::powf(10.0, -12.0);

    // This while loop advances through the domain until we reach upper_limit.
    while current_position < upper_limit {
        // Calculate the size of the next integration segment.
        let segment_size = get_segment_size(
            &func,
            current_position,
            min_segment_size,
            max_err,
            upper_limit,
        );

        // Calculate the integral of the next integration segment.
        let segment_integral =
            get_segment_integral(&func, current_position, segment_size);

        integral += segment_integral;
        current_position += segment_size;
    }

    integral
}

/// # Get segment size
///
/// Description: this function calculates how big the next integration segment
/// should be, and returns this value.
/// Initially, the size of the segment is set to init_segment_size.
/// The size of the segment is then doubled until the error in the integral is
/// greater than max_err, or the segment exceeds the edge of the domain.
/// The size of the largest segment which has an error less than max_err is
/// returned.
///
/// Example use case: todo: add example use case.
///
fn get_segment_size<F>(
    func: F,
    start_position: f64,
    init_segment_size: f64,
    max_err: f64,
    upper_limit: f64,
) -> f64
where
    F: Fn(f64) -> f64,
{
    // Error handling - start_position must be less than or equal to
    // upper_limit, init_segment_size must be greater than 0, and max_err must
    // be greater than 0.
    if start_position >= upper_limit
        || init_segment_size <= 0.0
        || max_err <= 0.0
    {
        panic!("Invalid inputs! Ensure start_position less than upper_limit, init_segment_size greater than 0 and max_err greater than 0.");
    }

    // Handles the edge case when upper_limit - start_position is less than or equal to init_segment_size.
    if upper_limit - start_position <= init_segment_size {
        return upper_limit - start_position;
    }

    let mut segment_size: f64 = init_segment_size;

    // Calculate the value of func at the start, middle and end of the segment.
    let f_0 = func(start_position);
    let f_1 = func(start_position + 0.5 * init_segment_size);
    let f_2 = func(start_position + init_segment_size);

    // a, b, c are quadratic interpolation coefficients.
    // Todo: move this step into a separate function.
    let a =
        (f_0 - 2.0 * f_1 + f_2) / (0.5 * init_segment_size * init_segment_size);
    let b = (-3.0 * f_0 + 4.0 * f_1 - f_2) / init_segment_size;
    let c = f_0;

    // f_temp stores the value of func at the end of the integration segment.
    let mut f_temp = func(start_position + segment_size);

    // error stores an estimate of the error in the integral.
    let mut error =
        ((f_temp - a * segment_size * segment_size - b * segment_size - c)
            * segment_size)
            .abs();

    // This while loop doubles segment_size until either the error saturates,
    // or the segment exceeds the edge of the domain.
    while error < max_err && start_position + segment_size <= upper_limit {
        segment_size *= 2.0;
        f_temp = func(start_position + segment_size);
        error =
            ((f_temp - a * segment_size * segment_size - b * segment_size - c)
                * segment_size)
                .abs();
    }

    // If the segment exceeds the edge of the domain and the error hasn't
    // saturated, segment_size is set so that the segment completely fills the
    // remaining portion of the domain. If the error has been saturated, we
    // return the second largest value of segment_size.
    if start_position + segment_size > upper_limit && error < max_err {
        upper_limit - start_position
    } else {
        segment_size *= 0.5;
        segment_size
    }
}

/// # Get segment integral
///
/// Description: uses Simpson's rule to integrate a real-valued function of a
/// real variable, func, over an integration segment.
/// First, func is approximated as a quadratic function using a quadratic
/// interpolation.
/// The integral of this quadratic function is then calculated.
/// Finally, the value of this integral is returned.
///
/// Example use case: todo: add example use case.
///
fn get_segment_integral<F>(
    func: F,
    start_position: f64,
    segment_size: f64,
) -> f64
where
    F: Fn(f64) -> f64,
{
    // Calculate the value of func at the start, middle and end of the segment.
    let f_0 = func(start_position);
    let f_1 = func(start_position + (segment_size / 2.0));
    let f_2 = func(start_position + segment_size);

    // Uses Simpson's rule to approximate the integral over the segment.
    let segment_integral = segment_size * (f_0 + 4.0 * f_1 + f_2) / 6.0;
    segment_integral
}
