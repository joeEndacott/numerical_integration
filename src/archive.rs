/// # Variable step size Simpson's rule algorithm v1.
///
fn variable_step_simpson<F>(
    lower_limit: f64,
    upper_limit: f64,
    step_size: f64,
    max_segment_size: f64,
    max_err: f64,
    func: F,
) -> f64
where
    F: Fn(f64) -> f64,
{
    // variable_step_simpson approximates the definite integral of some function, func, between lower_limit and upper_limit
    // The domain is split up into segments of different sizes, and func is approximated as being parabolic over the range of each segment
    // The integral over each segment is computed, and all of these integrals are summed to approximate the total integral

    // Currently, this function is a bit clunky and you need to fiddle with the parameters to get the correct result
    // The function is also less accurate than the fixed step Simpson algorithm and the Riemann sum function
    // I am going to rewrite parts of the algorithm and change the way that segment_size is increased or decreased in the next iteration of the function

    // Inputs:
    // lower_limit -> lower limit of integration
    // upper_limit -> upper limit of integration
    // step_size -> size of each step. A segment is composed of two or more steps
    // max_segment_size -> maximum segment size
    // max_err -> related to the maximum allowed error for each segment integral
    // func -> function to be integrated

    // Output:
    // integral -> approximate value of definite integral

    // Error handling - init_step_size must be > 0, and upper_limit must be > lower limit
    if step_size <= 0.0 || lower_limit >= upper_limit {
        panic!(
            "Invalid inputs! Ensure num_steps > 0, lower_limit < upper_limit."
        );
    }

    let mut integral = 0.0;
    let mut current_position = lower_limit; // current_position is the current position in the domain of integration
    let mut current_segment_size = 2.0 * step_size; // A segment is composed of two or more steps
    let float_err: f64 = f64::powf(10.0, -12.0); // float_err is used to compensate for the error in addition of floating point numbers

    // let mut total_steps = 0; // For testing - remove when function is complete

    // While loop advances through the domain until we reach upper_limit
    while current_position < upper_limit {
        let (segment_integral, segment_size) = adaptive_simpson_step(
            current_position,
            upper_limit,
            step_size,
            current_segment_size,
            max_segment_size,
            max_err,
            float_err,
            &func,
        ); // adaptive_simpson_step calculates the size and integral of the next integration segment
        integral += segment_integral;
        current_segment_size = segment_size; // The current segment size is used as the initial guess for the size of the next segment
        current_position += segment_size;
        // total_steps += 1;
    }

    //println!("{total_steps}");
    integral
}

fn adaptive_simpson_step<F>(
    start_position: f64,
    upper_limit: f64,
    step_size: f64,
    init_segment_size: f64,
    max_segment_size: f64,
    max_err: f64,
    float_err: f64,
    func: F,
) -> (f64, f64)
where
    F: Fn(f64) -> f64,
{
    // adaptive_simpson_step finds the size and integral of the next integration segment

    // Inputs:
    // start_position -> starting coordinate of segment
    // upper_limit -> upper limit of integration
    // step_size -> size of each step. A segment is composed of two or more steps
    // init_segment_size -> initial guess for the size of the segment
    // max_segment_size -> maximum segment size
    // func -> function to be integrated

    // Outputs:
    // segment_integral -> approximate integral of function over the segment
    // segment_size -> size of segment

    // Error handling for inputs to function
    if start_position >= upper_limit
        || step_size <= 0.0
        || init_segment_size <= 0.0
        || max_err <= 0.0
    {
        panic!("Invalid inputs! Ensure start_position < upper_limit, step_size > 0, init_segment_size > 0, max_err > 0.");
    }

    // Handles edge case when start_position is less than 2 steps away from upper_limit
    if start_position + 2.0 * step_size > upper_limit - float_err {
        let segment_integral = 0.5
            * (func(start_position) + func(upper_limit))
            * (upper_limit - start_position); // Applies the trapezium rule to calculate the integral
        let segment_size = upper_limit - start_position;
        return (segment_integral, segment_size);
    }

    // Calculate the function values at the start_position, and the two adjacent points to the right.
    let f_0 = func(start_position);
    let f_1 = func(start_position + step_size);
    let f_2 = func(start_position + 2.0 * step_size);

    // Quadratic interpolation coefficients
    let a = (f_0 - 2.0 * f_1 + f_2) / (2.0 * step_size * step_size);
    let b = (-3.0 * f_0 + 4.0 * f_1 - f_2) / (2.0 * step_size);
    let c = f_0;

    let mut segment_size = init_segment_size.min(upper_limit - start_position); // If init_segment_size is too big, shrinks segment_size

    // f_temp stores the value of func at the rightmost point of the integration segment
    let mut f_temp = func(start_position + segment_size);

    // error stores the absolute value of the difference between the value of the true function and the quadratic interpolation at the rightmost point of the integration segment
    let mut error =
        (f_temp - a * segment_size * segment_size - b * segment_size - c).abs();

    // If error is smaller/bigger than max_err, step_direction is positive/negative, and so the integration segment increases/decreases in size
    let step_direction = if error < max_err {
        step_size
    } else {
        -step_size
    };

    // Increases/decreases the size of the integration segment until the error is saturated, or the segment grows too small or too big
    while (error < max_err) == (step_direction > 0.0)
        && segment_size >= 2.0 * step_size
        && segment_size <= max_segment_size
        && start_position + segment_size <= upper_limit + float_err
    {
        segment_size += step_direction;
        f_temp = func(start_position + segment_size);
        error =
            (f_temp - a * segment_size * segment_size - b * segment_size - c)
                .abs();
    }

    // Corrects for segment_size overshooting when step_direction > 0
    segment_size -= (step_direction + step_size) / 2.0;

    // Final check to ensure that segment_size is within the correct limits
    segment_size =
        segment_size.clamp(2.0 * step_size, upper_limit - start_position);

    // Approximates the function as a quadratic and integrates over the integration segment
    let segment_integral = (a * segment_size * segment_size * segment_size
        / 3.0)
        + (b * segment_size * segment_size / 2.0)
        + c * segment_size;

    (segment_integral, segment_size)
}

///
///
///
/// # Variable step size Simpson's rule algorithm v2.
///
/// # Simpson's rule algorithm with variable step size
///
/// Details:
/// Author -> Joe Endacott (https://github.com/joeEndacott).
/// Date of last modification -> 19/12/2024.
///
/// Description:
/// variable_step_simpson approximates the definite integral of some function,
/// func, between lower_limit and upper_limit.
/// The domain is split up into segments of different sizes, and func is approximated as being quadratic over the range of each segment.
/// The integral over each segment is computed, and all of these integrals are
/// summed to approximate the total integral.
///
/// Inputs:
/// func -> function to be integrated.
/// lower_limit -> lower limit of integration.
/// upper_limit -> upper limit of integration.
/// min_segment_size -> smallest allowed segment size.
/// max_err -> the approximate error in the integral over a segment.
///
/// Output: integral -> approximate value of definite integral.
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
        // Calculate the size and integral of the next integration segment.
        let segment_size = get_segment_size(
            &func,
            current_position,
            min_segment_size,
            max_err,
            upper_limit,
        );
        let segment_integral =
            get_segment_integral(&func, current_position, segment_size);
        integral += segment_integral;
        current_position += segment_size;
    }

    integral
}

/// # Get segment size
///
/// Details:
/// Author -> Joe Endacott (https://github.com/joeEndacott).
/// Date of last modification -> 19/12/2024.
///
/// Description:
/// Todo: add description.
///
/// Inputs:
/// func -> function to be integrated.
/// start_position -> starting coordinate of segment.
/// init_segment_size -> starting segment size.
/// max_err -> approximate error in the integral over the segment.
/// upper_limit -> upper limit of integration.
///
/// Outputs:
/// segment_size -> size of the segment.
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
    // Error handling for inputs to function
    if start_position >= upper_limit
        || init_segment_size <= 0.0
        || max_err <= 0.0
    {
        panic!("Invalid inputs! Ensure start_position less than upper_limit, init_segment_size greater than 0 and max_err greater than 0.");
    }

    // Handles edge case when upper_limit - start_position is less than or equal to init_segment_size
    if upper_limit - start_position <= init_segment_size {
        return upper_limit - start_position;
    }

    let mut segment_size: f64 = init_segment_size;

    // Calculate the function values at the start_position, and the two adjacent points to the right.
    let f_0 = func(start_position);
    let f_1 = func(start_position + 0.5 * init_segment_size);
    let f_2 = func(start_position + init_segment_size);

    // a, b, c are quadratic interpolation coefficients.
    // Todo: move this step into a separate function.
    let a =
        (f_0 - 2.0 * f_1 + f_2) / (0.5 * init_segment_size * init_segment_size);
    let b = (-3.0 * f_0 + 4.0 * f_1 - f_2) / init_segment_size;
    let c = f_0;

    // f_temp stores the value of func at the rightmost point of the integration segment.
    let mut f_temp = func(start_position + segment_size);

    // error stores an estimate of the error in the integral
    let mut error =
        ((f_temp - a * segment_size * segment_size - b * segment_size - c)
            * segment_size)
            .abs();

    while error < max_err && start_position + segment_size <= upper_limit {
        segment_size *= 2.0;
        f_temp = func(start_position + segment_size);
        error =
            ((f_temp - a * segment_size * segment_size - b * segment_size - c)
                * segment_size)
                .abs();
    }

    if start_position + segment_size > upper_limit && error < max_err {
        upper_limit - start_position
    } else {
        segment_size *= 0.5;
        segment_size
    }

    // Approximates the function as a quadratic and integrates over the integration segment
    //let segment_integral = (a * segment_size * segment_size * segment_size
    // 3.0)
    //+ (b * segment_size * segment_size / 2.0)
    //+ c * segment_size;
}

/// # Get segment integral
///
/// Details:
/// Author -> Joe Endacott (https://github.com/joeEndacott).
/// Date of last modification -> 20/12/2024.
///
/// Description:
/// Uses Simpson's rule to integrate func over an integration segment.
/// First, func is approximated using a quadratic interpolation.
/// The integral of this quadratic function is then calculated, and this integral is returned.
///
///
/// Inputs:
/// func -> function to be integrated.
/// start_position -> starting coordinate of segment.
/// segment_size -> segment size.
///
/// Output:
/// segment_integral -> approximate value of the integral of func over the segment.
///
fn get_segment_integral<F>(
    func: F,
    start_position: f64,
    segment_size: f64,
) -> f64
where
    F: Fn(f64) -> f64,
{
    let f_0 = func(start_position);
    let f_1 = func(start_position + (segment_size / 2.0));
    let f_2 = func(start_position + segment_size);

    let segment_integral = segment_size * (f_0 + 4.0 * f_1 + f_2) / 6.0;
    segment_integral
}

///
///
///
/// # Riemann sum algorithm v1.
///
/// # Riemann sum numerical integration algorithm.
///
/// Details:
/// Author -> Joe Endacott (https://github.com/joeEndacott).
/// Date of last modification -> 16/12/2024.
///
/// Description:
/// fixed_step_riemann approximates the definite integral of some function,
/// func, between lower_limit and upper_limit.
/// The domain is split up into num_steps equal size steps, and func is
/// approximated as being flat over the range of each step.
/// The integral over each step is computed, and all of these integrals are
/// summed to approximate the total integral.
///
/// Inputs:
/// func -> function to be integrated.
/// lower_limit -> lower limit of integration.
/// upper_limit -> upper limit of integration.
/// num_steps -> number of steps that the domain is divided up into (to get a
/// better approximation, make num_steps bigger).
///
/// Output: integral -> approximate value of definite integral.
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

///
///
///
/// # Fixed step size Simpson's rule algorithm v1.
///
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
/// func -> function to be integrated.
/// lower_limit -> lower limit of integration.
/// upper_limit -> upper limit of integration.
/// num_steps -> half the number of steps that the domain is divided up into
/// (to get a better approximation, make num_steps bigger).
/// The Simpson algorithm requires an even number of steps, hence why the number
/// of steps used is twice num_steps.
///
/// Output:
/// integral -> approximate value of definite integral.
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
