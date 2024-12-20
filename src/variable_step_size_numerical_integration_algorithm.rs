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
