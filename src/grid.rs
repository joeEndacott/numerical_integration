/// # Grid
///
/// Description: this struct represents a grid of points in 1D. The
/// coordinate of each point corresponds to an element in the vector
/// grid_points.
///
/// Example use case: todo: add example use case.
///
#[derive(Debug, Clone)]
pub struct Grid {
    pub grid_points: Vec<f64>,
}

impl Grid {
    /// # New uniform grid
    ///
    /// Description: creates a uniform grid of num_points points between
    /// start_point and end_point inclusive.
    ///
    /// Example use case: todo: add example use case.
    ///
    pub fn new_uniform_grid(
        start_point: f64,
        end_point: f64,
        num_points: usize,
    ) -> Self {
        // Error handling for when start_point is greater than or equal to end_point or num_points is less than or equal to 1.
        if start_point >= end_point || num_points <= 1 {
            let grid_points: Vec<f64> = vec![start_point];
            return Grid { grid_points };
        }

        // step_size is the distance between adjacent grid points
        let step_size = (end_point - start_point) / (num_points as f64 - 1.0);

        // Creates a vector containing the grid points
        let grid_points: Vec<f64> = (0..num_points)
            .map(|i| start_point + (i as f64) * step_size)
            .collect();

        Grid { grid_points }
    }
}

/// # Grid function
///
/// Description: this struct represents a real-valued function of a real
/// variable sampled on a grid of 1D points. The value of the function at each
/// sampling point corresponds to an element in the vector function_values. The
/// sampling points are contained in the Grid grid.
///
#[derive(Debug)]
pub struct GridFunction {
    pub grid: Grid,
    pub function_values: Vec<f64>,
}

impl GridFunction {
    /// # New grid function
    ///
    /// Description: this function generates a grid function, given a
    /// real-valued function of a real variable func and a Grid grid.
    ///
    pub fn new_grid_function<F>(grid: Grid, func: F) -> Self
    where
        F: Fn(f64) -> f64,
    {
        // Creates a vector containing the value of func at each grid point.
        let function_values: Vec<f64> =
            grid.grid_points.iter().map(|&x| func(x)).collect();

        GridFunction {
            grid,
            function_values,
        }
    }
}
