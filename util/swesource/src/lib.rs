use libm::exp;
use ndarray::Array2;

// Function: gaussian2d
// Returns a 2D Gaussian of the form exp(-1*i((x-x0)^2 + (y-y0)^2)/(2*sigma^2))
//
// Inputs: sigma
//          x0 - distribution center, x
//          y0 - distribution center, y 
//          dx - resolution, x
//          dy - resolution, y
//          nx - number of cells in x-direction
//          ny - number of cells in y-direction
// Returns: Array2 - 2D array of size nx x ny
//
pub fn gaussian2d(sigma: f64, x0: f64, y0: f64, dx: f64, dy: f64, nx: usize, ny: usize) -> Array2<f64>{
    let mut gaussian2d = Array2::<f64>::zeros((nx,ny));
    for ((i,j), elt) in gaussian2d.indexed_iter_mut(){
       *elt = exp(-1.0*((((i as f64)*dx -x0).powf(2.0))/(2.0*sigma.powf(2.0)) +
               (((j as f64)*dy -y0).powf(2.0))/(2.0*sigma.powf(2.0))));
    }
    gaussian2d
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gaussian2d_dims() {
        let g2d = gaussian2d(1.0, 0.0, 0.0, 0.1, 0.1, 10, 5);
        assert_eq!(g2d.raw_dim(), [10,5]);
    }
}
