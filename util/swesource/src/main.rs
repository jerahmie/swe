use swesource::gaussian2d;
use netcdf;


fn main(){
    let nx = 101;
    let ny = 101;
    let g2 = gaussian2d(0.2, 5., 5., 0.1, 0.1, nx, ny);

    // create netcdf file
    let mut file = netcdf::create("gaussian2d.nc").expect("Could not create NetCDF file.");

    // declare dimensions
    file.add_dimension("nx", nx).expect("Could not add dimension 'nx'");
    file.add_dimension("ny", ny).expect("Could not add dimension 'ny'");

    // declare variables
    let mut var = file.add_variable::<f64>("h", &["nx", "ny"]).unwrap();
    var.put_values(&g2.into_raw_vec().as_slice(), ..).unwrap();
}
