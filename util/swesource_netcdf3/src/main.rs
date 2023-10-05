use std::env;
//use swesource_netcdf3::gaussian2d;
use swesource_netcdf3::uniform2d;
use netcdf3::{FileWriter, DataSet, Version};

fn main()
{
    let nx: usize = 200;
    let ny: usize = 200;
//    let filename = "gaussian2d.nc";
    let filename = "uniform.nc";
    let dx = 2.0;
    let dy = 2.0;
//    let args: Vec<String> = env::args().collect();
//    let sigma = args[1].parse::<f64>().unwrap();
//    let x0 = args[2].parse::<f64>().unwrap();
//    let y0 = args[3].parse::<f64>().unwrap();
//    let dx = args[4].parse::<f64>().unwrap();
//    let dy = args[5].parse::<f64>().unwrap();
//    let nx = args[6].parse::<usize>().unwrap();
//    let ny = args[7].parse::<usize>().unwrap();
//
//    let g2 = gaussian2d(sigma, x0*dx, y0*dy, dx, dy, nx, ny);
    let g2 = uniform2d(1.0, nx, ny);
    let mut data_set: DataSet = DataSet::new();
        
    let _res = data_set.add_fixed_dim("nx", nx).expect("Unable add nx to data set.");
    let _res = data_set.add_fixed_dim("ny", ny).expect("Unable to add ny to data set.");
    data_set.add_var_f32("h", &["nx", "ny"]).expect("Unable to add h to data_set");
    let mut file_writer: FileWriter = FileWriter::open(filename).unwrap();
    file_writer.set_def(&data_set, Version::Classic, 0).unwrap();
    file_writer.write_var_f32("h", g2.into_raw_vec().as_slice()).unwrap();
    file_writer.close().unwrap();
}

