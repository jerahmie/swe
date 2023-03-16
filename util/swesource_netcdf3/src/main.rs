use swesource_netcdf3::gaussian2d;
use netcdf3::{FileWriter, DataSet, Version};

fn main()
{
    let nx: usize = 100;
    let ny: usize = 100;
    let filename = "gaussian2d.nc";
    let g2 = gaussian2d(0.2, 5.0, 5.0, 0.1, 0.1, nx, ny);

    //let data_set: DataSet = {
    let mut data_set: DataSet = DataSet::new();
        
    let _res = data_set.add_fixed_dim("nx", nx).expect("Unable add nx to data set.");
    let _res = data_set.add_fixed_dim("ny", ny).expect("Unable to add ny to data set.");
    data_set.add_var_f32("h", &["nx", "ny"]).expect("Unable to add h to data_set");
    let mut file_writer: FileWriter = FileWriter::open(filename).unwrap();
    file_writer.set_def(&data_set, Version::Classic, 0).unwrap();
    file_writer.write_var_f32("h", g2.into_raw_vec().as_slice()).unwrap();
    file_writer.close().unwrap();

}

