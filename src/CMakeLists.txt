cmake_minimum_required(VERSION 3.17)


set (NETCDF_F90 "YES")
find_package(MPI REQUIRED)
find_package(NetCDF MODULE REQUIRED Fortran)
add_executable(swe
	equations.f90
	sources.f90
	observer.f90
	swe.f90)
target_link_libraries(swe NetCDF::NetCDF_Fortran)
add_executable(swe_mpi halos.f90)
target_link_libraries(swe_mpi NetCDF::NetCDF_Fortran)
