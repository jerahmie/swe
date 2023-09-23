cmake_minimum_required(VERSION 3.17)


enable_language(Fortran)
set (NETCDF_F90 "YES")
find_package(NetCDF MODULE REQUIRED Fortran)
find_package(PnetCDF REQUIRED C Fortran)
find_package(MPI REQUIRED)
add_executable(swe
	equations.f90
	sources.f90
	observer.f90
	swe.f90)
target_link_libraries(swe NetCDF::NetCDF_Fortran)
add_executable(swe_mpi init_par.f90 halos.f90)
message( STATUS "MPI_FOUND: ${MPI_FOUND}, MPI_FORTRAN_FOUND: ${MPI_FORTRAN_FOUND}")
message( STATUS "MPI_INCLUDE_DIRS: ${MPI_Fortran_INCLUDE_DIRS}")
message( STATUS "MPI_LIBRARY_DIRS: ${MPI_Fortran_LIBRARIES}")
message( STATUS "PnetCDF_Fortran_FOUND: ${PnetCDF_Fortran_FOUND}")
message( STATUS "PnetCDF_Fortran_INCLUDE_DIRS: ${PnetCDF_Fortran_INCLUDE_DIRS}")
message( STATUS "PnetCDF_Fortran_LIBRARIES: ${PnetCDF_Fortran_LIBRARIES}")
include_directories(${MPI_Fortran_INCLUDE_DIRS})
include_directories(${PnetCDF_Fortran_INCLUDE_DIRS})
target_link_libraries(swe_mpi ${PnetCDF_Fortran_LIBRARIES}
	                      ${MPI_Fortran_LIBRARIES})


