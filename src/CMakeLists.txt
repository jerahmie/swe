cmake_minimum_required(VERSION 3.20)

set (NETCDF_F90 "YES")
#find_package (netCDF REQUIRED)
find_path (NETCDF_INCLUDE_DIR netcdf.mod PATHS /usr/include /usr/include)
find_library(NETCDF_LIBRARIES NAMES netcdff HINTS NETCDFF_DIR)
mark_as_advanced(NETCDF_LIBRARIES_F90)

include_directories(${NETCDF_INCLUDE_DIR})

add_executable(swe sources.f90 observer.f90 equations.f90 swe.f90)
target_link_libraries(swe ${NETCDF_LIBRARIES})
