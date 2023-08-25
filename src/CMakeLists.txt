cmake_minimum_required(VERSION 3.20)

set (NETCDF_F90 "YES")
find_package (netCDF REQUIRED)
find_path (NETCDFF_INCLUDES_F90 netcdf.mod HINTS NETCDFF_DIR NETCDFF_DIR)
find_library(CETCDFF_LIBRARIES_F90 NAMES netcdff HINTS NETCDFF_DIR)
mark_as_advanced(NETCDF_LIBRARIES_F90)

#include_directories(${NETCDF_INCLUDES})

add_executable(swe sources.f90 observer.f90 equations.f90 swe.f90)
target_link_libraries(swe ${NETCDF_INCLUDES} ${NETCDF_LIBRARIES})
