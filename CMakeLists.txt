cmake_minimum_required(VERSION 3.18)

enable_language(Fortran)
project(swe VERSION 0.1)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED True)
find_package(NetCDF 4.7 REQUIRED)

file(GLOB_RECURSE sources src/*.f90)
add_executable(prog ${sources})
