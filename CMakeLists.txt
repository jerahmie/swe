cmake_minimum_required(VERSION 3.20)

enable_language(Fortran)
project(swe VERSION 0.1)

get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

add_subdirectory(src)
