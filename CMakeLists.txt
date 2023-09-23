cmake_minimum_required(VERSION 3.18)

enable_language(Fortran)
project(swe VERSION 0.1)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

list( APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules") 

add_subdirectory(src)
add_subdirectory(util)
add_subdirectory(test)

