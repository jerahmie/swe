# SWE: Shallow Water Equations

A simple model to gain experience with Fortran by creating a solver for the shallow water equation.

## Building with CMake
`$ mkdir buildtree && cd buildtree`
`$ cmake .. && make -j`nproc`

## Running Code
`$ ./src/swe`

## View Results
`$ ncview observer.nc`


## Specifying specific compilers
cmake -D CMAKE_C_COMPILER=gcc -D CMAKE_CXX_COMPILER=g++ -D CMAKE_Fortran_COMPILER=gfortran .
