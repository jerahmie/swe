Testing of Shallow Water Equations solver components.

This consists of a C++ main testing loop that targets c-library implementation of SWE functions
contained in libswe.


gcc -o swe_c swe_c.c -I../src/libswe -L../src -lswe
LD_LIBRARY_PATH=../src ./swe_c
