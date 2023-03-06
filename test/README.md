gcc -o swe_c swe_c.c -I../src/libswe -L../src -lswe
LD_LIBRARY_PATH=../src ./swe_c
