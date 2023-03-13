FC=gfortran
FC_PARALLEL=mpifort
FCFLAGS=-std=f2003
CC=gcc
CC_PARALLEL=mpicc
CBUILDFLAGS=-shared
CFLAGS=-Wall -fPIC

SRCS=swelib/equationsc.c
FSRCS=equations.f90 observer.f90 sources.f90 swe.f90
OBJS=$(subst .f90,.o,$(FSRCS))
INCLUDE_DIR=/usr/include
LDFLAGS=-lnetcdff
LDLIBS=/usr/lib/x86_64-linux-gnu
TARGET=swe
SWELIB=libswe.so

all: $(TARGET) $(SWELIB)

swe: $(OBJS)
	$(FC) -o $(TARGET) $(OBJS) $(LDFLAGS) -L$(NETCDF)/lib

$(OBJS): %.o : %.f90
	$(FC) $(FCFLAGS) -c $< -I$(INCLUDE_DIR) -I$(NETCDF)/include

$(SWELIB): equationsc.o
	$(CC) $(CBUILDFLAGS) -o $(SWELIB) equationsc.o

equationsc.o: swelib/equationsc.c
	$(CC) $(CFLAGS) -c swelib/equationsc.c -I$(INCLUDE_DIR)

clean:
	rm -f *.o $(TARGET) $(SWELIB)

distclean: clean
	rm -r $(TARGET) 
