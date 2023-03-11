FC=gfortran
FC_PARALLEL=mpifort
FCFLAGS=-std=f2003
CC=gcc
CC_PARALLEL=mpicc
CBUILDFLAGS=-shared
CFLAGS=-Wall -fPIC
CSRCS=$(wildcard libswe/*.c)
CHEADERS=$(wildcard libswe/*.h)
FSRCS=equations.f90 observer.f90 sources.f90 swe.f90
OBJS=$(subst .f90,.o,$(FSRCS))
COBJS=$(subst .c,.o,$(CSRCS))
INCLUDES=-I/usr/include -I${HOME}/local/PnetCDF/include
MPI_INCLUDE_DIR=/usr/include/x86_64-linux-gnu/openmpi
LDFLAGS=-lnetcdff
LDFLAGS_PARALLEL=-lpnetcdf
LDLIBS=L/usr/lib/x86_64-linux-gnu
TARGET=swe
HALOS=halos
SWELIB=libswe.so

all: $(TARGET) $(SWELIB) $(HALOS)

$(HALOS): halos.o
	$(FC_PARALLEL) -o $(HALOS) halos.o $(LDFLAGS_PARALLEL)

halos.o: halos.f90
	$(FC_PARALLEL) -c halos.f90 $(INCLUDES)

swe: $(OBJS)
	$(FC) -o $(TARGET) $(OBJS) $(LDFLAGS)

$(OBJS): %.o : %.f90
	$(FC) $(FCFLAGS) -c $< $(INCLUDES)

libswe/%.o: libswe/%.c
	$(CC) $(CFLAGS) -c $< -o $@ $(INCLUDES) -Ilibswe

swe_f: test/swe_f.o
	$(FC) -o test/swe_f test/swe_f.o $(LDFLAGS)

test/swe_f.o: test/swe_f.f90
	$(FC) $(FCFLAGS) -c test/swe_f.f90 -o test/swe_f.o -I$(INCLUDES)
	
libswe.so: $(COBJS)
	$(CC) -shared -o libswe.so $(COBJS)

clean:
	rm -f *.o libswe/*.o $(TARGET) $(SWELIB) $(HALOS)

distclean: clean
	rm -r $(TARGET) 
