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
COBJS=$(subst .c,.o,$(CRSCS))
INCLUDE_DIR=/usr/include
MPI_INCLUDE_DIR=/usr/include/x86_64-linux-gnu/openmpi
LDFLAGS=-lnetcdff
LDLIBS=L/usr/lib/x86_64-linux-gnu
TARGET=swe
HALOS=halos
SWELIB=libswe.so

all: $(TARGET) $(SWELIB) $(HALOS)

$(HALOS): halos.o
	$(FC_PARALLEL) -o $(HALOS) halos.o

halos.o: halos.f90
	$(FC_PARALLEL) -c halos.f90 

swe: $(OBJS)
	$(FC) -o $(TARGET) $(OBJS) $(LDFLAGS)

$(OBJS): %.o : %.f90
	$(FC) $(FCFLAGS) -c $< -I$(INCLUDE_DIR)

$(SWELIB): equationsc.o
	$(CC) $(CBUILDFLAGS) -o $(SWELIB) equationsc.o

libswe/%.o: libswe/%.c
	$(CC) $(CFLAGS) -c $< -o $@ -I$(INCLUDE_DIR) -Ilibswe

swe_f: test/swe_f.o
	$(FC) -o test/swe_f test/swe_f.o $(LDFLAGS)

test/swe_f.o: test/swe_f.f90
	$(FC) $(FCFLAGS) -c test/swe_f.f90 -o test/swe_f.o -I$(INCLUDE_DIR)
	
libswe.so: $(COBJS)
	$(CC) -shared -o libswe.so $(COBJS)

clean:
	rm -f *.o libswe/*.o $(TARGET) libswe.so

distclean: clean
	rm -r $(TARGET) 
