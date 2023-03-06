FC=gfortran
FCFLAGS='-std=f2003'
CC=gcc
CFLAGS=-fPIC

SRCS=equations.f90 observer.f90 sources.f90 swe.f90
#CSRCS=equationsc.c observerc.c
CSRCS=$(wildcard libswe/*.c)
CHEADERS=$(subst .c,.h,$(CSRCS))
COBJS=$(subst .c,.o,$(CSRCS))
OBJS=$(subst .f90,.o,$(SRCS))
INCLUDE_DIR=/usr/include 
LDFLAGS=-lnetcdff
LDLIBS=/usr/lib/x86_64-linux-gnu
TARGET=swe

#all: ; $(info $$(COBJS) is [$(COBJS)])
all: 	libswe.so swe

swe: $(OBJS)
	$(FC) -o $(TARGET) $(OBJS) $(LDFLAGS)

swe.o: swe.f90
	$(FC) $(FCFLAGS) -c swe.f90 -I$(INCLUDE_DIR)

equations.o: equations.f90
	$(FC) $(FCFLAGS) -c equations.f90 -I$(INCLUDE_DIR)

observer.o: observer.f90
	$(FC) $(FCFLAGS) -c observer.f90 -I$(INCLUDE_DIR)

sources.o: sources.f90
	$(FC) $(FCFLAGS) -c sources.f90 -I$(INCLUDE_DIR)

libswe/%.o: libswe/%.c
	$(CC) $(CFLAGS) -c $< -o $@ -I$(INCLUDE_DIR) -Ilibswe
	
libswe.so: $(COBJS)
	$(CC) -shared -o libswe.so $(COBJS)

clean:
	rm -f *.o libswe/*.o $(TARGET) libswe.so

distclean: clean
	rm -r $(TARGET) 
