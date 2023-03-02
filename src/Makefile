FC=gfortran

SRCS=equations.f90 observer.f90 sources.f90 swe.f90
OBJS=$(subst .f90,.o,$(SRCS))
INCLUDE_DIR=/usr/include
LDFLAGS=-lnetcdff
LDLIBS=/usr/lib/x86_64-linux-gnu
TARGET=swe

all: swe

swe: $(OBJS)
	$(FC) -o $(TARGET) $(OBJS) $(LDFLAGS)

swe.o: swe.f90
	$(FC) -c swe.f90 -I$(INCLUDE_DIR)

equations.o: equations.f90
	$(FC) -c equations.f90 -I$(INCLUDE_DIR)

observer.o: observer.f90
	$(FC) -c observer.f90 -I$(INCLUDE_DIR)

sources.o: sources.f90
	$(FC) -c sources.f90 -I$(INCLUDE_DIR)


clean:
	rm -f $(OBJS) $(TARGET)

distclean: clean
	rm -r $(TARGET) 
