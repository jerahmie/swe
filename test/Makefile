FC=gfortran
FCFLAGS='-std=f2003'
CC=gcc
CFLAGS=-fPIC

INCLUDE_DIR=/usr/include 
LDLIBS=/usr/lib/x86_64-linux-gnu
TARGET=swe_f

#all: ; $(info $$(COBJS) is [$(COBJS)])
all: $(TARGET)


$(TARGET): swe_f.o
	$(FC) -o $(TARGET) swe_f.o ../src/sources.o $(LDFLAGS) -L../src -lswe

swe_f.o: swe_f.f90
	$(FC) $(FCFLAGS) -c swe_f.f90 -I$(INCLUDE_DIR) -I../src
	
clean:
	rm -f *.o $(TARGET)

distclean: clean
	rm -r $(TARGET) 
