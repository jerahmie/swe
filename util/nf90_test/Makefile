FC=gfortran
TARGET=nf90_test

all: $(TARGET)

$(TARGET): src/nf90_test.o
	$(FC) -o release/nf90_test src/nf90_test.o -lnetcdff

src/nf90_test.o: src/nf90_test.f90
	$(FC) -c src/nf90_test.f90 

clean: 
	rm -f src/*.o release/nf90_test
