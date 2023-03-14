FC=gfortran
FCPARALLEL=mpif90
TARGET=nf90_test
TARGET_MPI_TEST=nf90mpi_test
TARGET_MPI_CREATE=nf90mpi_create

all: $(TARGET) $(TARGET_MPI_TEST) $(TARGET_MPI_CREATE)

$(TARGET): nf90_test.o
	$(FC) -o nf90_test nf90_test.o -L$(NETCDF)/lib -lnetcdff

nf90_test.o: src/nf90_test.f90
	$(FC) -c src/nf90_test.f90 -I$(NETCDF)/include

$(TARGET_MPI_TEST): nf90mpi_test.o
	$(FCPARALLEL) -o nf90mpi_test nf90mpi_test.o -L$(PNETCDF)/lib -lpnetcdf

nf90mpi_test.o: src/nf90mpi_test.f90
	$(FCPARALLEL) -c src/nf90mpi_test.f90 -I$(PNETCDF)/include

$(TARGET_MPI_CREATE): nf90mpi_create.o
	$(FCPARALLEL) -o nf90mpi_create nf90mpi_create.o -L$(PNETCDF)/lib -lpnetcdf

nf90mpi_create.o: src/nf90mpi_create.f90
	$(FCPARALLEL) -c src/nf90mpi_create.f90 -I$(PNETCDF)/include



clean: 
	rm -f *.o nf90mpi_test nf90_test
