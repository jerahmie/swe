!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: nf90_test
!
! Test read data from NetCDF file produced by another program.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program nf90_test

  use mpi
  use pnetcdf

  implicit none 

  character(*), parameter :: ncfilename = "gaussian2d.nc"
  integer :: ierror, rank, size, tag, status(MPI_STATUS_SIZE)
  integer :: ncid, ncstatus
  integer :: ndims, nvars, nglobalatts,  unlimdimid, formatnum

  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierror)

  ncstatus = nf90mpi_open(mpi_comm=MPI_COMM_WORLD,  path=trim(ncfilename), mpi_info=0, &
                          omode=nf90_nowrite, ncid=ncid)
  if (ncstatus .ne. 0) then
    print *, "[ERROR] nf90_open, status: ", ncstatus
    call exit(1)
  end if

  print *, "Openeded ", trim(ncfilename)

  ncstatus = nf90mpi_inquire(ncid, ndims, nvars, nglobalatts, unlimdimid, formatnum)
  if (ncstatus .ne. 0) then
  print *, "[ERROR] nf90_inquire, status: ", ncstatus
  call exit(1)
  end if
  print *, "Found number of dimensions: ", ndims
  print *, "Found number of variables: ", nvars
  print *, "Found number of global attributes: ", nglobalatts
  print *, "unimited Dimension ID: ", unlimdimid
  print *, "Found format nubmer: ", formatnum

  ncstatus = nf90mpi_close(ncid)
  if (ncstatus .ne. 0) then
    print *, "[ERROR] nf90_close, status: ", ncstatus
  end if


  if (rank .eq. 0) then
    print *, "nf90_test finished with no errors." 
  end if

  call MPI_Finalize(ierror)

end program nf90_test
