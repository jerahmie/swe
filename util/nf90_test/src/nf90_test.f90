!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: nf90_test
!
! Test read data from NetCDF file produced by another program.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program nf90_test

  use netcdf

  implicit none 

  character(*), parameter :: ncfilename = "gaussian2d.nc"
  integer :: ncid, status
  integer :: ndims, nvars, nglobalatts,  unlimdimid, formatnum

  status = nf90_open(path=trim(ncfilename), mode=nf90_nowrite, ncid=ncid)
  if (status .ne. 0) print *, "[ERROR] nf90_open, status: ", status
  print *, "Openeded ", trim(ncfilename)

  status = nf90_inquire(ncid, ndims, nvars, nglobalatts, unlimdimid, formatnum)
  if (status .ne. 0) print *, "[ERROR] nf90_inquire, status: ", status
  print *, "Found number of dimensions: ", ndims
  print *, "Found number of variables: ", nvars
  print *, "Found number of global attributes: ", nglobalatts
  print *, "unimited Dimension ID: ", unlimdimid
  print *, "Found format nubmer: ", formatnum

  status = nf90_close(ncid)
  if (status .ne. 0) print *, "[ERROR] nf90_close, status: ", status

  print *, "nf90_test finished with no errors." 

end program nf90_test
