!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: nf90_test
!
! Test read data from NetCDF file produced by another program.
!
! heavily influenced by 
! https://github.com/Parallel-NetCDF/PnetCDF/blob/master/examples/F90/put_var.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check(err, message)
    use mpi
    use pnetcdf
    implicit none
    integer :: err
    character (len=*) :: message

    if (err .ne. NF90_NOERR) then
      write(6,*) trim(message), trim(nf90mpi_strerror(err))
      call MPI_Abort(MPI_COMM_WORLD, -1, err)
    end if
    
  end subroutine check


program nf90_test

  use mpi
  use pnetcdf

  implicit none 

  character(len=256) :: arg, ncfilename
  integer :: err, i
  integer :: ierr, rank, nprocs, tag, status(MPI_STATUS_SIZE)
  integer :: ncid, ncstatus
  integer :: ndims, nvars, nglobalatts,  unlimdimid, formatnum
  logical :: verbose

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

 ! take filename from command line 
 if (rank .eq. 0) then
  ncfilename = "testfile.nc"
  i = 0
  do
    call get_command_argument(i, arg)
    if (len_trim(arg) == 0) exit
      ncfilename = trim(arg)
    i = i + 1
  end do
 end if

 call MPI_Bcast(ncfilename, 256, MPI_CHARACTER, 0, MPI_COMM_WORLD, err)

  ! set parameters 
  ! nx
  ! ny

  err = nf90mpi_create(MPI_COMM_WORLD, ncfilename, NF90_NOCLOBBER, MPI_INFO_NULL, ncid)
  call check(err, "nf90mpi_create")

  err = nf90mpi_close(ncid)
  call check(err, "nf90mpi_close")

  if (rank .eq. 0) then
    print *, "nf90mpi_create finished with no errors." 
    print *, "total number of MPI processes ", nprocs
  end if


  call MPI_Finalize(ierr)

end program nf90_test
