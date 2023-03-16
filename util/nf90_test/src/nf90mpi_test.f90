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
  integer :: ierror, rank, nprocs, tag, status(MPI_STATUS_SIZE)
  integer :: ncid, ncstatus
  integer :: ndims, nvars, nglobalatts, unlimdimid, formatnum, nvardims
  integer :: i, j, vartype, num
  integer, dimension(1024) :: vardimids ! maximum number of dimensions
  character(len=32), dimension(1024) :: vardimnames ! maximum number of dimensions
  integer(kind=8) :: dimval, nx, ny, local_nx, local_ny, subgridn
  character(len=32), allocatable :: dimname(:), varname(:)
  character(len=32) :: dimnamei
  real(kind=4), dimension(:,:), allocatable :: h
  integer(kind=MPI_OFFSET_KIND) :: starts(1,2), counts(1,2)

  call MPI_Init(ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierror)

  subgridn = int(sqrt(real(nprocs)))  ! number of subgrids

  ncstatus = nf90mpi_open(mpi_comm=MPI_COMM_WORLD,  path=trim(ncfilename), mpi_info=0, &
                          omode=nf90_nowrite, ncid=ncid)

  if (ncstatus .ne. 0) then
    print *, "[ERROR] nf90_open, status: ", ncstatus
    call exit(1)
  end if

  print *, "Opened ", trim(ncfilename)

  ncstatus = nf90mpi_inquire(ncid, ndims, nvars, nglobalatts, unlimdimid, formatnum)
  call check_result(ncstatus, "nf90mpi_inquire", rank)
  print *, "Found number of dimensions: ", ndims
  print *, "Found number of variables: ", nvars
  print *, "Found number of global attributes: ", nglobalatts
  print *, "unimited Dimension ID: ", unlimdimid
  print *, "Found format number: ", formatnum

  ! get number and names of dimensions
  ncstatus = nfmpi_inq_ndims(ncid, ndims)
  call check_result(ncstatus, "nf_inq_ndims", rank)
  allocate(dimname(ndims))
  print *, "Found number of dimensions: ", ndims
  if (ndims .gt. 0) then
    do i=1,ndims
      ncstatus = nfmpi_inq_dimname(ncid, i, dimname(i)) 
      call check_result(ncstatus, "nfmpi_inq_dimname", rank)
      ncstatus = nfmpi_inq_dim(ncid, i, dimnamei, dimval)
      call check_result(ncstatus, "nfmpi_inq_dim", rank)
      print *, "found dimension: ", trim(dimname(i)), ": ", dimnamei, dimval
    end do
  endif

  ! get number and names of variables
  ncstatus = nfmpi_inq_nvars(ncid, nvars)
  allocate(varname(nvars))
  
  call check_result(ncstatus, "nfmpi_inq_nvars", rank)
  print *, "==== Variables ====" 
  do i=1,nvars
    ncstatus = nfmpi_inq_varname(ncid, i, varname(i))
    call check_result(ncstatus, "nfmpi_inq_varname", rank)
    ncstatus = nfmpi_inq_vartype(ncid, i, vartype)
    call check_result(ncstatus, "nfmpi_inq_vartype", rank)
    ncstatus = nfmpi_inq_varndims(ncid, i, nvardims)
    call check_result(ncstatus, "nfmpi_inq_varndims", rank)
    ncstatus = nfmpi_inq_vardimid(ncid, i, vardimids)
    call check_result(ncstatus, "nfmpi_inq_vardimid", rank)
    print *, "  Name: ", varname(i)
    print *, "  Type: ", vartype
    print *, "  N Dims: ", nvardims
    print *, "  Dimension IDs: ", vardimids(1:nvardims)
    print *, "  Dimension Names: "
    do j=1,nvardims
      print *, "    ", dimname(vardimids(j))
    enddo
  end do
  
  ! read data 

  ncstatus = nfmpi_inq_dim(ncid, i, dimnamei, ny)
  call check_result(ncstatus, "nfmpi_inq_dimid", rank)
  print*, "i: ", i, ny
  ncstatus = nfmpi_inq_dim(ncid, i, dimnamei, nx)
  call check_result(ncstatus, "nfmpi_inq_dimid", rank)
  print*, "i: ", i, nx
  
  local_nx = nx/subgridn
  local_ny = ny/subgridn

  print*, "rank: ", rank, nprocs
  starts(1,1) = local_nx*(mod(rank, subgridn)) + 1
  starts(1,2) = local_ny*(mod(rank, subgridn)) + 1
  counts(1,1) = local_nx
  counts(1,2) = local_ny
  print *, "starts: ", starts
  print *, "counts: ", counts
  print *, "local_nx, local_ny: ", local_nx, local_ny
  allocate(h(local_nx,local_ny))
  !print *, "rank, nprocs", rank, nprocs
  !ncstatus = nfmpi_get_var_real_all(ncid, 1, start, count, h)
  !ncstatus = nfmpi_get_var_real(ncid, 1, h)
  i=1
  num=2
  !print*, "h before", h
  ncstatus = nfmpi_get_vara_real_all(ncid, i, starts, counts, h)
  call check_result(ncstatus, "nfmpi_get_vara_float", rank)
  if ( rank .eq. 0 ) then
    print*, "h after", h
  endif

  ncstatus = nf90mpi_close(ncid)
  if (ncstatus .ne. 0) then
    print *, "[ERROR] nf90_close, status: ", ncstatus
  end if

  if (rank .eq. 0) then
    print *, "nf90_test finished with no errors." 
  end if

  call MPI_Finalize(ierror)
  
  contains

    ! check if pnetcdf function returns non zero error code.
    ! if error, display user message and exit program
    subroutine check_result(ncstatus, message, rank)
      
      implicit none
      integer :: ncstatus, rank
      character(*) :: message

      if (ncstatus /= 0) then
        write (6, *) "NetCDF procedure reported error: (", ncstatus, "), rank ", rank,": ", trim(message)
        call exit(ncstatus)
      end if
    end subroutine

end program nf90_test
