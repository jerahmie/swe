program swe_mpi
    use mpi
    use pnetcdf

    implicit none

    type Neighbors
        integer :: right, left, below, above
    end type Neighbors

    type(Neighbors) :: cell_neighbors

    integer :: rank, size, ierror, tag, status(MPI_STATUS_SIZE)
    integer :: ncid
    character(255) :: cwd, ncfilename
    character(*), parameter :: ncfile_rel = "/../util/swesource/gaussian2d.nc"
    logical :: res

    ! check for input file
    call getcwd(cwd)
    ncfilename=trim(cwd)//ncfile_rel
    inquire(file=trim(ncfilename), exist=res)
    if (.not. res) then
        print *, "Could not find file: ", trim(ncfilename)
        call exit(1)
    end if

    call MPI_Init(ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierror)

    ! find nearest neighbors
    cell_neighbors = nearest_neighbors(rank, size) 

    ! open netcdf file for read
    ierror = nf90mpi_open(mpi_comm=MPI_COMM_WORLD, path=ncfilename, omode=NF_NOWRITE, mpi_info=MPI_INFO_NULL, ncid=ncid)
    if (ierror .ne. 0) print *, "nf90mpi_open error: ", ierror, ". Rank: ", rank

    ! close netcdf file
    ierror = nf90mpi_close(ncid)
    if (ierror .ne. 0) print *, "nf90mpi_close: ", ierror, ". Rank: ", rank

    if (rank .eq. 0) then
        print *, "Program finished with ", size, " threads."
    end if
    
    call MPI_Finalize(ierror)
    
    contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function: nearest_neighbors
!
! Input:   i -- index of compute cell
!          N -- total number of cells in square grid
! Returns: derived-type that contains the index of the cell neighbor
!          to the right, left, above and below 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
pure function nearest_neighbors(i,N)
    
    implicit none

    integer, intent(in) :: i, N
    integer :: ni

    ! return value: nearest neighbors
    type(Neighbors) :: nearest_neighbors

    ni = int(sqrt(real(N,kind=8)))

    ! note: (i/ni)*ni = integer equivalent of floor(i/ni)*ni
    nearest_neighbors%right = mod(i+1,ni) + (i/ni)*ni 
    nearest_neighbors%left = mod(i+ni-1,ni) + (i/ni)*ni
    nearest_neighbors%below = mod(i+N-ni, N)
    nearest_neighbors%above = mod(i+N+ni, N)

            
end function nearest_neighbors
    
subroutine derivative(nx, ny, x, dx)

    implicit none 

    integer, intent(in) :: nx, ny            ! field dims

    ! Actually, this MPI task only "owns" x(1:nx,1:ny)
    real, dimension(0:nx+1,0:ny+1), intent(in) :: x     ! field
    real, dimension(0:nx+1,0:ny+1), intent(out) :: dx   ! x-deriv of field

    integer :: i, j


    ! In setup routine
    ! *) NPX - num tasks in x-direction (initially, sqrt(comm_size))
    ! *) NPY - num tasks in y-direction (initially, sqrt(comm_size))
    ! *) PX  - x-coord in decomp = (comm_rank % NPX)
    ! *) PY  - y-coord in decomp = (comm_rank / NPY)

    !
    ! Do "halo exchange" here to ensure that
    ! x(0,:) and x(nx+1,:) are updated from neighboring processors
    !
    ! 1) Every MPI task needs to know its location in the overall decomposition (PX, PY)
    ! 2) For the halo exchange, updated left column from processor PX-1
    !                           updated right column from processor PX+1

    ! call MPI_Irecv(..., x(0,:), ...)     ! recv from task to left
    ! call MPI_Irecv(..., x(nx+1,:), ...)  ! recv from task to right

    ! call MPI_Isend(..., x(nx,:), ...)    ! send to task to right
    ! call MPI_Isend(..., x(1,:), ...)     ! send to task to left

    do j = 1, ny
    do i = 1, nx
       dx(i,j) = x(i+1,j) - x(i-1,j)
    end do
    end do
end subroutine derivative

end program

