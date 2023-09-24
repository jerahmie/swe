program swe_mpi
    use mpi
    use pnetcdf
    use init_par
    use swe_mpi_help

    implicit none

    type(Neighbors) :: cell_neighbors
    integer :: num_args, err 
    integer :: rank, comm_size, tag, status(MPI_STATUS_SIZE)
    integer :: ncid, ncstatus
    integer :: i, j, ndims, dimids(2), varid
    integer :: nx, ny, local_nx, local_ny, nsubgrid
    integer, dimension(1024) :: vardimids
    character(255) :: ncfile_input, ncfile_output
    character(len=32), dimension(:), allocatable :: varname, dimnamei
    integer(kind=MPI_OFFSET_KIND), dimension(:), allocatable :: dimval
    logical :: res
    real(kind=4), dimension(:,:), allocatable :: h, dhx
    integer(kind=MPI_OFFSET_KIND) :: starts(2), counts(2)
  
    ! get filename from command line
    num_args = command_argument_count()
    if ( num_args .ne. 2 ) then
        call help()    
        call exit(1)
    end if

    call get_command_argument(1, ncfile_input)
    call get_command_argument(2, ncfile_output)

    
    inquire(file=trim(ncfile_input), exist=res)
    if (.not. res) then
        print *, "Could not find file: ", trim(ncfile_input)
        call exit(1)
    end if

    call MPI_Init(ncstatus)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ncstatus)
    call MPI_Comm_size(MPI_COMM_WORLD, comm_size, ncstatus)

    ! number of subgrids on side 
    nsubgrid = nlen_subgrids(comm_size)
    if (comm_size .ne. (nsubgrid*nsubgrid)) then
        print *, "Number of MPI processes not square."
        call exit(1)
    end if

    ! find nearest neighbors
    cell_neighbors = nearest_neighbors(rank, comm_size) 

    ! TODO: move load data calls to module
    ! open netcdf file for read
    ncstatus = nf90mpi_open(mpi_comm=MPI_COMM_WORLD, path=ncfile_input, &
                            omode=NF_NOWRITE, mpi_info=MPI_INFO_NULL, ncid=ncid)
    call ncdf_check(ncstatus, "nf90mpi_open", rank)

    ! get number, names, and values of dimensions
    ncstatus = nfmpi_inq_ndims(ncid, ndims)  ! number of dimensions
    call ncdf_check(ncstatus, "nfmpi_inq_ndims")
    
    allocate(dimval(ndims))
    allocate(dimnamei(ndims))
    
    do i=1,ndims
        ncstatus = nfmpi_inq_dim(ncid, i, dimnamei(i), dimval(i))
        call ncdf_check(ncstatus, "nfmpi_inq_dim")
        if (dimnamei(i) .eq. "nx") then
            nx = dimval(i)
        else if (dimnamei(i) .eq. "ny") then
            ny = dimval(i)
        endif
    end do

    
    local_nx = nx/nsubgrid
    local_ny = ny/nsubgrid

    ! read data subregion
    starts(1) = local_nx*mod(rank,nsubgrid)+1
    starts(2) = local_ny*(rank*local_ny/nx) + 1
    counts(1) = local_nx
    counts(2) = local_ny

    allocate(h(local_nx+2, local_ny+2))
    allocate(dhx(local_nx+2, local_ny+2))

    i=1 ! h is only variable in source netcdf
    ncstatus = nfmpi_get_vara_real_all(ncid, i , starts, counts, h)
    call ncdf_check(ncstatus, "nfmpi_get_vara_real_all", rank)
   
    ! close netcdf file
    ncstatus = nf90mpi_close(ncid)
    call ncdf_check(ncstatus, "nf90mpi_close", rank)

    ! execute solver
    call derivative(local_nx, local_ny, h, dhx)

    ! write results
    ! TODO: move to output module 
    ! create new netcdf file for read
    ncstatus = nf90mpi_create(mpi_comm=MPI_COMM_WORLD, path=trim(ncfile_output), &
                              cmode=NF_CLOBBER, mpi_info=MPI_INFO_NULL, ncid=ncid)
    call ncdf_check(ncstatus, "nf90mpi_open for write", rank)
    ncstatus = nfmpi_def_dim(ncid, "ny", int(ny, kind=MPI_OFFSET_KIND), dimids(2))
    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", rank)
    ncstatus = nfmpi_def_dim(ncid, "nx", int(nx, kind=MPI_OFFSET_KIND), dimids(1))
    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", rank)
    
    ncstatus = nfmpi_def_var(ncid, "var", NF90_FLOAT, 2, dimids, varid)
    call ncdf_check(ncstatus, "nf90mpi_def_var for write", rank)

    ncstatus = nf90mpi_enddef(ncid)
    call ncdf_check(ncstatus, "nf90mpi_enddef for write", rank)
    if (rank .eq. 0) then
        print*, "dimids: ", dimids, " varid: ", varid
    endif
    
    ncstatus = nfmpi_put_vara_real_all(ncid, i, starts, counts, h)
    call ncdf_check(ncstatus, "nfmpi_put_vara_real_all for write", rank)
    
    ! close netcdf file
    ncstatus = nf90mpi_close(ncid)
    call ncdf_check(ncstatus, "nf90mpi_close for write", rank)

    if (rank .eq. 0) then
        print *, "Program finished with ", comm_size, " threads."
    end if
    
    call MPI_Finalize(ncstatus)
    
    contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Function: nlen_subgrids
        ! Description: returns number of subgrids on side of in square grids
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        pure function nlen_subgrids(nprocs)
!
!            implicit none
!            
!            integer, intent(in) :: nprocs
!            integer :: nlen_subgrids
!
!            nlen_subgrids = int(sqrt(real(nprocs)))
!
!        end function nlen_subgrids



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Function: nearest_neighbors
        ! Description: Returns struct with nearest neighbors on a cyclic grid
        !
        ! Input:   i -- index of compute cell
        !          N -- total number of cells in square grid
        ! Returns: derived-type that contains the index of the cell neighbor
        !          to the right, left, above and below 
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        pure function nearest_neighbors(i,N)
!            
!            implicit none
!
!            integer, intent(in) :: i, N
!            integer :: ni
!
!            ! return value: nearest neighbors
!            type(Neighbors) :: nearest_neighbors
!
!            ni = int(sqrt(real(N,kind=8)))
!
!            ! note: (i/ni)*ni = integer equivalent of floor(i/ni)*ni
!            nearest_neighbors%right = mod(i+1,ni) + (i/ni)*ni 
!            nearest_neighbors%left = mod(i+ni-1,ni) + (i/ni)*ni
!            nearest_neighbors%below = mod(i+N-ni, N)
!            nearest_neighbors%above = mod(i+N+ni, N)
!
!        end function nearest_neighbors


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Function: ncdfcheck
        ! Description:
        !
        ! Input: ncstatus - ncstatus from PnetCDF function call 
        !        message - user informational message
        !        rank (optional) - mpi processor rank 
        !         
        ! Returns: Exits program if ncstatus is not zero
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine ncdf_check(ncstatus, message, rank)
       
            implicit none
            integer :: ncstatus
            integer, optional :: rank
            character(*) :: message
            
            if (ncstatus .ne. 0) then
                if (present(rank)) then
                    write (6, *) "PnetCDF procedure reported error: (" , ncstatus, "), rank ", rank, ": ", trim(message)   
                else
                    write (6, *) "PnetCDF procedure reported error: (", ncstatus, "): ", trim(message)
                end if
            end if
        
        end subroutine ncdf_check

    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Function: derivative
        ! Description:
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
               ! dx(i,j) = x(i+1,j) - x(i-1,j)
               dx(i,j) = x(i,j)
            end do
            end do
        end subroutine derivative

end program

