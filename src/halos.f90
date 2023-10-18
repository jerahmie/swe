program swe_mpi
    use mpi
    use pnetcdf
    use init_par
    use swe_mpi_help
    use equations_par
    use observer_par

    implicit none

    type(Neighbors) :: cell_neighbors
    integer :: num_args
    integer :: comm_rank, comm_size, ncid, ncstatus
    integer :: i, j, it,  ndims, dimids(2), varid, nxid, nyid, timedimid
    integer :: local_nx, local_ny, nsubgrid
    character(255) :: ncfile_input, ncfile_output
    character(len=32), dimension(:), allocatable :: dimnamei ! varname
    integer(kind=MPI_OFFSET_KIND),dimension(:), allocatable :: dimval
    logical :: res
    real(kind=4), allocatable, dimension(:,:) :: h, u, v
    integer(kind=MPI_OFFSET_KIND) :: starts(2), counts(2)
 
    ! Simulation variables
    NX = 200
    NY = 200
    NT = 1
    !dt = 0.1

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
    call MPI_Comm_rank(MPI_COMM_WORLD, comm_rank, ncstatus)
    call MPI_Comm_size(MPI_COMM_WORLD, comm_size, ncstatus)

    ! number of subgrids on side 
    nsubgrid = nlen_subgrids(comm_size)
    if (comm_size .ne. (nsubgrid*nsubgrid)) then
        print *, "Number of MPI processes not square."
        call exit(1)
    end if

    ! find nearest neighbors
    cell_neighbors = nearest_neighbors(comm_rank, comm_size) 
! TODO: move load data calls to module
    ! open netcdf file for read
    ncstatus = nf90mpi_open(mpi_comm=MPI_COMM_WORLD, path=ncfile_input, &
                            omode=NF_NOWRITE, mpi_info=MPI_INFO_NULL, ncid=ncid)
    call ncdf_check(ncstatus, "nf90mpi_open", comm_rank)

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
    starts(1) = local_nx*mod(comm_rank,nsubgrid) + 1
    starts(2) = local_ny*(comm_rank*local_ny/nx) + 1
    counts(1) = local_nx
    counts(2) = local_ny

    allocate(h(0:local_nx+1, 0:local_ny+1))
    allocate(u(0:local_nx+1, 0:local_ny+1))
    allocate(v(0:local_nx+1, 0:local_ny+1))
    h = 0.0
    u = 0.0
    v = 0.0

    i=1 ! h is only variable in source netcdf
    ncstatus = nfmpi_get_vara_real_all(ncid, i , starts, counts, h(1:local_nx,1:local_ny))
    call ncdf_check(ncstatus, "nfmpi_get_vara_real_all", comm_rank)
   
    ! close netcdf file
    ncstatus = nf90mpi_close(ncid)
    call ncdf_check(ncstatus, "nf90mpi_close", comm_rank)


    ! write results
    ! TODO: move to output module 
    ! create new netcdf file for read
!    ncstatus = nf90mpi_create(mpi_comm=MPI_COMM_WORLD, path=trim(ncfile_output), &
!                              cmode=NF_CLOBBER, mpi_info=MPI_INFO_NULL, ncid=ncid)
!    call ncdf_check(ncstatus, "nf90mpi_open for write", comm_rank)
!    ncstatus = nfmpi_def_dim(ncid, "ny", int(ny, kind=MPI_OFFSET_KIND), nyid)
!    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", comm_rank)
!    ncstatus = nfmpi_def_dim(ncid, "nx", int(nx, kind=MPI_OFFSET_KIND), nxid)
!    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", comm_rank)
!    ncstatus = nfmpi_def_dim(ncid, "time", nf90mpi_unlimited, timedimid)
!    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", comm_rank)
!    
!    ncstatus = nfmpi_def_var(ncid, "var", NF90_FLOAT, 2, (/nxid, nyid, timedimid/), varid)
!    call ncdf_check(ncstatus, "nf90mpi_def_var for write", comm_rank)
!
!    ncstatus = nf90mpi_enddef(ncid)
!    call ncdf_check(ncstatus, "nf90mpi_enddef for write", comm_rank)
!    if (comm_rank .eq. 0) then
!        print*, "dimids: ", nxid, nyid, " varid: ", varid
!    endif
    observer_par_init(ncfile_output, )
    ! execute solver
    do it=1,NT
        u = u + dt*derivative(h, u, v, dx, dy, du)
        v = v + dt*derivative(h, u, v, dx, dy, dv)
        h = h + dt*derivative(h, u, v, dx, dy, dh)
        ncstatus = nfmpi_put_vara_real_all(ncid, i, starts, counts, h(1:local_nx,1:local_ny))
        call ncdf_check(ncstatus, "nfmpi_put_vara_real_all for write", comm_rank)
    end do
    
    ! close netcdf file
    ncstatus = nf90mpi_close(ncid)
    call ncdf_check(ncstatus, "nf90mpi_close for write", comm_rank)

    if (comm_rank .eq. 0) then
        print *, "Program finished with ", comm_size, " threads."
    end if
    
    call MPI_Finalize(ncstatus)
    deallocate(h)
    deallocate(u)
    deallocate(v)
    
    contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Function: mpi_check
        ! Description:
        !
        ! Input: ierr - MPI error code
        !        message - user informational message
        !        comm_rank (optional) - mpi processor rank 
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine check_mpi_status(ierr, message, comm_rank)

            implicit none

            integer :: ierr
            integer, optional :: comm_rank
            character(*) :: message

            if (ierr.ne. 0) then
                if (present(comm_rank)) then
                    write (6, *) "MPI procedure reported error: (" , &
                        ierr, "), comm_rank ", comm_rank, ": ", trim(message)   
                else
                    write (6, *) "MPI procedure reported error: (", &
                        ierr, "): ", trim(message)
                end if
            end if

        end subroutine ! check_mpi_status

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Function: ncdfcheck
        ! Description:
        !
        ! Input: ncstatus - ncstatus from PnetCDF function call 
        !        message - user informational message
        !        comm_rank (optional) - mpi processor rank 
        !         
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine ncdf_check(ncstatus, message, comm_rank)
       
            implicit none
            integer :: ncstatus
            integer, optional :: comm_rank
            character(*) :: message
            
            if (ncstatus .ne. 0) then
                if (present(comm_rank)) then
                    write (6, *) "PnetCDF procedure reported error: (" , &
                        ncstatus, "), comm_rank ", comm_rank, ": ", trim(message)   
                else
                    write (6, *) "PnetCDF procedure reported error: (", &
                        ncstatus, "): ", trim(message)
                end if
            end if
        
        end subroutine ncdf_check

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Function: du
        ! Description:
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
         function derivative(h, u, v, dx, dy, df)
           
            implicit none 
            
             real, dimension(:,:), intent(inout) :: u, v, h 
             real, intent(in) :: dx, dy
             integer :: nx, ny  ! field dims
             ! return type
             real, dimension(0:size(h,dim=1)-1,0:size(h,dim=2)-1) :: derivative
            
             interface
             function df(h,u,v,dx,dy)
                real, dimension(:,:), intent(in) :: h, u, v
                real, intent(in) :: dx, dy
                real, dimension(0:size(h,dim=1)-1,0:size(h,dim=2)-1) :: df
             end function df
             end interface

             nx = size(h,dim=1)-2 ! length of size is nx + 2 halo cells
             ny = size(h,dim=2)-2

             call update_halos(u, nx, ny)
             call update_halos(v, nx, ny)
             call update_halos(h, nx, ny)
             derivative = df(h,u,v,dx,dy)

         end function derivative

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Function: update_halos 
        ! Description:
        ! Do "halo exchange" here to ensure that
        ! x(0,:) and x(nx+1,:) are updated from neighboring processors
        !
        ! 1) Every MPI task needs to know its location in the overall decomposition (PX, PY)
        ! 2) For the halo exchange, updated left column from processor PX-1
        !                           updated right column from processor PX+1
        !
        ! call MPI_Irecv(..., x(0,:), ...)     ! recv from task to left
        ! call MPI_Irecv(..., x(nx+1,:), ...)  ! recv from task to right
        !
        ! call MPI_Isend(..., x(nx,:), ...)    ! send to task to right
        ! call MPI_Isend(..., x(1,:), ...)     ! send to task to left
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
         subroutine update_halos(f, nx, ny)
           
            implicit none 
            
             type(Neighbors) :: cell_neighbors
             integer, intent(in) :: nx, ny          ! field dims
             real, dimension(0:nx+1,0:ny+1), intent(inout) :: f
             integer :: PX, PY, pleft, pright, pabove, pbelow
             integer  :: stat(8*MPI_STATUS_SIZE)
             integer :: i, j, tag, ierr
             integer :: rreq(8)
             real, dimension(nx) :: data_buf_left, data_buf_right
             real, dimension(ny) :: data_buf_upper, data_buf_lower

             ! In setup routine
             ! *) NPX - num tasks in x-direction (initially, sqrt(comm_size))
             ! *) NPY - num tasks in y-direction (initially, sqrt(comm_size))
             ! *) PX  - x-coord in decomp = (comm_rank % NPX)
             ! *) PY  - y-coord in decomp = (comm_rank / NPY)
 
             PX = nlen_subgrids(comm_size)
             PY = PX ! square MPI grid
             
            cell_neighbors = nearest_neighbors(comm_rank, comm_size) 
             
            pright = cell_neighbors%right
            pleft = cell_neighbors%left
            pbelow = cell_neighbors%below
            pabove = cell_neighbors%above

            ! halo update from left cell
            call MPI_Irecv(data_buf_left, ny, MPI_REAL, pleft, 0, &
                           MPI_COMM_WORLD, rreq(1), ierr)
            call check_mpi_status(ierr, "MPI_Irecv, left", comm_rank)

            ! halo update from right cell
            call MPI_Irecv(data_buf_right, ny, MPI_REAL, pright, 1, &
                           MPI_COMM_WORLD, rreq(2), ierr)
            call check_mpi_status(ierr, "MPI_Irecv, right", comm_rank)

            ! halo update from above cell
            call MPI_Irecv(data_buf_upper, ny, MPI_REAL, pabove, 2, &
                           MPI_COMM_WORLD, rreq(3), ierr)
            call check_mpi_status(ierr, "MPI_Irecv, right", comm_rank)

            ! halo update from below cell
            call MPI_Irecv(data_buf_lower, ny, MPI_REAL, pbelow, 3, &
                           MPI_COMM_WORLD, rreq(4), ierr)
            call check_mpi_status(ierr, "MPI_Irecv, right", comm_rank)

            call MPI_Isend(f(nx,1:ny), ny, MPI_REAL, pright, 0, &
                           MPI_COMM_WORLD, rreq(5), ierr)
            call check_mpi_status(ierr, "MPI_Send, right", comm_rank)
            call MPI_Isend(f(1,1:ny), ny, MPI_REAL, pleft, 1, &
                           MPI_COMM_WORLD, rreq(6), ierr)
            call check_mpi_status(ierr, "MPI_Send, left", comm_rank)
            call MPI_Isend(f(1:nx,1), nx, MPI_REAL, pbelow, 2, &
                           MPI_COMM_WORLD, rreq(7), ierr)
            call check_mpi_status(ierr, "MPI_Send, right", comm_rank)
            call MPI_Isend(f(1:nx,ny), nx, MPI_REAL, pabove, 3, &
                           MPI_COMM_WORLD, rreq(8), ierr)
            call check_mpi_status(ierr, "MPI_Send, left", comm_rank)
           
            !call MPI_Wait(rreq, MPI_STATUS_IGNORE, ierr)
            call MPI_Waitall(8, rreq, MPI_STATUSES_IGNORE, ierr)
            call check_mpi_status(ierr, "MPI_Waitall", comm_rank)

            ! update halos from buffers
            f(0,1:ny) = data_buf_left
            f(nx+1,1:ny) = data_buf_right
            
            f(1:nx,0) = data_buf_lower
            f(1:nx,ny+1) = data_buf_upper

         end subroutine update_halos
  
end program

