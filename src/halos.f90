program swe_mpi
    use mpi
    use pnetcdf
    use init_par
    use swe_mpi_help
    use equations_par

    implicit none

    type(Neighbors) :: cell_neighbors
    integer :: num_args
    integer :: comm_rank, comm_size, ncid, ncstatus
    integer :: i, j, it,  ndims, dimids(2), varid
    integer :: local_nx, local_ny, nsubgrid
    character(255) :: ncfile_input, ncfile_output
    character(len=32), dimension(:), allocatable :: dimnamei ! varname
    integer(kind=MPI_OFFSET_KIND), dimension(:), allocatable :: dimval
    logical :: res
    real(kind=4), allocatable, dimension(:,:) :: h, dhx, u, v
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

    ! execute solver
    do it=1,NT
        call update_h(h,u,v,local_nx,local_ny,dx,dy,dt)
    end do

    ! write results
    ! TODO: move to output module 
    ! create new netcdf file for read
    ncstatus = nf90mpi_create(mpi_comm=MPI_COMM_WORLD, path=trim(ncfile_output), &
                              cmode=NF_CLOBBER, mpi_info=MPI_INFO_NULL, ncid=ncid)
    call ncdf_check(ncstatus, "nf90mpi_open for write", comm_rank)
    ncstatus = nfmpi_def_dim(ncid, "ny", int(ny, kind=MPI_OFFSET_KIND), dimids(2))
    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", comm_rank)
    ncstatus = nfmpi_def_dim(ncid, "nx", int(nx, kind=MPI_OFFSET_KIND), dimids(1))
    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", comm_rank)
    
    ncstatus = nfmpi_def_var(ncid, "var", NF90_FLOAT, 2, dimids, varid)
    call ncdf_check(ncstatus, "nf90mpi_def_var for write", comm_rank)

    ncstatus = nf90mpi_enddef(ncid)
    call ncdf_check(ncstatus, "nf90mpi_enddef for write", comm_rank)
    if (comm_rank .eq. 0) then
        print*, "dimids: ", dimids, " varid: ", varid
    endif
    
    ncstatus = nfmpi_put_vara_real_all(ncid, i, starts, counts, h(1:local_nx,1:local_ny))
    call ncdf_check(ncstatus, "nfmpi_put_vara_real_all for write", comm_rank)
    
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
        ! Function: ncdfcheck
        ! Description:
        !
        ! Input: ncstatus - ncstatus from PnetCDF function call 
        !        message - user informational message
        !        comm_rank (optional) - mpi processor rank 
        !         
        ! Returns: Exits program if ncstatus is not zero
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
        ! Function: update_h
        ! Description:
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
         subroutine update_h(h,u,v,nx,ny,dx,dy,dt)
           
            implicit none 
            
             type(Neighbors) :: cell_neighbors
             integer, intent(in) :: nx, ny          ! field dims
             real, intent(in) :: dx, dy             ! grid resolution
             real, intent(in) :: dt                 ! temporal resolution
             integer :: PX, PY, pleft, pright, pabove, pbelow
             real, dimension(0:nx+1,0:ny+1), intent(inout) :: h, u, v
             !real, dimension(0:nx+1,0:ny+1), intent(in) :: u, v 
             integer  :: stat(MPI_STATUS_SIZE)
             integer :: i, j, tag, ierr, rreq
             
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
            !print *, comm_rank, "right: ", cell_neighbors%right, "pleft: ", cell_neighbors%left, &
            !                     "pbelow: ", cell_neighbors%below, "pabove: ", cell_neighbors%above
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
            ! Halo cells u  
            call MPI_Irecv(u(0,:), ny+2, MPI_REAL, pleft, pright, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Irecv(u(nx+1,:), ny+2, MPI_REAL, pright, pleft, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Irecv(u(:,0), nx+2, MPI_REAL, pbelow, pabove, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Irecv(u(:,ny+1), nx+2, MPI_REAL, pabove, pbelow, &
                           MPI_COMM_WORLD, rreq, ierr)

            call MPI_Isend(u(nx,:), ny+2, MPI_REAL, pright, pright, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Isend(u(1,:), ny+2, MPI_REAL, pleft, pleft, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Isend(u(:,ny), nx+2, MPI_REAL, pabove, pabove, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Isend(u(:,1), nx+2, MPI_REAL, pbelow, pbelow, &
                           MPI_COMM_WORLD, rreq, ierr)
            ! Halo cells v
            call MPI_Irecv(v(0,:), ny+2, MPI_REAL, pleft, pright, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Irecv(v(nx+1,:), ny+2, MPI_REAL, pright, pleft, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Irecv(v(:,0), nx+2, MPI_REAL, pbelow, pabove, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Irecv(v(:,ny+1), nx+2, MPI_REAL, pabove, pbelow, &
                           MPI_COMM_WORLD, rreq, ierr)

            call MPI_Isend(v(nx,:), ny+2, MPI_REAL, pright, pright, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Isend(v(1,:), ny+2, MPI_REAL, pleft, pleft, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Isend(v(:,ny), nx+2, MPI_REAL, pabove, pabove, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Isend(v(:,1), nx+2, MPI_REAL, pbelow, pbelow, &
                           MPI_COMM_WORLD, rreq, ierr)
            ! Halo cells h
            call MPI_Irecv(h(0,:), ny+2, MPI_REAL, pleft, pright, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Irecv(h(nx+1,:), ny+2, MPI_REAL, pright, pleft, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Irecv(h(:,0), nx+2, MPI_REAL, pbelow, pabove, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Irecv(h(:,ny+1), nx+2, MPI_REAL, pabove, pbelow, &
                           MPI_COMM_WORLD, rreq, ierr)

            call MPI_Isend(h(nx,:), ny+2, MPI_REAL, pright, pright, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Isend(h(1,:), ny+2, MPI_REAL, pleft, pleft, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Isend(h(:,ny), nx+2, MPI_REAL, pabove, pabove, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Isend(h(:,1), nx+2, MPI_REAL, pbelow, pbelow, &
                           MPI_COMM_WORLD, rreq, ierr)
            call MPI_Wait(rreq, stat, ierr)
            
            u = u + euler_forward(h,u,v,dx,dy,dt,du)
            v = v + euler_forward(h,u,v,dx,dy,dt,dv)
            h = h + euler_forward(h,u,v,dx,dy,dt,dh)

        end subroutine update_h
   
end program

