module observer_par

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Module: observer_par
! Obvserve height data from shallow water equation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    private

    integer, private :: comm_rank, comm_size, ncstatus
    integer, private :: nxid, nyid, timedimid, varid, ncid

    public :: observer_par_init, observer_par_finalize, ncdf_check

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine observer_par_init
! Initialize PNetCDF output file. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine observer_par_init(ncfile_output, nx, ny)

    use mpi
    use pnetcdf 

    implicit none

    character(255), intent(in) :: ncfile_output
    integer, intent(in) :: nx, ny

    ncstatus = nf90mpi_create(mpi_comm=MPI_COMM_WORLD, path=trim(ncfile_output), &
                              cmode=NF90_CLOBBER, mpi_info=MPI_INFO_NULL, ncid=ncid)
    call ncdf_check(ncstatus, "nf90mpi_open for write", comm_rank)
    ncstatus = nfmpi_def_dim(ncid, "ny", int(ny, kind=MPI_OFFSET_KIND), nyid)
    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", comm_rank)
    ncstatus = nfmpi_def_dim(ncid, "nx", int(nx, kind=MPI_OFFSET_KIND), nxid)
    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", comm_rank)
    ncstatus = nfmpi_def_dim(ncid, "time", nf90mpi_unlimited, timedimid)
    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", comm_rank)
    
    !ncstatus = nfmpi_def_var(ncid, "var", NF90_FLOAT, 2, (/nxid, nyid, timedimid/), varid)
    ncstatus = nfmpi_def_var(ncid, "var", NF90_DOUBLE, 3, (/nxid, nyid, timedimid/), varid)
    call ncdf_check(ncstatus, "nf90mpi_def_var for write", comm_rank)

    ncstatus = nf90mpi_enddef(ncid)
    call ncdf_check(ncstatus, "nf90mpi_enddef for write", comm_rank)
    if (comm_rank .eq. 0) then
        print*, "dimids: ", nxid, nyid, " varid: ", varid
    endif
end subroutine observer_par_init

! subroutine observer_par_write(ncid)
!

subroutine observer_par_finalize()
    use mpi
    use pnetcdf

    implicit none

    ncstatus = nf90mpi_close(ncid)
    call ncdf_check(ncstatus, "nf90mpi_close for write", comm_rank)

end subroutine observer_par_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine: ncdfcheck
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

end module observer_par