module observer_par

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Module: observer_par
! Obvserve height data from shallow water equation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    private

    !integer, private :: ncid



contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Function observer_par_init
! Initialize PNetCDF output file. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function observer_par_init(ncfile_name)

    use mpi
    use pnetcdf 

    implicit none

    ncstatus = nf90mpi_create(mpi_comm=MPI_COMM_WORLD, path=trim(ncfile_output), &
                              cmode=NF_CLOBBER, mpi_info=MPI_INFO_NULL, ncid=ncid)
    call ncdf_check(ncstatus, "nf90mpi_open for write", comm_rank)
    ncstatus = nfmpi_def_dim(ncid, "ny", int(ny, kind=MPI_OFFSET_KIND), nyid)
    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", comm_rank)
    ncstatus = nfmpi_def_dim(ncid, "nx", int(nx, kind=MPI_OFFSET_KIND), nxid)
    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", comm_rank)
    ncstatus = nfmpi_def_dim(ncid, "time", nf90mpi_unlimited, timedimid)
    call ncdf_check(ncstatus, "nf90mpi_def_dim for write", comm_rank)
    
    ncstatus = nfmpi_def_var(ncid, "var", NF90_FLOAT, 2, (/nxid, nyid, timedimid/), varid)
    call ncdf_check(ncstatus, "nf90mpi_def_var for write", comm_rank)

    ncstatus = nf90mpi_enddef(ncid)
    call ncdf_check(ncstatus, "nf90mpi_enddef for write", comm_rank)
    if (comm_rank .eq. 0) then
        print*, "dimids: ", nxid, nyid, " varid: ", varid
    endif
end function observer_par_init

! subroutine observer_par_write(ncid)
!

! function observer_par_finalize
! end function observer_par_finalize

end module observer_par