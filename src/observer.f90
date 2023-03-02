module observer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Module: observe
!
! The module is used to observe height and state data and to write them 
! out to a netcdf file.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  private 

  integer, private :: ncid, nxid, nyid, hvarid, timedimid

  public :: observer_init, observer_write, observer_finalize

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! 
  ! Name: observer_init
  ! 
  ! Input: NX -- size of array in x-direction
  !        NY -- size of array in y-direction
  ! Output: none
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine observer_init(NX, NY)

    use netcdf

    implicit none

    character (len=*), parameter :: filename = "observer.nc"
    integer :: istatus, NX, NY

    istatus = nf90_create(filename, nf90_clobber, ncid)
    if (istatus /= 0) stop "Unable to open netcdf file."

    ! Define dimensions
    istatus = nf90_def_dim(ncid, "NX", NX, nxid)
    if (istatus /= 0) stop "Unable to create NX dimension."
    istatus = nf90_def_dim(ncid, "NY", NY, nyid)
    if (istatus /= 0) stop "Unable to create NY dimension."
    istatus = nf90_def_dim(ncid, "time", nf90_unlimited, timedimid)
    if (istatus /= 0) stop "Unable to create time dimension."

    ! Define variables
!    istatus = nf90_def_var(ncid, "u", nf90_double, &
!                                (/ NX, NY, timedimid /), uvarid)
!    if (istatus /= 0) stop "Unable to create state parameter variable id."
!    istatus = nf90_def_var(ncid, "v", nf90_double, &
!                                (/ NX, NY, timedimid /), vvarid)
!    if (istatus /= 0) stop "Unable to create state parameter variable id."
    istatus = nf90_def_var(ncid, "h", nf90_double, &
                                (/ nxid, nyid, timedimid /), hvarid)
    if (istatus /= 0) stop "Unable to create state parameter variable id."

    istatus = nf90_enddef(ncid)
    if (istatus /= 0) stop "Could not exit define mode."

  end subroutine observer_init

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! 
  ! Name: observer_write
  ! 
  ! Input: u --
  !         v --
  !         h --
  !         i -- time index
  ! Output: none
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine observer_write(h, i)

    use netcdf

    implicit none

    real, dimension(:,:), intent(in) :: h
    integer :: istatus, i

    istatus = nf90_put_var(ncid, hvarid, h, start=(/ 1, 1, i/))
    if (istatus /= 0) then
      print *, "i: ", i
      print *, "h: ", h
      print *, "Status: ", istatus
      stop "Unable to write variable to file."
    end if

  end subroutine observer_write

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! 
  ! Name: observer_finalize
  ! 
  ! Input: 
  ! Output:
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine observer_finalize()

    use netcdf

    implicit none

    integer :: istatus

    istatus = nf90_close(ncid)
    if (istatus /= 0) stop "Unable to close netcdf file."

  end subroutine observer_finalize

end module observer

