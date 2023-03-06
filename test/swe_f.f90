!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! File: swe_f.f90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program swe_f
  use iso_c_binding

  implicit none

  integer :: i, j
  integer :: istat
  integer(c_int32_t) :: icount
  integer, parameter :: NX = 10, NY=5
  type(c_ptr), pointer :: arr2d_ptr
  real(c_double), dimension(NX,NY), target :: arr2d

  interface
    function say_hi() bind(c)
      use, intrinsic :: iso_c_binding
      integer(c_int) :: say_hi
    end function say_hi

    subroutine print_mat(arr2d, ni, nj ) bind(c, name="print_mat")
      use, intrinsic  :: iso_c_binding
      integer(c_int32_t), value :: ni, nj
      real(c_double), dimension(:,:), pointer :: arr2d
      end subroutine
  end interface


  icount = 0
  do i=1,NX
    do j=1,NY
      arr2d(i,j) = icount
      icount = icount + 1
    end do
  end do

  istat = say_hi() 

  !print *, arr2d
  !arr2d_ptr = c_loc(arr2d)
  !call print_mat(arr2d_ptr, NX, NY)
 
end program swe_f

