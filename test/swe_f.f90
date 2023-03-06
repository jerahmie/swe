!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! File: swe_f.f90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program swe_f
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_ptr

  implicit none

  interface
    function say_hi() bind(c)
      import
      integer(c_int) :: say_hi
    end function say_hi

    integer(c_int) function twice_arr(arr, arr_size) bind(c, name="twice_arr")
      import
      implicit none
      integer(c_int), value :: arr_size
      real(c_double), dimension(arr_size) :: arr
      
      end function twice_arr

    subroutine print_matf(arr, ni, nj ) bind(c, name="print_matf")
      import
      integer(c_int), value :: ni
      integer(c_int), value :: nj
      real(c_double), dimension(*) :: arr
      end subroutine print_matf

    integer(c_int) function twice_arr_2d(arr2d, ni, nj) bind(c, name="twice_arr_2df")
      import 
      integer(c_int), value :: ni, nj
      real(c_double), dimension(*) :: arr2d
    end function twice_arr_2d
 
  end interface

  integer :: i, j
  integer :: istat
  integer(c_int) :: icount
  integer(c_int), parameter :: NX = 10, NY=5
  real(c_double), dimension(NX) :: test1
  real(c_double), dimension(NX,NY) :: arr2d

  ! If this doesn't work, nothing will
  istat = say_hi()

  ! Populate and pass a 1D array to C-function 
  do i = 1, NX
    test1(i) = real(i, kind=c_double)
  end do
  
  call print_matf(test1, NX, 1) 
  istat = twice_arr(test1, NX )
  call print_matf(test1, NX, 1)


  ! Populate and pass a 2D array to C-function 
  icount = 1
  do i=1,NX
    do j=1,NY
      arr2d(i,j) = real(icount, kind=c_double)
      icount = icount + 1
    end do
  end do
  call print_matf(arr2d, NX, NY)
  istat = twice_arr_2d(arr2d, NX, NY)  
  call print_matf(arr2d, NX, NY)
 

end program swe_f

