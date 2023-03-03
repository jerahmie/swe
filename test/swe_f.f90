!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! File: swe_f.f90
! Test file to call c routines from fortran
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program swe_f
  use iso_c_binding
  implicit none

  interface
    function twice_arr_c(arr, n) bind(c, name='twice_arr')
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), value :: n
      real(c_double), dimension(n) :: arr
    end function twice_arr_c

    function avg_arr_c(arr, n) bind(c, name='avg_arr')
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), value :: n
      real(c_double), dimension(n) :: arr
    end function avg_arr_c

    function twice_arr_2d_c(arr, ni, nj) bind(c, name='twice_arr_2d')
      use, intrinsic :: iso_c_binding
      integer(c_int32_t), value :: ni, nj
      real(c_double), dimension(ni,nj) :: arr
      end function twice_arr_2d_c
  end interface

  integer :: i, j
  integer(c_int32_t) :: stat
  integer, parameter :: arr_len = 10, arr2i = 5, arr2j = 3
  real(kind=8), dimension(arr_len), target :: arr_fort
  real(kind=8), dimension(arr2i, arr2j) :: arr2_fort

  do i=1,arr_len
    arr_fort(i) = real(i, kind=8)**2
  end do

  print *, "arr_fort: ", arr_fort

  stat = avg_arr_c(arr_fort, arr_len) 
  print *, "avg_fort: ", arr_fort

!  do i=1,arr2i
!    do j=1,arr2j
!      arr2_fort(i,j) = real(i+j, kind=8)
!    end do
!  end do
!      
!  stat = twice_arr_2d_c(arr2_fort, arr2i, arr2j)
!  do i=1,arr2i
!    do j=1,arr2j
!      print *, real(i+j), arr2_fort(i,j)
!    end do
!  end do
  end program swe_f
