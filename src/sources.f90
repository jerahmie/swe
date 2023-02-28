module sources
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module: sources
!
! Functions to initialize h 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  private

  public :: gaussian2d

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !
  ! gaussian2d -- generate a 2D gaussian distribution on cartesian 
  !               coordinate system
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  pure function gaussian2d(sigma, x0, y0, dx, dy, NX, NY)

    implicit none

    integer, intent(in) :: NX, NY
    real, intent(in) :: sigma, x0, y0, dx, dy
    integer :: i, j

    ! return value
    real, dimension(NX,NY) :: gaussian2d

    do i = 1, NX
      do j = 1, NY
        gaussian2d(i,j) = exp(-1.0*(((real(i)*dx-x0)**2)/(2.0*sigma**2) + &
                                   ((real(j)*dy-y0)**2)/(2.0*sigma**2)))
      end do !j
    end do !i
  end function gaussian2d
  
end module sources
