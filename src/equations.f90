module equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module equations.f90
!
! Ordinary Differential Equation Solver for use with Shallow Water
! Equation simulations.
! 
! Jerahmie Radder, 24Feb2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  private

  public :: diff_center_x, diff_center_y, diff_forward

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Function: diff_center_x
  ! Description: Returns central difference in x-direction of array 2d 
  !              with periodic boundary conditions.
  ! Input: x -- 2d vector of values
  ! Returns: dx 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  pure function diff_center_x(x)

    implicit none
  
    real, dimension(:,:), intent(in) :: x
    integer :: xdim

    ! return value
    real, dimension(size(x,dim=1),size(x,dim=2)) :: diff_center_x
  
    xdim = size(x,dim=1)
    diff_center_x(1,:) = 0.5*(x(2,:) - x(xdim,:))

    diff_center_x(xdim,:) = 0.5*(x(1,:) - x(xdim-1,:))
    diff_center_x(2:xdim-1,:) = 0.5*(x(3:xdim,:) - x(1:xdim-2,:))
    
  end function diff_center_x


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Function: diff_center_y
  ! Description: Returns central difference in y-direction of array 2d 
  !              with periodic boundary conditions.
  ! Input: x -- 2d vector of values
  ! Returns: dx 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  pure function diff_center_y(x)

    implicit none
  
    real, dimension(:,:), intent(in) :: x
    integer :: xdim

    ! return value
    real, dimension(size(x,dim=1),size(x,dim=2)) :: diff_center_y
  
    xdim = size(x, dim=2) 
    diff_center_y(:,1) = 0.5*(x(:,2) - x(:,xdim))
    diff_center_y(:,xdim) = 0.5*(x(:,1) - x(:,xdim-1))
    diff_center_y(:,2:xdim-1) = 0.5*(x(:,3:xdim) - x(:,1:xdim-2))
    
  end function diff_center_y


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Function: euler_forward
  ! Description: Returns forward difference of 2d array
  !              with periodic boundary conditions.
  ! Input: x -- 2d vector of values
  ! Returns: dx 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  function euler_forward(x)
     
     implicit none
 
     real, dimension(:,:), intent(in) :: x
     
     ! return value
     real, dimension(size(x,dim=1),size(x,dim=2)) :: diff_forward
     !
     ! code
     !
     euler_forward = 0.0
 
   end function euler_forward

end module equations 

