module equations_par
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module equations.f90
!
! Ordinary Differential Equation Solver for use with Shallow Water
! Equation simulations.
! 
! Jerahmie Radder, 24Feb2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  private

    real, parameter :: g = -9.8 ! acceleration due to gravity (m/s)

  public :: diff_center_x, diff_center_y, euler_forward, du, dv, dh, set_initial_state

  real, public :: dx = 2.0, dy = 2.0, dt = 0.2  ! grid resolution
  integer, public :: NX, NY, NT  ! grid dimensions

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function: du
  ! Description: Returns the time derivative of state variable u
  !              
  ! Input: h -- fluid thickness
  !        u -- x-directed velocity
  !        v -- y-directed velocity
  !        dx -- x grid resolution
  !        dy -- y grid resolution
  ! Returns: du
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  pure function du (h, u, v, dx, dy)
    
    implicit none

    real, dimension(:,:), intent(in) :: h, u, v 
    real, intent(in) :: dx, dy

    ! return value du
    real, dimension(size(u,dim=1), size(u,dim=2)) :: du

    du = g*diff_center_x(h,dx) - & 
         u*diff_center_x(u,dx) - &
         v*diff_center_y(u,dy)

  end function du
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function: dv
  ! Description: Returns the time derivative of the state variable v
  !              
  ! Input: h -- fluid thickness
  !        u -- x-directed velocity
  !        v -- y-directed velocity
  !        dx -- x grid resolution
  !        dy -- y grid resolution
  ! Returns: dv
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  pure function dv(h, u, v, dx, dy)

    implicit none 
    
    real, dimension(:,:), intent(in) :: h, u, v
    real, intent(in) :: dx, dy

    ! return value
    real, dimension(size(v,dim=1),size(v,dim=2)) :: dv
    dv = g*diff_center_y(h,dy) - &
         u*diff_center_x(v,dx) - &
         v*diff_center_y(v,dy)

  end function dv

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Function: dh
  ! Description: Returns the time derivative of fluid thickness h
  !              
  ! Input: h -- fluid thickness
  !        u -- x-directed velocity 
  !        v -- y-directed velocity
  !        dx -- x grid resolution
  !        dy -- y grid resolution
  ! Returns: dh
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pure function dh(h, u, v, dx, dy)

    implicit none
    real, dimension(:,:), intent(in) :: h, u, v
    real, intent(in) :: dx, dy
    
    ! return value
    real, dimension(size(h,dim=1),size(h,dim=2)) :: dh 
    dh = -1.0*(diff_center_x(u*h,dx) + diff_center_y(v*h,dy))

  end function dh

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Function: diff_center_x
  ! Description: Returns central difference in x-direction of array 2d 
  !              with periodic boundary conditions.
  ! Input: x -- 2d vector of values
  ! Returns: dx 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  pure function diff_center_x(r,dx)

    implicit none
  
    real, dimension(:,:), intent(in) :: r
    real, intent(in) :: dx
    integer :: xdim

    ! return value
    real, dimension(size(r,dim=1),size(r,dim=2)) :: diff_center_x
  
    xdim = size(r,dim=1)
    diff_center_x(1,:) = (0.5*dx) * (r(2,:) - r(xdim,:))
    diff_center_x(xdim,:) = (0.5*dx) * (r(1,:) - r(xdim-1,:))
    diff_center_x(2:xdim-1,:) = (0.5*dx) * (r(3:xdim,:) - r(1:xdim-2,:))
    
  end function diff_center_x


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Function: diff_center_y
  ! Description: Returns central difference in y-direction of array 2d 
  !              with periodic boundary conditions.
  ! Input: r -- 2d vector of values
  !        dy -- grid y-resolution
  ! Returns: diff_center_y -- center difference of  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  pure function diff_center_y(r, dy)

    implicit none
  
    real, dimension(:,:), intent(in) :: r
    real, intent(in) :: dy
    integer :: ydim

    ! return value
    real, dimension(size(r,dim=1),size(r,dim=2)) :: diff_center_y
  
    ydim = size(r, dim=2) 
    diff_center_y(:,1) = (0.5*dy) * (r(:,2) - r(:,ydim))
    diff_center_y(:,ydim) = (0.5*dy) * (r(:,1) - r(:,ydim-1))
    diff_center_y(:,2:ydim-1) = (0.5*dy) * (r(:,3:ydim) - r(:,1:ydim-2))
    
  end function diff_center_y


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Function: euler_forward
  ! Description: Returns forward difference of 2d array
  !              with periodic boundary conditions.
  ! Input: h -- fluid thickness
  !        u -- x-directed velocity 
  !        v -- y-directed velocity
  !        dx -- x grid resolution
  !        dy -- y grid resolution
  !        dt -- delta t 
  ! Returns: i
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !function euler_forward(h, u, v, dx, dy, dt, dh)
  pure function euler_forward(h, u, v, dx, dy, dt, df)
     
    implicit none
 
    real, dimension(:,:), intent(in) :: h, u, v
    real, intent(in) :: dx, dy, dt
    ! return value
    real, dimension(size(h,dim=1),size(h,dim=2)) :: euler_forward

    interface 
     pure function df(h,u,v,dx,dy)
        real, dimension(:,:), intent(in) :: h, u, v
        real, intent(in) :: dx, dy
        real, dimension(size(h,dim=1),size(h,dim=2)) :: df
     end function df
    end interface
    euler_forward = dt*df(h,u,v,dx,dy)
  
  end function euler_forward

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  ! Name: set_initial_state
  !
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  subroutine set_initial_state(h, u, v)
    use sources, only: gaussian2d 
    implicit none

    real, dimension(:,:) :: h, u, v
  
    h = 1.0 + 0.1 * gaussian2d(10.0, 0.5*real(NX)*dx, 0.5*real(NY)*dy, dx, dy, NX, NY)  
    u = 0.0    
    v = 0.0 
  end subroutine set_initial_state

end module equations_par 

