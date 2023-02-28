!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! swe - shallow wave equation solver
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program swe
  use equations
  use observer
  use sources, only : gaussian2d

  implicit none
  integer, parameter :: NX=201, NY=201  ! grid dimensions, NxM
  real, parameter :: dx = 1,  dy = 1 ! grid resolution
  integer, parameter :: NT=1000 ! number of time steps
  real, parameter :: dt = 0.2   ! time step
  real, parameter :: g = 9.8 ! gravitational constant, Earth surface (m/s)
  real, parameter :: h0 = 1.0 ! mean fluid thickness
  real, parameter :: h1 = 0.1 ! thickness perturbation magnitude
  real, dimension(NX,NY) :: u, v, h
  integer :: i

  print *, "shallow wave equations"

  ! initialize h, v
  h = h0 + h1*gaussian2d(10.0, 0.25*real(NX)*dx/2.0, 0.75*real(NY)*dy/2.0, dx, dy, NX, NY)
  v = 0.0
  u = 0.0

  ! initialize observer
  call observer_init(NX, NY)

  ! update h, v
  do i=1,NT
    u = u - dt*(g*(diff_center_x(h)/dx) + &
               u*diff_center_x(u)/dx + &
               v*diff_center_y(u)/dy)
    v = v - dt*(g*(diff_center_y(h)/dy) + &
               u*diff_center_x(v)/dx + &
               v*diff_center_y(v)/dy)
    h = h - dt*(diff_center_x(u*(h))/dx + &
               diff_center_y(v*(h))/dy)

    ! write wave to file
    call observer_write(h, i)
  end do
     
  ! finalize observer
  call observer_finalize() 

end program swe

