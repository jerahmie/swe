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
  real, parameter :: dx = 2.0,  dy = 2.0 ! grid resolution
  integer, parameter :: NT=500 ! number of time steps
  real, parameter :: dt = 0.2   ! time step
  real, parameter :: g = 9.8 ! gravitational constant, Earth surface (m/s)
  real, parameter :: h0 = 1.0 ! mean fluid thickness
  real, parameter :: h1 = 0.1 ! thickness perturbation magnitude
  real, dimension(NX,NY) :: u, v, h
  integer :: i

  print *, "shallow wave equations"

  ! initialize h, v, u
  h = h0 + h1*gaussian2d(10.0, 0.5*real(NX)*dx, 0.5*real(NY)*dy, dx, dy, NX, NY)
  v = 0.0
  u = 0.0
  !call set_initial_state()

  ! initialize observer
  call observer_init(NX, NY)

  ! update h, v
  do i=1,NT
    u = u + dt*du(h,u,v,dx,dy)
    !u = euler_forward(u, dt, du)
    v = v + dt*dv(h,u,v,dx,dy)
    !v = euler_forward(v, dv, dt)
    h = h + dt*dh(h,u,v,dx,dy)
    !h = euler_forward(h, dh)

    ! write wave to file
    call observer_write(h, i)
  end do
     
  ! finalize observer
  call observer_finalize() 

end program swe

