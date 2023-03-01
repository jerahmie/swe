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
  real, dimension(NX,NY) :: u, v, h
  integer :: i

  print *, "shallow wave equations"

  ! initialize state variables 
  call set_initial_state(h, u, v)

  ! initialize observer
  call observer_init(NX, NY)

  ! update h, v
  do i=1,NT

    u = u + euler_forward(h,u,v,dx,dy,dt,du)
    v = v + euler_forward(h,u,v,dx,dy,dt,dv)
    h = h + euler_forward(h,u,v,dx,dy,dt,dh)

    ! write wave to file
    call observer_write(h, i)
  end do
     
  ! finalize observer
  call observer_finalize() 

end program swe

