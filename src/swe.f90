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
  character(len=32) :: fmt1

  print *, "shallow wave equations"
  fmt1 =  '(a, "->", i4, "/", i4)'

  ! initialize state variables 
  call set_initial_state(h, u, v)

  ! initialize observer
  call observer_init(NX, NY)

  ! update h, v
  do i=1,NT
    !write (*,'(a, i3, " /", i3)', advance='no') achar(13), i, NT
    write (*,fmt1, advance='no') achar(13), i, NT
    u = u + euler_forward(h,u,v,dx,dy,dt,du)
    v = v + euler_forward(h,u,v,dx,dy,dt,dv)
    h = h + euler_forward(h,u,v,dx,dy,dt,dh)

    ! write wave to file
    call observer_write(h, i)
  end do
     
  ! finalize observer
  write (*,*) achar(10)
  call observer_finalize() 
end program swe

