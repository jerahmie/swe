!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module: init_par
!
! Parallel initialization routines for Parallel SWE solver.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module init_par

  private

  public :: nlen_subgrids, Neighbors, nearest_neighbors

  type Neighbors
    integer :: right, left, below, above
  end type Neighbors


contains

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Function: nlen_subgrids
        ! Description: returns number of subgrids on side of in square grids
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pure function nlen_subgrids(nprocs)

            implicit none

            integer, intent(in) :: nprocs
            integer :: nlen_subgrids

            nlen_subgrids = int(sqrt(real(nprocs)))

        end function nlen_subgrids

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Function: nearest_neighbors
        ! Description: Returns struct with nearest neighbors on a cyclic grid
        !
        ! Input:   i -- index of compute cell
        !          N -- total number of cells in square grid
        ! Returns: derived-type that contains the index of the cell neighbor
        !          to the right, left, above and below
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pure function nearest_neighbors(i,N)

            implicit none

            integer, intent(in) :: i, N
            integer :: ni

            ! return value: nearest neighbors
            type(Neighbors) :: nearest_neighbors

            ni = int(sqrt(real(N,kind=8)))

            ! note: (i/ni)*ni = integer equivalent of floor(i/ni)*ni
            nearest_neighbors%right = mod(i+1,ni) + (i/ni)*ni
            nearest_neighbors%left = mod(i+ni-1,ni) + (i/ni)*ni
            nearest_neighbors%below = mod(i+N-ni, N)
            nearest_neighbors%above = mod(i+N+ni, N)

        end function nearest_neighbors

end module init_par
