module regrid_grid1d_class

implicit none ; private

public :: grid1d_t, grid1d_init, grid1d_destroy

integer, parameter :: max_grids = 100, max_nbcells = 1000

! -----------------------------------------------------------------------------
! Definition of the one-dimensional grid structure
! -----------------------------------------------------------------------------
type grid1d_t
  ! Number of grid cells
  integer                           :: nb_cells;
  
  ! Cell widths (size = nb_cells)
  real, dimension(max_nbcells)   :: h;
  
  ! Coordinates of nodes (size = nb_cells+1)
  real, dimension(max_nbcells)   :: x;   

end type grid1d_t 



!type(grid1d_t), target, dimension(max_grids) :: grid1d

!integer :: num_grids = 0
logical :: is_initialized = .false.

contains

!------------------------------------------------------------------------------
! grid_init
! -----------------------------------------------------------------------------
subroutine grid1d_init ( grid, nb_cells, ngrid)
!------------------------------------------------------------------------------
! Initialization (memory allocation) of a grid
!------------------------------------------------------------------------------

  type(grid1d_t), intent(inout)  :: grid;
  integer, intent(in)            :: nb_cells;
  integer, intent(out)           :: ngrid


         
!  ngrid=num_grids+1
!  num_grids=num_grids+1
!  ngrid=num_grids

!  grid=>grid1d(ngrid)
  
  ! Set the number of cells
  grid%nb_cells = nb_cells;

!  Memory allocation
!  allocate ( grid%h(nb_cells) );
!  allocate ( grid%x(nb_cells+1) );

  ! Set all entries to zero
  grid%h = 0.0
  grid%x = 0.0

end subroutine grid1d_init

!------------------------------------------------------------------------------
! grid_destroy
! -----------------------------------------------------------------------------
subroutine grid1d_destroy ( grid )
!------------------------------------------------------------------------------
! Reclaim previously allocated memory for a grid
!------------------------------------------------------------------------------

  type(grid1d_t), intent(inout) :: grid;

  
!  deallocate ( grid%x );

end subroutine grid1d_destroy

end module regrid_grid1d_class
