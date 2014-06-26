module MIDAS_remapping

use regrid_grid1d_class ! see 'regrid_grid.F90'
use regrid_ppoly_class  ! see 'regrid_ppoly.F90'
use regrid_polynomial   ! see 'regrid_polynomial.F90'
use regrid_edge_values  ! see 'regrid_edge_values.F90'
use regrid_edge_slopes  ! see 'regrid_edge_slopes.F90'
use regrid_pcm          ! see 'regrid_pcm.F90'
use regrid_plm          ! see 'regrid_plm.F90'
use regrid_ppm          ! see 'regrid_ppm.F90'
use regrid_pqm          ! see 'regrid_pqm.F90'
use regrid_p1m          ! see 'regrid_p1m.F90'
use regrid_p3m          ! see 'regrid_p3m.F90'
use regrid_defs         ! see 'regrid_defs.F90' (contains types and parameters)


implicit none ; private


public remapping_integration, inflate_vanished_layers

  contains

  

subroutine remapping_integration ( grid0, u0, ppoly0, grid1, u1, method )
  
  ! Arguments
  type(grid1d_t), intent(in)      :: grid0;   ! source grid
  real, dimension(:), intent(in)    :: u0;      ! source cell averages
  type(ppoly_t), intent(inout) :: ppoly0;  ! source piecewise polynomial
  type(grid1d_t), intent(inout)        :: grid1;   ! target grid
  real, dimension(:), intent(inout) :: u1;      ! target cell averages
  integer                           :: method;  ! remapping scheme to use
  
  ! Local variables
  integer       :: i, j, k
  integer       :: n0, n1;      ! nb of cells in grid0 and grid1, respectively
  integer       :: j0, j1;      ! indexes of source cells containing target 
                                ! cell edges
  real          :: x0, x1;      ! coordinates of target cell edges  
  real          :: q0, q1;      ! partially integrated quantities in source 
                                ! cells j0 and j1
  real          :: q;           ! complete integration
  real          :: a, b;        ! interval of integration (global coordinates)
  real          :: xi0, xi1;    ! interval of integration (local -- normalized 
                                ! -- coordinates)

  ! A priori, both grids contains the same number of cells but, who knows...
  n0 = grid0%nb_cells
  n1 = grid1%nb_cells

  ! Loop on cells in target grid (grid1). For each target cell, we need to find
  ! in which source cells the target cell edges lie. The associated indexes are 
  ! noted j0 and j1.
  do i = 1,grid1%nb_cells
    ! Determine the coordinates of the target cell edges
    x0 = grid1%x(i)
    x1 = grid1%x(i+1)

    ! ============================================================
    ! Check whether target cell is vanished. If it is, the cell
    ! average is simply the interpolated value at the location
    ! of the vanished cell. If it isn't, we need to integrate the
    ! quantity within the cell and divide by the cell width to
    ! determine the cell average.
    ! ============================================================
    ! 1. Cell is vanished
    if ( abs(x0 - x1) .EQ. 0.0 ) then
    
      j0 = -1
    
      do j = 1,grid0%nb_cells
        ! Left edge is found in cell j
        if ( ( x0 .GE. grid0%x(j) ) .AND. ( x0 .LE. grid0%x(j+1) ) ) then
          j0 = j
          exit; ! once target grid cell is found, exit loop
        end if
      end do

      ! If, at this point, j0 is equal to -1, it means the vanished
      ! cell lies outside the source grid. In other words, it means that
      ! the source and target grids do not cover the same physical domain
      ! and there is something very wrong !
      if ( j0 .EQ. -1 ) then
        write(*,*) i 
        print *, 'The location of the vanished cell could &
             not be found in "remapping_integration"'
        call abort()
      end if

      ! We check whether the source cell (i.e. the cell in which the
      ! vanished target cell lies) is vanished. If it is, the interpolated 
      ! value is set to be mean of the edge values (which should be the same).
      ! If it isn't, we simply interpolate.
      if ( grid0%h(j0) .EQ. 0.0 ) then
        u1(i) = 0.5 * ( ppoly0%E(j0,1) + ppoly0%E(j0,2) )
      else
        xi0 = x0 / grid0%h(j0) - grid0%x(j0) / grid0%h(j0)
    
        select case ( method )
        
          case ( INTEGRATION_PCM )   
            u1(i) = ppoly0%coefficients(j0,1)
            
          case ( INTEGRATION_PLM )  
            
            u1(i) = evaluation_polynomial ( ppoly0%coefficients(j0,:), 2, xi0 )
        
          case ( INTEGRATION_PPM )
            
            u1(i) = evaluation_polynomial ( ppoly0%coefficients(j0,:), 3, xi0 )
        
          case ( INTEGRATION_PQM )
            
            u1(i) = evaluation_polynomial ( ppoly0%coefficients(j0,:), 5, xi0 )
            
          case default
        
          print *,'The selected integration method is invalid'
          call abort()    
        
        end select   
        
      end if ! end checking whether source cell is vanished
    
    ! 2. Cell is not vanished
    else

    ! Find the cells in source grid containing the target cell edges
    j0 = -1
    j1 = -1
    
    do j = 1,grid0%nb_cells
    
        ! Left edge is found in cell j
        if ( ( x0 .GE. grid0%x(j) ) .AND. ( x0 .LE. grid0%x(j+1) ) ) then
            j0 = j
            exit;   ! once target grid cell is found, exit loop
        end if
        
    end do  

    do j = 1,grid0%nb_cells
        
        ! Right edge is found in cell j
        if ( ( x1 .GE. grid0%x(j) ) .AND. ( x1 .LE. grid0%x(j+1) ) ) then
            j1 = j
            exit;   ! once target grid cell is found, exit loop
        end if
        
    end do ! end loop on source grid cells

    if (j0 == -1 .or. j1 == -1) then
        print *, 'Unable to find target grid edges '
        print *, 'grid0%x = ',grid0%x
        print *, 'x0, x1 = ',x0, x1
        print *, 'j0,j1= ',j0, j1
        call abort
    endif
    
    ! Here, we make sure that the boundary edges of boundary cells
    ! coincide
    if ( i .EQ. 1 ) then
      j0 = 1
    end if  
    
    if ( i .EQ. grid1%nb_cells ) then
      j1 = grid0%nb_cells  !!!! mjh changed grid1% to grid0% !!!!!
    end if  

    ! To integrate, two cases must be considered: (1) the target cell is
    ! entirely comtained within a cell of the source grid and (2) the target
    ! cell spans at least two cells of the source grid.


    if ( j0 .EQ. j1 ) then
    ! The target cell is entirely contained within a cell of the source
    ! grid. This situation is represented by the following schematic, where
    ! the cell in which x0 and x1 are located has index j0=j1 :
    ! 
    ! ----|-----o--------o----------|-------------
    !           x0       x1
    !
      ! Determine normalized coordinates
      xi0 = x0 / grid0%h(j0) - grid0%x(j0) / grid0%h(j0)
      xi1 = x1 / grid0%h(j0) - grid0%x(j0) / grid0%h(j0)

      ! Depending on which polynomial is used, integrate quantity
      ! between x0 and xi1. Integration is carried out in normalized
      ! coordinates, hence: \int_x0^x1 p(x) dx = h \int_xi0^xi1 p(xi) dxi
      select case ( method )
    
        case ( INTEGRATION_PCM )     
          q = ppoly0%coefficients(j0,1) * ( x1 - x0 )
        
        case ( INTEGRATION_PLM )    
          q = grid0%h(j0) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 1 )
        
        case ( INTEGRATION_PPM )
          q = grid0%h(j0) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 2 )
    
        case ( INTEGRATION_PQM )
          q = grid0%h(j0) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 4 )

        case default
            print *,'The selected integration method is invalid'
            call abort()
      end select     
    
    else
    ! The target cell spans at least two cells of the source grid.
    ! This situation is represented by the following schematic, where
    ! the cells in which x0 and x1 are located have indexes j0 and j1,
    ! respectively :
    ! 
    ! ----|-----o---|--- ... --|---o----------|-------------
    !           x0                 x1
    !
    ! We first integrate from x0 up to the right boundary of cell j0, then
    ! add the integrated amounts of cells located between j0 and j1 and then
    ! integrate from the left boundary of cell j1 up to x1

      q = 0.0
      ! Integrate from x0 up to right boundary of cell j0
      !xi0 = x0 / grid0%h(j0) - grid0%x(j0) / grid0%h(j0)
      xi0 = (x0 - grid0%x(j0)) / grid0%h(j0)
      xi1 = 1.0
      select case ( method )
    
      case ( INTEGRATION_PCM )
!          PRINT *,'CALLING INTEGRATION_POLY',j0,j1
          q = q + ppoly0%coefficients(j0,1) * ( grid0%x(j0+1) - x0 )
        case ( INTEGRATION_PLM )    
          q = q + grid0%h(j0) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 1 )
    
        case ( INTEGRATION_PPM )
          q = q + grid0%h(j0) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 2 )
    
        case ( INTEGRATION_PQM )
          q = q + grid0%h(j0) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j0,:), 4 )

        case default
            print *, 'The selected integration method is invalid'
            call abort()
      end select     
    
      ! Integrate contents within cells strictly comprised between j0 and j1
      if ( j1 .GT. (j0+1) ) then
        do k = j0+1,j1-1
          q = q + grid0%h(k) * u0(k)
        end do
      end if

      ! Integrate from left boundary of cell j1 up to x1
      xi0 = 0.0
      !xi1 = x1 / grid0%h(j1) - grid0%x(j1) / grid0%h(j1)

      xi1 = (x1 - grid0%x(j1)) / grid0%h(j1)
      
      select case ( method )
    
        case ( INTEGRATION_PCM )     
          q = q + ppoly0%coefficients(j1,1) * ( x1 - grid0%x(j1) )
        
        case ( INTEGRATION_PLM )    
          q = q + grid0%h(j1) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j1,:), 1 )
    
        case ( INTEGRATION_PPM )
          q = q + grid0%h(j1) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j1,:), 2 )
    
        case ( INTEGRATION_PQM )
          q = q + grid0%h(j1) * &
              integration_polynomial ( xi0, xi1, ppoly0%coefficients(j1,:), 4 )

        case default
            print *,'The selected integration method is invalid'
            call abort()
      end select     
      
    end if ! end integration for non-vanished cells 
    
    ! The cell average is the integrated value divided by the cell width
    u1(i) = q / grid1%h(i)
    
    end if ! end if clause to check if cell is vanished
    
  end do ! end i loop on target grid cells

end subroutine remapping_integration


!------------------------------------------------------------------------------
! Inflate vanished layers to finite (nonzero) width
!------------------------------------------------------------------------------
subroutine inflate_vanished_layers ( grid, min_thickness )

  ! Argument
  type(grid1d_t), intent(inout)       :: grid
  real, intent(in) :: min_thickness
    
  ! Local variable
  integer   :: N
  integer   :: k
  integer   :: k_found
  integer   :: count_nonzero_layers
  real      :: delta
  real      :: correction
  real      :: max_thickness

  
  N = grid%nb_cells
  
  ! Count number of nonzero layers
  count_nonzero_layers = 0
  do k = 1,N
    if ( grid%h(k) .GT. min_thickness ) then
      count_nonzero_layers = count_nonzero_layers + 1
    end if
  end do

  ! If all layer thicknesses are greater than the threshold, exit routine
  if ( count_nonzero_layers .eq. N ) return

  ! If all thicknesses are zero, inflate them all and exit
  if ( count_nonzero_layers .eq. 0 ) then  
    do k = 1,N
      grid%h(k) = min_thickness
    end do
    return
  end if    
  
  ! Inflate zero layers
  correction = 0.0
  do k = 1,N
    if ( grid%h(k) .le. min_thickness ) then
      delta = min_thickness - grid%h(k)
      correction = correction + delta
      grid%h(k) = grid%h(k) + delta
    end if  
  end do
  
  ! Modify thicknesses of nonzero layers to ensure volume conservation
  max_thickness = grid%h(1)
  k_found = 1
  do k = 1,grid%nb_cells
    if ( grid%h(k) .gt. max_thickness ) then
      max_thickness = grid%h(k)
      k_found = k
    end if  
  end do
  
  grid%h(k_found) = grid%h(k_found) - correction
  
  ! Redefine grid coordinates according to new layer thicknesses
  grid%x(1) = 0.0
  do k = 1,N
    grid%x(k+1) = grid%x(k) + grid%h(k)
  end do    
  
end subroutine inflate_vanished_layers

end module MIDAS_remapping
