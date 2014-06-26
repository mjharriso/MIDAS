module pyale_mod

  use regrid_defs
  use regrid_grid1d_class
  use regrid_ppoly_class
  use regrid_edge_values
  use regrid_edge_slopes  
  use regrid_pcm                                     ! see 'regrid_pcm.F90'
  use regrid_plm                                     ! see 'regrid_pcm.F90'
  use regrid_ppm                                     ! see 'regrid_pcm.F90'
  use regrid_pqm                                     ! see 'regrid_pcm.F90'
  
  use MIDAS_remapping

  implicit none

  private


  integer, parameter :: max_grids=10
  type(grid1d_t), dimension(0:max_grids), target :: grid
  
  integer :: ngrids = 0

  public :: pyale_grid_init, pyale_grid_set, remap, pyale_grid_destroy
  
contains


  function pyale_grid_init(nb_cells)

    integer, intent(in) :: nb_cells
    integer :: pyale_grid_init
    type(grid1d_t), pointer :: grid_ptr

    grid_ptr=>grid(ngrids)
    call grid1d_init(grid_ptr,nb_cells,ngrids)

    pyale_grid_init=ngrids

    ngrids=ngrids+1

    return

  end function pyale_grid_init

  subroutine pyale_grid_destroy()

    integer :: n
    type(grid1d_t), pointer :: grid_ptr
    
    do n=0,ngrids
       grid_ptr => grid(n)
       call grid1d_destroy(grid_ptr)
    enddo

    ngrids=0

  end subroutine pyale_grid_destroy
  
    
  subroutine pyale_grid_set(n, x)
    integer, intent(in) :: n
    real(kind=8), dimension(:), intent(in) :: x


    integer :: nz, k

    nz=size(x,1)-1

    if (nz /= grid(n)%nb_cells) then
        print *,'improper size array in call to grid_set'
        print *, 'ngrids= ',ngrids
        print *,'nz= ',nz, ' n= ',n,' nb_cells= ',grid(n)%nb_cells
        call abort
    endif        

    do k=1,nz
       grid(n)%x(k)=x(k)
       grid(n)%h(k)=x(k+1)-x(k)
    enddo

    grid(n)%x(nz+1) = x(nz+1)

    return

  end subroutine pyale_grid_set


  subroutine remap(u0,u1,zi,zo,method,bndy_extrapolation,missing)
    real(kind=8), dimension(:,:,:), intent(in) :: u0
    real(kind=8), dimension(:,:,:), intent(inout) :: u1
    real(kind=8), dimension(:,:,:), intent(in) :: zi,zo
    logical, intent(in) :: bndy_extrapolation
    real(kind=8), intent(in), optional :: missing
    
    character(len=*), intent(in) :: method
    
    type(ppoly_t) :: ppoly

    real, parameter :: epsln=1.e-10
    real, parameter :: min_thickness=1.e-9
    
    integer :: nz,ni,nj,i,j,k,imethod, degree,nz2, n1, n2
    real, dimension(size(u0,3))    :: uin
    
    ni=size(u0,1);nj=size(u0,2);nz=size(u0,3);nz2=size(u1,3)

    
    select case (method)
    case ('pcm')
        imethod=INTEGRATION_PCM
        degree=0
    case('plm')
        imethod=INTEGRATION_PLM
        degree=1
    case('ppm_h4','ppm_ih4')
        imethod=INTEGRATION_PPM
        degree=2
    case('pqm_ih4ih3','pqm_ih6ih5')
        imethod=INTEGRATION_PQM
        degree=4
    case default
        print *,'Invalid integration method'
        call abort()
    end select

    n1= pyale_grid_init(nz)
    n2= pyale_grid_init(nz2)    

    call ppoly_init(ppoly,nz,degree)

    do j=1,nj
       do i=1,ni
!          if (abs(zi(i,j,nz)-epsln)>1.e-2) then
              call pyale_grid_set(n1,zi(i,j,:))
              call pyale_grid_set(n2,zo(i,j,:))
!              call inflate_vanished_layers(grid(n1),min_thickness)
!              call inflate_vanished_layers(grid(n2),min_thickness)              

              uin(:)=u0(i,j,:)

              if (PRESENT(missing)) then
                  do k=2,nz-1
                     if (abs((uin(k)-missing)/missing) < 1.e-3) then
                         uin(k)=uin(k-1)
                     endif
                  enddo

              endif

              select case (method)
              case ('pcm')
                  call pcm_reconstruction(grid(n1),ppoly,uin)
              case('plm')
                  call plm_reconstruction(grid(n1),ppoly,uin)
                  if (bndy_extrapolation) call plm_boundary_extrapolation(grid(n1),ppoly,uin)
              case('ppm_h4')
!                  print *,grid
                  call edge_values_explicit_h4(grid(n1),uin,ppoly%E)
                  call ppm_reconstruction(grid(n1),ppoly,uin)
                  if (bndy_extrapolation) call ppm_boundary_extrapolation(grid(n1),ppoly,uin)
              case('ppm_ih4')
                  call edge_values_implicit_h4(grid(n1),uin,ppoly%E)
                  call ppm_reconstruction(grid(n1),ppoly,uin)
                  if (bndy_extrapolation) call ppm_boundary_extrapolation(grid(n1),ppoly,uin)
              case('pqm_ih4ih3')
                  call edge_values_implicit_h4(grid(n1),uin,ppoly%E)
                  call edge_slopes_implicit_h3(grid(n1),uin,ppoly%S)                  
                  call pqm_reconstruction(grid(n1),ppoly,uin)
                  if (bndy_extrapolation) call pqm_boundary_extrapolation_v1(grid(n1),ppoly,uin)
              case('pqm_ih6ih5')
                  call edge_values_implicit_h6(grid(n1),uin,ppoly%E)
                  call edge_slopes_implicit_h5(grid(n1),uin,ppoly%S)                  
                  call pqm_reconstruction(grid(n1),ppoly,uin)
                  if (bndy_extrapolation) call pqm_boundary_extrapolation_v1(grid(n1),ppoly,uin)                  
              case default
                  print *,'pyale error, the selected remapping method is invalid'
                  call abort()
              end select
              
              call remapping_integration(grid(n1),uin,ppoly,grid(n2),u1(i,j,:),imethod)

!          endif
       enddo
    enddo

    call pyale_grid_destroy()
    call ppoly_destroy(ppoly)
    
  end subroutine remap
  

  
end module pyale_mod
