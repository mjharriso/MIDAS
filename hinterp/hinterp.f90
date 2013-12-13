module hinterp_mod

  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_type, horiz_interp_init,horiz_interp_del
  use fms_mod, only : fms_init, fms_end
  use fms_io_mod, only : fms_io_init,fms_io_exit
  implicit none
  private

  

  public :: hinterp


  type(horiz_interp_type), dimension(:), allocatable :: Interp
  
contains

!
! method  0=conservative; 1=bilinear; 2=bicubic
  
!  subroutine hinterp( lon_in,lat_in,mask_in,data_in,lon_out,lat_out,mask_out,&
!       data_out,src_modulo,method,missing)

  subroutine hinterp( lon_in,lat_in,data_in,lon_out,lat_out,&
                        data_out,src_modulo,method,missing,verbose)    
    
    real(kind=8), intent(in),  dimension(:,:)      :: lon_in
    real(kind=8), intent(in),  dimension(:,:) :: lat_in
!    real(kind=8), intent(in),  dimension(:,:,:,:) :: mask_in
    real(kind=8), intent(in),  dimension(:,:,:,:) :: data_in            
    real(kind=8), intent(in),  dimension(:,:)      :: lon_out
!    real(kind=8), intent(inout),  dimension(:,:)      :: mask_out    
    real(kind=8), intent(in),  dimension(:,:) :: lat_out
    real(kind=8), intent(inout),  dimension(:,:,:,:) :: data_out
    real(kind=8), intent(in)                        :: missing
    
    logical, intent(in) :: src_modulo
    integer, intent(in), optional :: verbose
    
    integer(kind=4), intent(in) :: method

!    real,  dimension(size(data_in,1),size(data_in,2)) :: mask_in_
    real(kind=8),  dimension(size(data_in,1),size(data_in,2)) :: data_in_
    real(kind=8),  dimension(size(data_out,1),size(data_out,2)) :: data_out_

    real, dimension(size(lon_in,1)) :: xax_in
    real, dimension(size(lon_in,2)) :: yax_in    
    
    integer :: ierr,nk,nt,k,m,ni,nj,i,j


    integer :: iverbose, num_nbrs
    real    :: max_dist


    iverbose=0;num_nbrs=0;max_dist=0.0

    if (PRESENT(verbose)) iverbose = verbose
    
    nk=size(data_out,3);nt=size(data_out,4);ni=size(lon_out,1);nj=size(lon_out,2)
    
    ierr=-1

    
    call fms_init()
    call fms_io_init()
    call horiz_interp_init()
!    mask_in_(:,:)=mask_in(:,:,1,1)


    allocate(Interp(nk))
    


    do m=1,nt
       do k=1,nk

          if (m.eq.1) then
              if (method == 0) then
                  call horiz_interp_new(Interp(k),lon_in,lat_in,lon_out,lat_out,interp_method="conservative")    
              else if (method == 1) then
!        call horiz_interp_new(Interp,lon_in,lat_in,lon_out,lat_out,verbose,"bilinear",mask_in=mask_in_,mask_out=mask_out)
                  call horiz_interp_new(Interp(k),lon_in,lat_in,lon_out,lat_out,verbose,"bilinear")        
              else if (method == 2) then
                  xax_in=lon_in(:,1)
                  yax_in=lat_in(1,:)
                  call horiz_interp_new(Interp(k),xax_in,yax_in,lon_out,lat_out,interp_method="bicubic")
              else
                  print *,'Invalid interpolation method passed to hinterp: ',method
                  return
              endif
          endif
          
          data_in_(:,:)=data_in(:,:,k,m)
!          mask_in_(:,:)=mask_in(:,:,k,m)
          call horiz_interp(Interp(k),data_in_, data_out_,verbose,missing_value=missing,new_missing_handle=.true.)

          do j=1,nj
             do i=1,ni
                if (abs(data_out(i,j,k,m)-missing) .lt. abs(missing)*1.e-3) then
                    data_out(i,j,k,m)=missing
                else
                    data_out(i,j,k,m)=data_out_(i,j)
                endif
             enddo
          enddo
       enddo
    enddo
    

    deallocate(Interp)
!    call horiz_interp_del(Interp)
    
    ierr=0
    

  end subroutine hinterp
  

  
  
end module hinterp_mod
