module hinterp_mod

  use horiz_interp_mod, only : horiz_interp_new, horiz_interp, horiz_interp_type, horiz_interp_init,horiz_interp_del
  use fms_mod, only : fms_init, fms_end
  use fms_io_mod, only : fms_io_init,fms_io_exit
  implicit none
  private

  

  public :: hinterp


  type(horiz_interp_type) :: Interp
  
contains


  subroutine hinterp( lon_in,lat_in,mask_in,data_in,lon_out,lat_out,data_out,src_modulo,method,missing) 
    real(kind=8), intent(in),  dimension(:,:)      :: lon_in
    real(kind=8), intent(in),  dimension(size(lon_in,1),size(lon_in,2)) :: lat_in
    real(kind=8), intent(in),  dimension(:,:,:,:) :: mask_in
    real(kind=8), intent(in),  dimension(:,:,:,:) :: data_in            
    real(kind=8), intent(in),  dimension(:,:)      :: lon_out
    real(kind=8), intent(in),  dimension(size(lon_out,1),size(lon_out,2)) :: lat_out
    real(kind=8), intent(inout),  dimension(:,:,:,:) :: data_out
    real(kind=8), intent(in)                        :: missing
    
    logical, intent(in) :: src_modulo
    character(len=*), intent(in) :: method

    real(kind=8),  dimension(size(mask_in,1),size(mask_in,2)) :: data_in_, mask_in_
    real(kind=8),  dimension(size(data_out,1),size(data_out,2)) :: mask_out, data_out_

    real(kind=8), dimension(size(lon_in,1)) :: xax_in
    real(kind=8), dimension(size(lon_in,2)) :: yax_in    
    
    integer :: ierr,nk,nt,k,m,ni,nj,i,j


    integer(kind=4) :: verbose, num_nbrs
    real(kind=8)    :: max_dist

    verbose=0;num_nbrs=0;max_dist=0.0
    
    nk=size(data_out,3);nt=size(data_out,4);ni=size(lon_out,1);nj=size(lon_out,2)
    
    ierr=-1


    call fms_init()
    call fms_io_init()
    call horiz_interp_init()
    mask_in_(:,:)=mask_in(:,:,1,1)

    if (method == 'conservative') then
        call horiz_interp_new(Interp,lon_in,lat_in,lon_out,lat_out,interp_method="conservative")    
    else if (method == 'bilinear') then
        call horiz_interp_new(Interp,lon_in,lat_in,lon_out,lat_out,verbose,method,num_nbrs,max_dist,src_modulo,mask_in_,mask_out)
    else if (method == 'bicubic') then
        xax_in=lon_in(:,1)
        yax_in=lat_in(1,:)
        call horiz_interp_new(Interp,xax_in,yax_in,lon_out,lat_out,verbose,method,num_nbrs,max_dist,src_modulo,.true.,mask_in_,mask_out,.false.)
    else
        print *,'Invalid interpolation method passed to hinterp: ',trim(method)
        return
    endif
    
    do m=1,nt
       do k=1,nk
          data_in_(:,:)=data_in(:,:,k,m)
          mask_in_(:,:)=mask_in(:,:,k,m)
          call horiz_interp(Interp,data_in_, data_out_,verbose,mask_in_,mask_out,missing)

          do j=1,nj
             do i=1,ni
                if (mask_out(i,j).eq.0.0) then
                    data_out(i,j,k,m)=missing
                else
                    data_out(i,j,k,m)=data_out_(i,j)
                endif
             enddo
          enddo
       enddo
    enddo
    
    
    call horiz_interp_del(Interp)
    
    ierr=0
    

  end subroutine hinterp
  

  
  
end module hinterp_mod
