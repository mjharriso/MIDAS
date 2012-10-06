module remap_sfc_fluxes


  implicit none

  public remap

contains

   subroutine remap(flux,rho,rho_bound,flux_remap)

     real(kind=8), dimension(:) :: flux, rho
     real(kind=8), dimension(:) :: flux_remap
     real(kind=8), dimension(:) :: rho_bound

     integer :: i,j
     integer :: len, len_b

     len=size(flux,1)
     len_b=size(rho_bound,1)

     flux_remap=0.0
     do i=1,len
        if (rho(i) .lt. rho_bound(1)) then
            cycle
        endif
        do j=2,len_b
           if (rho(i) .ge. rho_bound(j-1) .and. rho(i) .lt. rho_bound(j)) then
               flux_remap(j-1) = flux_remap(j-1) + flux(i)
               exit
           endif
        enddo
     enddo

     return

   end subroutine remap

 end module remap_sfc_fluxes
