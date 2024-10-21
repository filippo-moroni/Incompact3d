
!This file is not part of standard Xcompact3d releases (xcompact3d.com).

!-----------------------------------------------------------------------------!
! DESCRIPTION: This file contains subroutines for the calculation of the 
!              following statistics:
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
! DESCRIPTION: 
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
!-----------------------------------------------------------------------------!

! to be completed 


subroutine stat_correlation_z(ux2,uy2,uz2,nx,nz,nr,nt,...)

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use post_processing, only : under_sampling
  
  implicit none
   
  ! Variables definition, velocity, (y-pencils)
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2          
  
  ! Number of points in homogeneous directions, number of snapshots and number of realizations
  integer,     intent(in) :: nx,nz,nt,nr
  
  ! Local work arrays
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3,ta3
  
  ! Correlation functions (first index: j (rows); second index: r (columns))
  real(mytype),intent(inout),dimension(zsize(2),zsize(3)) :: RuuzH1, RvvzH1, RwwzH1, RuvzH1, RppzH1
  
  real(mytype) :: den          ! denominator of the divisions
  integer      :: i,j,k,rr,kpr 

#ifdef TTBL_MODE 
  den = real(nx*nz*nr/under_sampling**2,mytype)
#else
  den = real(nx*nz*nr*nt/under_sampling**2,mytype)
#endif

  ! Transpose arrays along z
  call transpose_y_to_z(ux2,ux3)
  call transpose_y_to_z(uy2,uy3)
  call transpose_y_to_z(uz2,uz3)
  
  ! Correlation function calculation
  do k=1,zsize(3)
      do j=1,zsize(2)
          do i=1,zsize(1)
              do rr=1,zsize(3)
                  
                  ! Index for z-direction plus index for separation variable 'r'
                  kpr = k + rr - 1
                  
                  ! Shift to the beginning of the array if we go beyond its index range (periodic)
                  if (kpr > nz) kpr = kpr - nz
                  
                  !--- Streamwise fluctuations correlation ---!
                  
                  ! Product of fluctuations at distance 'r'
                  ta3(i,j,k) = ux3(i,j,k)*ux3(i,j,kpr)
                  
                  ! Accumulation inside the correlation function variable (at each subdomain)
                  RuuzH1(j,rr) = RuuzH1(j,rr) + ta3(i,j,k)/den
                                    
              enddo
          enddo
      enddo
  enddo

end subroutine stat_correlation_z




