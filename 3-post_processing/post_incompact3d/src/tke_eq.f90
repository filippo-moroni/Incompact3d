!------------------------------------------------------------!
! This file contains subroutines for the calculation of      !
! Turbulent Kinetic Energy (TKE) equation for Channels       !
! and TTBLs (statistical homogeneity in x and z directions). !
!------------------------------------------------------------!

! Turbulent transport of TKE without y-derivative 
subroutine turb_trans_tke(ux2,uy2,uz2,nr,nt,kvprime2)   

  use param
  use variables
  use decomp_2d
  
  implicit none
  
  ! Fluctuations
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2
  
  ! Number of flow realizations and number of snapshots
  integer,     intent(in) :: nr, nt                                                                                                  
  integer                 :: i,j,k
  
  ! Denominator of the divisions  
  real(mytype) :: den   
  
  ! Turbulent kinetic energy multiplied by v', y-pencils
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: kvprime2   
          
  ! Different denominator according to Channel or TTBL mode  
#ifdef TTBL_MODE 
  den = real(nr,mytype)
#else
  den = real(nr*nt,mytype)
#endif
  
  ! Turbulent transport of TKE without y-derivative [ 0.5 * (u'^2 + v'^2 + w'^2) * v']
  kvprime2 = kvprime2 + zpfive*(ux2**2 + uy2**2 + uz2**2)*uy2 / den
    
end subroutine stat_vorticity






