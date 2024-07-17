!------------------------------------------------------------!
! This file contains subroutines for the calculation of      !
! Turbulent Kinetic Energy (TKE) equation for Channels       !
! and TTBLs (statistical homogeneity in x and z directions). !
!------------------------------------------------------------!

! Calculation of TKE terms that cannot be obtained from mean flow statistics data.
! All results obtained here are in y-pencils.
! Results must be later derived.

subroutine extra_terms_tke(ux2,uy2,uz2,nr,nt,kvprime_mean,pprimevprime_mean,)   

  use param
  use variables
  use decomp_2d
  
  implicit none
  
  ! Fluctuations
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2,pre2
  
  ! Number of flow realizations and number of snapshots
  integer,     intent(in) :: nr, nt                                                                                                  
  integer                 :: i,j,k
  
  ! Arrays for derivative calculations
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,di2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,di3
  
  ! Denominator of the divisions  
  real(mytype) :: den   
  real(mytype) :: lind
  
  ! Turbulent transport of TKE by v'
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: kvprime_mean
  
  ! Pressure-strain (or coupling) term, y-direction
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: pprimevprime_mean
  
  ! Pseudo-dissipation for TKE
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: pseudo_eps_tke_mean   
          
  ! Different denominator according to Channel or TTBL mode  
#ifdef TTBL_MODE 
  den = real(nr,mytype)
#else
  den = real(nr*nt,mytype)
#endif
  
  ! Turbulent transport of TKE by v' [ 0.5 * (u'^2 + v'^2 + w'^2) * v']
  kvprime_mean = kvprime_mean + zpfive*(ux2**2 + uy2**2 + uz2**2)*uy2 / den
    
  ! Pressure-strain (or coupling) term, y-direction
  pprimevprime_mean = pprimevprime_mean + pre2*uy2 / den
  
  ! Transpose fluctuations in x-direction
  call transpose_y_to_x(ux2,ux1)
  call transpose_y_to_x(uy2,uy1)
  call transpose_y_to_x(uz2,uz1)
        
  ! x-derivatives
  call derx (ta1,ux1,di1,sx,ffx, fsx, fwx, xsize(1),xsize(2),xsize(3),0,lind)
  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,lind)
  call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,lind)
  
  ! y-derivatives
  call transpose_x_to_y(ux1,td2)
  call transpose_x_to_y(uy1,te2)
  call transpose_x_to_y(uz1,tf2)
  call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,lind)
  call dery (tb2,te2,di2,sy,ffy, fsy, fwy, ppy,ysize(1),ysize(2),ysize(3),0,lind)
  call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,lind)
  
  ! z-derivatives
  call transpose_y_to_z(td2,td3)
  call transpose_y_to_z(te2,te3)
  call transpose_y_to_z(tf2,tf3)
  call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,lind)
  call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,lind)
  call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,lind)
  
  ! all back to x-pencils
  call transpose_z_to_y(ta3,td2)
  call transpose_z_to_y(tb3,te2)
  call transpose_z_to_y(tc3,tf2)
  
  call transpose_y_to_x(td2,tg1)
  call transpose_y_to_x(te2,th1)
  call transpose_y_to_x(tf2,ti1)
  call transpose_y_to_x(ta2,td1)
  call transpose_y_to_x(tb2,te1)
  call transpose_y_to_x(tc2,tf1)
  
  !du/dx=ta1 du/dy=td1 and du/dz=tg1
  !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
  !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
  
  ! Pseudo-dissipation for TKE 
  di1 = xnu * (ta1**2 + td1**2 + tg1**2 + &
               tb1**2 + te1**2 + th1**2 + &
               tc1**2 + tf1**2 + ti1**2 )
  
  ! Transpose array along y and sum
  call transpose_x_to_y(di1,di2)
  vortxmean2 = vortxmean2 + di2/den
  
  
  pseudo_eps_tke_mean 
    
end subroutine extra_terms_tke






