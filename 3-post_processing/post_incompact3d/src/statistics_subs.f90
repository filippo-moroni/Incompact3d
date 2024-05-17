!----------------------------------------------------------!
!   This file contains subroutines for post-processing     !
!                 of Incompact3d Snapshots.                !
!      Adapted from original Incompact3d file (v2.0)       !
!                    of R. Corsini                         !
!----------------------------------------------------------!

! This file contains subroutines for the calculation of the following statistics:

! One-point statistics:
! - average, variance, skewness and kurtosis of velocity components
! - Reynolds stresses
! - average and variance of pressure, scalar field and mixed fluctuations
! - vorticity and mean gradient
! - total dissipation rate

! Two-points statistics:
! - (to be implemented, e.g. correlations) 
!----------------------------------------------------------!

! Mean statistics (average, variance, skewness, kurtosis)
subroutine stat_mean(ux2,uy2,uz2,pre2,phi2,nr,nt, &
                     u1mean,v1mean,w1mean,u2mean,v2mean,w2mean, &
                     u3mean,v3mean,w3mean,u4mean,v4mean,w4mean, &
                     uvmean,uwmean,vwmean,pre1mean,pre2mean,phi1mean, &
                     phi2mean,uphimean,vphimean,wphimean)

  use param
  use variables
  use decomp_2d
  use decomp_2d_io
  
  implicit none
  
  ! Variables definition
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2,pre2     ! velocity components and pressure
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3)) :: phi2                 ! scalar field
  integer,     intent(in) :: nr                                                         ! number of flow realizations
  integer,     intent(in) :: nt                                                         ! number of snapshots
  
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2                             ! temporary array (local)
  
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: u1mean,v1mean,w1mean        ! 1st order moment (average)
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: u2mean,v2mean,w2mean        ! 2nd order moment (variance)
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: u3mean,v3mean,w3mean        ! 3rd order moment (skewness)
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: u4mean,v4mean,w4mean        ! 4th order moment (kurtosis)
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: uvmean,uwmean,vwmean        ! Reynolds stresses
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: pre1mean,pre2mean           ! average and variance of pressure
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: phi1mean,phi2mean           ! average and variance of scalar field
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: uphimean,vphimean,wphimean  ! average of mixed fluctuations
  
  real(mytype) :: den  ! denominator of the divisions
  
#ifdef TTBL_MODE 
  den = real(nr,mytype)
#else
  den = real(nr*nt,mytype)
#endif
                                                                                                                                           
  !---x-component---!
  ! average
  u1mean=u1mean+ux2/den 

  ! variance
  ta2=ux2*ux2
  u2mean=u2mean+ta2/den

  ! skewness
  ta2=ux2*ux2*ux2
  u3mean=u3mean+ta2/den

  ! kurtosis 
  ta2=ux2*ux2*ux2*ux2
  u4mean=u4mean+ta2/den
  
  !---y-component---!
  ! average
  v1mean=v1mean+uy2/den 

  ! variance
  ta2=uy2*uy2
  v2mean=v2mean+ta2/den

  ! skewness
  ta2=uy2*uy2*uy2
  v3mean=v3mean+ta2/den

  ! kurtosis 
  ta2=uy2*uy2*uy2*uy2
  v4mean=v4mean+ta2/den
  
  !---z-component---!
  ! average
  w1mean=w1mean+uz2/den

  ! variance
  ta2=uz2*uz2
  w2mean=w2mean+ta2/den

  ! skewness
  ta2=uz2*uz2*uz2
  w3mean=w3mean+ta2/den

  ! kurtosis 
  ta2=uz2*uz2*uz2*uz2
  w4mean=w4mean+ta2/den


  !---Reynolds stresses---!
  !<uv>
  ta2=ux2*uy2
  uvmean=uvmean+ta2/den

  !<uw>
  ta2=ux2*uz2
  uwmean=uwmean+ta2/den

  !<vw>
  ta2=uy2*uz2
  vwmean=vwmean+ta2/den

  !---pressure---!
  pre1mean=pre1mean+pre2/den

  ta2=pre2*pre2
  pre2mean=pre2mean+ta2/den

  !---scalar---!
  if (iscalar==1) then

        ! average: phi
        phi1mean=phi1mean+phi2/den

        ! variance: phi
        ta2=phi2*phi2
        phi2mean=phi2mean+ta2/den

        ! mixed fluctuations: u',phi'
        ta2=ux2*phi2
        uphimean=uphimean+ta2/den

        ! mixed fluctuations: v',phi'
        ta2=uy2*phi2
        vphimean=vphimean+ta2/den

        ! mixed fluctuations: w',phi'
        ta2=uz2*phi2
        wphimean=wphimean+ta2/den

  endif

end subroutine stat_mean
!********************************************************************
! Vorticity and mean gradient sqrt[ (du/dy)**2 + (dw/dy)**2 ]
subroutine stat_vorticity(ux1,uy1,uz1,nr,nt,vortxmean2,vortymean2,vortzmean2,mean_gradient2)   

  use param
  use variables
  use decomp_2d
  use decomp_2d_io
  use MPI
  use dbg_schemes, only : sqrt_prec
  
  implicit none
  
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
  integer,     intent(in) :: nr                                                 ! number of flow realizations
  integer,     intent(in) :: nt                                                 ! number of snapshots
  integer                 :: i,j,k
  
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,di2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,di3
  
  real(mytype) :: den  ! denominator of the divisions 
  real(mytype) :: lind
  
  ! Vorticity 
  real(mytype),            dimension(xsize(1),xsize(2),xsize(3)) :: vortxmean1,vortymean1,vortzmean1  ! average vorticity components, x-pencils
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: vortxmean2,vortymean2,vortzmean2  ! average vorticity components, y-pencils
  
  ! Mean gradient
  real(mytype),            dimension(xsize(1),xsize(2),xsize(3)) :: mean_gradient1                    ! mean gradient, x-pencils
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: mean_gradient2                    ! mean gradient, y-pencils
      
  ! x-derivatives
  call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,lind)
  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,lind)
  call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,lind)
  ! y-derivatives
  call transpose_x_to_y(ux1,td2)
  call transpose_x_to_y(uy1,te2)
  call transpose_x_to_y(uz1,tf2)
  call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,lind)
  call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,lind)
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
  
  di1 = zero
  
#ifdef TTBL_MODE 
  den = real(nr,mytype)
#else
  den = real(nr*nt,mytype)
#endif
  
  !---- Mean gradient ----!
  do i=1,xsize(1)
    do j=1,xsize(2)
      do k=1,xsize(3)
        di1(i,j,k) = sqrt_prec(td1(i,j,k)**2 + tf1(i,j,k)**2) ! sqrt[ (du/dy)**2 + (dw/dy)**2 ]
      enddo
    enddo
  enddo
  
  mean_gradient1 = mean_gradient1 + di1/den
  
  !---Vorticity average---!
  
  ! Vorticity along x 
  di1 = tf1 - th1  !dw/dy - dv/dz
  vortxmean1 = vortxmean1 + di1/den
    
  ! Vorticity along y
  di1 = tg1 - tc1  !du/dz - dw/dx
  vortymean1 = vortymean1 + di1/den
     
  ! Vorticity along z
  di1 = tb1 - td1  !dv/dx - du/dy
  vortzmean1 = vortzmean1 + di1/den
   
  ! Transpose arrays along y
  call transpose_x_to_y(vortxmean1,vortxmean2)
  call transpose_x_to_y(vortymean1,vortymean2)
  call transpose_x_to_y(vortzmean1,vortzmean2)
  call transpose_x_to_y(mean_gradient1,mean_gradient2)

end subroutine stat_vorticity

!********************************************************************
! Calculate total dissipation
subroutine stat_dissipation(ux1,uy1,uz1,nr,nt,epsmean2)
  
  use param
  use variables
  use decomp_2d
  use decomp_2d_io
  use MPI
  use dbg_schemes, only : sqrt_prec
  
  implicit none
  
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
  integer,     intent(in) :: nr                                                 ! number of flow realizations
  integer,     intent(in) :: nt                                                 ! number of snapshots
  integer                 :: i,j,k
  
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,di2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,di3
  
  real(mytype) :: den  ! denominator of the divisions 
  real(mytype) :: lind
  
  ! Pseudo-dissipation 
  real(mytype),            dimension(xsize(1),xsize(2),xsize(3)) :: epsmean1  ! average total dissipation, x-pencils
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: epsmean2  ! average total dissipation, y-pencils
  
  ! x-derivatives
  call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,lind)
  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,lind)
  call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,lind)
  ! y-derivatives
  call transpose_x_to_y(ux1,td2)
  call transpose_x_to_y(uy1,te2)
  call transpose_x_to_y(uz1,tf2)
  call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,lind)
  call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,lind)
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

  di1 = zero
  
#ifdef TTBL_MODE 
  den = real(nr,mytype)
#else
  den = real(nr*nt,mytype)
#endif
  
  !---- Instantaneous total dissipation rate ----!
  
  di1 = (2*ta1**2 + 2*te1**2 + 2*ti1**2 + &  !(2*du/dx**2 + 2*dv/dy**2 + 2*dw/dz**2 +
         (tb1+td1)**2 +                   &  !(dv/dx+du/dy)**2 +
         (tf1+th1)**2 +                   &  !(dw/dy+dv/dz)**2 +
         (tg1+tc1)**2)*xnu                   !(du/dz+dw/dx)**2)*nu
  
  epsmean1 = epsmean1 + di1/den
  
  ! Transpose array along y
  call transpose_x_to_y(epsmean1,epsmean2)

end subroutine stat_dissipation



