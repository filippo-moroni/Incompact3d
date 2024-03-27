!----------------------------------------------------------!
! This file contains subroutines for post-processing       !
!                 of Incompact3d Snapshots.                !
!      Adapted from original Incompact3d file (v2.0)       !
!                    of R. Corsini                         !
!----------------------------------------------------------!

! This file contains subroutines for the calculation of the following statistics:

! One-point statistics:
! - average, variance, skewness and kurtosis of velocity components
! - Reynolds stresses
! - average and variance of pressure, scalar field and mixed fluctuations
! - vorticity

! Two-points statistics:
! - (to be implemented, e.g. correlations) 

! Flow parameters: 
! - delta_99, displacement thickness, momentum thickness
! - shear velocity
! - related Re numbers
!----------------------------------------------------------!

! Mean statistics (average, variance, skewness, kurtosis)
subroutine stat_mean(ux2,uy2,uz2,pre2,phi2,nr, &
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
  
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2                             ! temporary array (local)
  
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: u1mean,v1mean,w1mean        ! 1st order moment (average)
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: u2mean,v2mean,w2mean        ! 2nd order moment (variance)
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: u3mean,v3mean,w3mean        ! 3rd order moment (skewness)
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: u4mean,v4mean,w4mean        ! 4th order moment (kurtosis)
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: uvmean,uwmean,vwmean        ! Reynolds stresses
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: pre1mean,pre2mean           ! average and variance of pressure
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: phi1mean,phi2mean           ! average and variance of scalar field
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: uphimean,vphimean,wphimean  ! average of mixed fluctuations
  
                                                                                                                                            
  !---x-component---!
  ! average
  u1mean=u1mean+ux2/real(nr,mytype) 

  ! variance
  ta2=ux2*ux2
  u2mean=u2mean+ta2/real(nr,mytype)

  ! skewness
  ta2=ux2*ux2*ux2
  u3mean=u3mean+ta2/real(nr,mytype)

  ! kurtosis 
  ta2=ux2*ux2*ux2*ux2
  u4mean=u4mean+ta2/real(nr,mytype)
  
  !---y-component---!
  ! average
  v1mean=v1mean+uy2/real(nr,mytype) 

  ! variance
  ta2=uy2*uy2
  v2mean=v2mean+ta2/real(nr,mytype)

  ! skewness
  ta2=uy2*uy2*uy2
  v3mean=v3mean+ta2/real(nr,mytype)

  ! kurtosis 
  ta2=uy2*uy2*uy2*uy2
  v4mean=v4mean+ta2/real(nr,mytype)
  
  !---z-component---!
  ! average
  w1mean=w1mean+uz2/real(nr,mytype) 

  ! variance
  ta2=uz2*uz2
  w2mean=w2mean+ta2/real(nr,mytype)

  ! skewness
  ta2=uz2*uz2*uz2
  w3mean=w3mean+ta2/real(nr,mytype)

  ! kurtosis 
  ta2=uz2*uz2*uz2*uz2
  w4mean=w4mean+ta2/real(nr,mytype)


  !---Reynolds stresses---!
  !<uv>
  ta2=ux2*uy2
  uvmean=uvmean+ta2/real(nr,mytype)

  !<uw>
  ta2=ux2*uz2
  uwmean=uwmean+ta2/real(nr,mytype)

  !<vw>
  ta2=uy2*uz2
  vwmean=vwmean+ta2/real(nr,mytype)

  !---pressure---!
  pre1mean=pre1mean+pre2/real(nr,mytype)

  ta2=pre2*pre2
  pre2mean=pre2mean+ta2/real(nr,mytype)

  !---scalar---!
  if (iscalar==1) then

        ! average: phi
        phi1mean=phi1mean+phi2/real(nr,mytype)

        ! variance: phi
        ta2=phi2*phi2
        phi2mean=phi2mean+ta2/real(nr,mytype)

        ! mixed fluctuations: u',phi'
        ta2=ux2*phi2
        uphimean=uphimean+ta2/real(nr,mytype)

        ! mixed fluctuations: v',phi'
        ta2=uy2*phi2
        vphimean=vphimean+ta2/real(nr,mytype)

        ! mixed fluctuations: w',phi'
        ta2=uz2*phi2
        wphimean=wphimean+ta2/real(nr,mytype)

  endif

end subroutine stat_mean
!********************************************************************
! Vorticity 
subroutine stat_vorticity(ux1,uy1,uz1,nr,vortxmean2,vortymean2,vortzmean2,mean_gradient2)   

  use param
  use variables
  use decomp_2d
  use decomp_2d_io
  use MPI
  
  implicit none
  
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
  integer,intent(in) :: nr
  
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,di2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,di3
  
  real(mytype) :: lind
  
  ! Vorticity 
  real(mytype),            dimension(xsize(1),xsize(2),xsize(3)) :: vortxmean1,vortymean1,vortzmean1  ! average vorticity components, x-pencils
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: vortxmean2,vortymean2,vortzmean2  ! average vorticity components, y-pencils
  
  ! Mean gradient
  real(mytype),            dimension(xsize(1),xsize(2),xsize(3)) :: mean_gradient1                    ! mean gradient dU/dy, x-pencils
  real(mytype),intent(out),dimension(ysize(1),ysize(2),ysize(3)) :: mean_gradient2                    ! mean gradient dU/dy, y-pencils
      
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
  
  !---- Mean gradient ----!
  
  di1 = td1  ! du/dy
  mean_gradient1 = mean_gradient1 + di1/real(nr,mytype)
  
  !---Vorticity average---!
  
  ! Vorticity along x 
  di1 = tf1 - th1  !dw/dy - dv/dz
  vortxmean1 = vortxmean1 + di1/real(nr,mytype)
    
  ! Vorticity along y
  di1 = tg1 - tc1  !du/dz - dw/dx
  vortymean1 = vortymean1 + di1/real(nr,mytype)
     
  ! Vorticity along z
  di1 = tb1 - td1  !dv/dx - du/dy
  vortzmean1 = vortzmean1 + di1/real(nr,mytype)
   
  ! Transpose arrays along y
  call transpose_x_to_y(vortxmean1,vortxmean2)
  call transpose_x_to_y(vortymean1,vortymean2)
  call transpose_x_to_y(vortzmean1,vortzmean2)
  call transpose_x_to_y(mean_gradient1,mean_gradient2)

end subroutine stat_vorticity
!********************************************************************
! Calculation of flow parameters: delta_99, displacement thickness, momentum thickness, 
! shear velocity and related Re numbers at low order only for a first quick estimation
subroutine stat_parameters(u1meanHT,ie,nt,delta_99,disp_t,mom_t,re_tau,re_ds,re_theta,sh_vel)

  use param
  use variables
  use decomp_2d
  use decomp_2d_io

  implicit none
  
  real(mytype),intent(in),dimension(ysize(2)) :: u1meanHT
  integer,     intent(in) :: ie,nt
    
  real(mytype),intent(out),dimension(nt)      :: delta_99,disp_t,mom_t
  real(mytype),intent(out),dimension(nt)      :: re_tau,re_ds,re_theta
  real(mytype),intent(out),dimension(nt)      :: sh_vel
  
  integer            :: j, iunit
  character(99)      :: filename
  
  ! Reading the coordinates in y of faces' elements
  write(filename,"('yp.dat')") 
               
  filename = adjustl(filename)
        
  open(newunit=iunit,file=trim(filename),form='formatted')
     
  do j = ystart(2),yend(2)
     
  read(iunit, *) yp(j)
    
  end do
     
  close(iunit)
     
  ! First-order accurate calculations of integral quantities
     
  ! delta_99
  do j = ystart(2),yend(2)
     
  delta_99(ie) = yp(j)  
        
  if(u1meanHT(j) < 0.01*u1meanHT(ystart(2))) exit  ! This condition is valid only for a temporal TBL (with streamwise translating wall)
               
  end do
     
  ! Modify the index in u1meanHT for:
  ! j    :  forward rectangular integration
  ! j + 1: backward rectangular integration
     
  ! displacement thickness, (O(1))
  do j = ystart(2),yend(2) - 1
     
  disp_t(ie) = disp_t(ie) + u1meanHT(j)*(yp(j+1) - yp(j)) 
                    
  end do
  
  disp_t(ie) = disp_t(ie)/uwall
     
  !disp_t(ie) = yp(ysize(2)) - disp_t(ie)  ! valid for a standard spatial BL (otherwise no further calculation for temporal BLs)
         
  ! momentum thickness, (O(1))
  do j = ystart(2),yend(2) - 1
     
  mom_t(ie) = mom_t(ie) + (u1meanHT(j)/uwall - (u1meanHT(j)/uwall)**2)*(yp(j+1) - yp(j))
                    
  end do
     
  ! shear velocity, (O(2)) 
  sh_vel(ie) = sqrt(xnu*u1meanHT(2)/yp(2))
     
  ! Reynolds numbers
  re_tau  (ie) = delta_99(ie)*sh_vel(ie)/xnu  ! friction Re number (or delta99^+)
  re_ds   (ie) = disp_t  (ie)*uwall/xnu       ! Re number based on displacement thickness delta star (ds)
  re_theta(ie) = mom_t   (ie)*uwall/xnu       ! Re number based on momentum thickness theta 
  
end subroutine stat_parameters
