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

!----------------------------------------------------------!

! Mean statistics (average, variance, skewness, kurtosis)
subroutine STAT_MEAN(ux2,uy2,uz2,pre2,phi2,ta2, &
                     u1mean,v1mean,w1mean,u2mean,v2mean,w2mean, &
                     u3mean,v3mean,w3mean,u4mean,v4mean,w4mean, &
                     uvmean,uwmean,vwmean,pre1mean,pre2mean,phi1mean, &
                     phi2mean,uphimean,vphimean,wphimean,nr)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io
  
  implicit none
  
  ! Variables definition
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2,pre2     ! velocity components and pressure
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3)) :: phi2                 ! scalar field
  
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2                             ! temporary array
  
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u1mean,v1mean,w1mean            ! 1st order moment (average)
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u2mean,v2mean,w2mean            ! 2nd order moment (variance)
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u3mean,v3mean,w3mean            ! 3rd order moment (skewness)
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: u4mean,v4mean,w4mean            ! 4th order moment (kurtosis)
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uvmean,uwmean,vwmean            ! Reynolds stresses
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: pre1mean,pre2mean               ! average and variance of pressure
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: phi1mean,phi2mean               ! average and variance of scalar field
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uphimean,vphimean,wphimean      ! average of mixed fluctuations
  
  integer :: is                                                                         ! index for the different scalar fields
  integer :: nr
                                                                                                                                            
  !---x-component---!
  ! average
  u1mean=u1mean+ux2/nr 

  ! variance
  ta2=ux2*ux2
  u2mean=u2mean+ta2/nr

  ! skewness
  ta2=ux2*ux2*ux2
  u3mean=u3mean+ta2/nr

  ! kurtosis 
  ta2=ux2*ux2*ux2*ux2
  u4mean=u4mean+ta2/nr
  
  !---y-component---!
  ! average
  v1mean=v1mean+uy2/nr 

  ! variance
  ta2=uy2*uy2
  v2mean=v2mean+ta2/nr

  ! skewness
  ta2=uy2*uy2*uy2
  v3mean=v3mean+ta2/nr

  ! kurtosis 
  ta2=uy2*uy2*uy2*uy2
  v4mean=v4mean+ta2/nr
  
  !---z-component---!
  ! average
  w1mean=w1mean+uz2/nr 

  ! variance
  ta2=uz2*uz2
  w2mean=w2mean+ta2/nr

  ! skewness
  ta2=uz2*uz2*uz2
  w3mean=w3mean+ta2/nr

  ! kurtosis 
  ta2=uz2*uz2*uz2*uz2
  w4mean=w4mean+ta2/nr


  !---Reynolds stresses---!
  !<uv>
  ta2=ux2*uy2
  uvmean=uvmean+ta2/nr

  !<uw>
  ta2=ux2*uz2
  uwmean=uwmean+ta2/nr

  !<vw>
  ta2=uy2*uz2
  vwmean=vwmean+ta2/nr

  !---pressure---!
  pre1mean=pre1mean+pre2/nr

  ta2=pre2*pre2
  pre2mean=pre2mean+ta2/nr

  !---scalar---!
  if (iscalar==1) then

        ! average: phi
        phi1mean(:,:,:)=phi1mean(:,:,:)+phi2(:,:,:)/nr

        ! variance: phi
        ta2=phi2(:,:,:)*phi2(:,:,:)
        phi2mean(:,:,:)=phi2mean(:,:,:)+ta2/nr

        ! mixed fluctuations: u',phi'
        ta2=ux2*phi2(:,:,:)
        uphimean(:,:,:)=uphimean(:,:,:)+ta2/nr

        ! mixed fluctuations: v',phi'
        ta2=uy2*phi2(:,:,:)
        vphimean(:,:,:)=vphimean(:,:,:)+ta2/nr

        ! mixed fluctuations: w',phi'
        ta2=uz2*phi2(:,:,:)
        wphimean(:,:,:)=wphimean(:,:,:)+ta2/nr

     enddo
  endif

end subroutine STAT_MEAN
!********************************************************************
! Vorticity 
subroutine STAT_VORTICITY(ux1,uy1,uz1,ifile,nr)   ! to be checked and completed

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io
  USE MPI
  
  implicit none

  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,di2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,di3
  integer :: ijk,nvect1,ifile
  character(len=30) :: filename
  
  integer :: nr
  real(mytype) :: lind
  
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: vortxmean,vortymean,vortzmean   ! average vorticity components
      
  nvect1=xsize(1)*xsize(2)*xsize(3)

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
  
  !---Vorticity average---!
  
  ! Variable declaration
  di1=0._mytype
    
  ! Vorticity along x    
  do ijk=1,nvect1 !dw/dy - dv/dz
     di1(ijk,1,1)=tf1(ijk,1,1)-th1(ijk,1,1)
  enddo  
  vortxmean = vortxmean + di1/nr
  
  ! Vorticity along y
  do ijk=1,nvect1 !du/dz - dw/dx
     di1(ijk,1,1)=tg1(ijk,1,1)-tc1(ijk,1,1)
  enddo
  vortymean = vortymean + di1/nr
  
  ! Vorticity along z
   do ijk=1,nvect1 !dv/dx - du/dy
      di1(ijk,1,1)=tb1(ijk,1,1)-td1(ijk,1,1)
   enddo
   vortzmean = vortzmean + di1/nr

end subroutine STAT_VORTICITY
!********************************************************************



