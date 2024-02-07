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
! - (to be implemented, correlations) 

!----------------------------------------------------------!

! Mean statistics (average, variance, skewness, kurtosis)
subroutine STAT_MEAN(ux1,uy1,uz1,pre1,phi1,ta1, &
                     u1mean,v1mean,w1mean,u2mean,v2mean,w2mean, &
                     u3mean,v3mean,w3mean,u4mean,v4mean,w4mean, &
                     uvmean,uwmean,vwmean,pre1mean,pre2mean,phi1mean, &
                     phi2mean,uphimean,vphimean,wphimean,nt)

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io
  
  implicit none

  ! Variables definition
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,pre1      ! velocity components and pressure
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1             ! scalar field
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1                              ! temporary array
  
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: u1mean,v1mean,w1mean             ! 1st order moment (average)
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: u2mean,v2mean,w2mean             ! 2nd order moment (variance)
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: u3mean,v3mean,w3mean             ! 3rd order moment (skewness)
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: u4mean,v4mean,w4mean             ! 4th order moment (kurtosis)
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uvmean,uwmean,vwmean             ! Reynolds stresses
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: pre1mean,pre2mean                ! average and variance of pressure
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: phi1mean,phi2mean           ! average and variance of scalar field
  real(mytype),dimension(xsize(1),xsize(2),xsize(3),nphi) :: uphimean,vphimean,wphimean  ! average and variance of mixed fluctuations
 
  integer :: is                                                                          ! index for the different scalar fields
  integer :: nr                                                                          ! number of realizations of the same flow case                                                           

  !---x-component---!
  !u1=ux1
  u1mean=u1mean+ux1/real(nr,mytype) 

  !u2=ux1*ux1
  ta1=ux1*ux1
  u2mean=u2mean+ta1/real(nr,mytype)

  !u3=ux1*ux1*ux1
  ta1=ux1*ux1*ux1
  u3mean=u3mean+ta1/real(nr,mytype)

  !u4=ux1*ux1*ux1*ux1 
  ta1=ux1*ux1*ux1*ux1
  u4mean=u4mean+ta1/real(nr,mytype)

  !---y-component---!
  !v1=uy1
  v1mean=v1mean+uy1/real(nr,mytype)

  !v2=uy1*uy1
  ta1=uy1*uy1
  v2mean=v2mean+ta1/real(nr,mytype)

  !v3=uy1*uy1*uy1 
  ta1=uy1*uy1*uy1
  v3mean=v3mean+ta1/real(nr,mytype)

  !v4=uy1*uy1*uy1*uy1 
  ta1=uy1*uy1*uy1*uy1
  v4mean=v4mean+ta1/real(nr,mytype)

  !---z-component---!
  !w1=uz1
  w1mean=w1mean+uz1/real(nr,mytype)

  !w2=uz1*uz1
  ta1=uz1*uz1
  w2mean=w2mean+ta1/real(nr,mytype)

  !w3=uz1*uz1*uz1 
  ta1=uz1*uz1*uz1
  w3mean=w3mean+ta1/real(nr,mytype)

  !w4=uz1*uz1*uz1*uz1 
  ta1=uz1*uz1*uz1*uz1
  w4mean=w4mean+ta1/real(nr,mytype)

  !---Reynolds stresses---!
  !uv=ux1*uy1
  ta1=ux1*uy1
  uvmean=uvmean+ta1/real(nr,mytype)

  !uw=ux1*uz1
  ta1=ux1*uz1
  uwmean=uwmean+ta1/real(nr,mytype)

  !vwmean=uy1*uz1
  ta1=uy1*uz1
  vwmean=vwmean+ta1/real(nr,mytype)

  !---pressure---!
  pre1mean=pre1mean+pre1/real(nr,mytype)

  ta1=pre1*pre1
  pre2mean=pre2mean+ta1/real(nr,mytype)

  !---scalar---!
  if (iscalar==1) then
     do is=1, nphi

        !phi1mean=phi1
        phi1mean(:,:,:,is)=phi1mean(:,:,:,is)+phi1(:,:,:,is)/real(nr,mytype)

        !phi2mean=phi1*phi1
        ta1=phi1(:,:,:,is)*phi1(:,:,:,is)
        phi2mean(:,:,:,is)=phi2mean(:,:,:,is)+ta1/real(nr,mytype)

        !uphimean=ux1*phi1
        ta1=ux1*phi1(:,:,:,is)
        uphimean(:,:,:,is)=uphimean(:,:,:,is)+ta1/real(nr,mytype)

        !vphimean=uy1*phi1
        ta1=uy1*phi1(:,:,:,is)
        vphimean(:,:,:,is)=vphimean(:,:,:,is)+ta1/real(nr,mytype)

        !wphimean=uz1*phi1
        ta1=uz1*phi1(:,:,:,is)
        wphimean(:,:,:,is)=wphimean(:,:,:,is)+ta1/real(nr,mytype)

     enddo
  endif

end subroutine STAT_MEAN
!********************************************************************
! Vorticity (instantaneous)
subroutine STAT_VORTICITY(ux1,uy1,uz1,ifile)

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
  
  ! saving options
  save_w = 1   ! total vorticity magnitude  
  save_w1 = 0  ! vorticity x-component
  save_w2 = 0  ! vorticity y-component
  save_w3 = 0  ! vorticity z-component
  
  nvect1=xsize(1)*xsize(2)*xsize(3)

  ! x-derivatives
  call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  ! y-derivatives
  call transpose_x_to_y(ux1,td2)
  call transpose_x_to_y(uy1,te2)
  call transpose_x_to_y(uz1,tf2)
  call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  ! z-derivatives
  call transpose_y_to_z(td2,td3)
  call transpose_y_to_z(te2,te3)
  call transpose_y_to_z(tf2,tf3)
  call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
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

  if (save_w.eq.1) then
     di1=0._mytype
     do ijk=1,nvect1
        di1(ijk,1,1)=sqrt((tf1(ijk,1,1)-th1(ijk,1,1))**2+&
             (tg1(ijk,1,1)-tc1(ijk,1,1))**2+(tb1(ijk,1,1)-td1(ijk,1,1))**2)
     enddo
     write(filename,"('./data/w',I3.3)") ifile
     call decomp_2d_write_one(1,di1,filename,1)
  endif

  if (save_w1.eq.1) then
     di1=0._mytype
     do ijk=1,nvect1 !dw/dy - dv/dz
        di1(ijk,1,1)=tf1(ijk,1,1)-th1(ijk,1,1)
     enddo
     write(filename,"('./data/w1',I3.3)") ifile
     call decomp_2d_write_one(1,di1,filename,1)
  endif

  if (save_w2.eq.1) then
     di1=0._mytype
     do ijk=1,nvect1 !du/dz - dw/dx
        di1(ijk,1,1)=tg1(ijk,1,1)-tc1(ijk,1,1)
     enddo
     write(filename,"('./data/w2',I3.3)") ifile
     call decomp_2d_write_one(1,di1,filename,1)
  endif

  if (save_w3.eq.1) then
     di1=0._mytype
     do ijk=1,nvect1 !dv/dx - du/dy
        di1(ijk,1,1)=tb1(ijk,1,1)-td1(ijk,1,1)
     enddo
     write(filename,"('./data/w3',I3.3)") ifile
     call decomp_2d_write_one(1,di1,filename,1)
  endif

end subroutine STAT_VORTICITY
!********************************************************************


