
!This file is not part of standard Xcompact3d releases (xcompact3d.com).

!-----------------------------------------------------------------------------!
! DESCRIPTION: This file contains subroutines for the calculation of the 
!              following statistics:
!              1) One-point statistics:
!                 - average, variance, skewness and kurtosis of velocity 
!                   components;
!                 - Reynolds stress;
!                 - average and variance of pressure, scalar field and mixed 
!                   fluctuations (velocity-scalar);
!                 - mean vorticity (x, y, z) and mean gradients 
!                   (x, z and scalar field);
!                 - total dissipation rate;
!                 - terms of TKE equation that cannot be obtained by mean 
!                   flow statistics. 
!              2) Two-points statistics:
!                 - spanwise correlation functions for velocity components 
!                   (Ruuz, Rvvz, Rwwz, Ruvz);
!                 - spanwise correlation function for scalar field (s: scalar)
!                   (Rssz).
!              Adapted from original Incompact3d file (v2.0) 
!              of Roberto Corsini. 
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
!              Roberto Corsini <roberto.corsini@unimore.it> 
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculate mean statistics (average, variance, skewness, 
!              kurtosis, Reynolds stresses and mixed fluctuations).
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
!              Roberto Corsini <roberto.corsini@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine stat_mean(ux2,uy2,uz2,pre2,phi2,                           &
                     u1mean,v1mean,w1mean,u2mean,v2mean,w2mean,       &
                     u3mean,v3mean,w3mean,u4mean,v4mean,w4mean,       &
                     uvmean,uwmean,vwmean,pre1mean,pre2mean,vpremean, &
                     phi1mean,phi2mean,uphimean,vphimean,wphimean)

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use param
  use variables
  use post_processing, only : nr,nt,den

  implicit none
  
  ! Variables definition (velocity components, pressure and scalar field)
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3))             :: ux2,uy2,uz2,pre2          
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3),1:numscalar) :: phi2 
    
  ! Local 
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2                                       ! temporary array (local)
  
  ! Output
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: u1mean,v1mean,w1mean        ! 1st order moment (average)
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: u2mean,v2mean,w2mean        ! 2nd order moment (variance)
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: u3mean,v3mean,w3mean        ! 3rd order moment (skewness)
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: u4mean,v4mean,w4mean        ! 4th order moment (kurtosis)
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: uvmean,uwmean,vwmean        ! Reynolds stresses
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: pre1mean,pre2mean           ! average and variance of pressure
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: vpremean                    ! average of pressure-strain term in y-direction
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: phi1mean,phi2mean           ! average and variance of scalar field
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: uphimean,vphimean,wphimean  ! average of mixed fluctuations for scalar field
                                                                                                                                               
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
  ! <uv>
  ta2=ux2*uy2
  uvmean=uvmean+ta2/den

  ! <uw>
  ta2=ux2*uz2
  uwmean=uwmean+ta2/den

  ! <vw>
  ta2=uy2*uz2
  vwmean=vwmean+ta2/den

  !---pressure---!
  ! average
  pre1mean=pre1mean+pre2/den
  
  ! variance 
  ta2=pre2*pre2
  pre2mean=pre2mean+ta2/den
  
  ! pressure-strain in y-direction
  ta2=uy2*pre2
  vpremean=vpremean+ta2/den 

  !---scalar---!
  if (iscalar==1) then

        ! average: phi
        phi1mean=phi1mean+phi2(:,:,:,1)/den

        ! variance: phi
        ta2=phi2(:,:,:,1)*phi2(:,:,:,1)
        phi2mean=phi2mean+ta2/den

        ! mixed fluctuations: u',phi'
        ta2=ux2*phi2(:,:,:,1)
        uphimean=uphimean+ta2/den

        ! mixed fluctuations: v',phi'
        ta2=uy2*phi2(:,:,:,1)
        vphimean=vphimean+ta2/den

        ! mixed fluctuations: w',phi'
        ta2=uz2*phi2(:,:,:,1)
        wphimean=wphimean+ta2/den

  endif

end subroutine stat_mean

!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculate mean vorticity components, mean velocity gradients,  
!              mean scalar gradient and total dissipation rate.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
!              Roberto Corsini <roberto.corsini@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine stat_gradients(ux1,uy1,uz1,phi1,                              &
                          vortxmean2,vortymean2,vortzmean2,              &
                          mean_gradientx2,mean_gradientz2,mean_gradphi2, &
                          epsmean2)   

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use param
  use variables
  
  use var,  only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  use var,  only : ta2,tb2,tc2,td2,te2,tf2,di2
  use var,  only : ta3,tb3,tc3,td3,te3,tf3,di3
  
  use post_processing, only : nr,nt,den
  
  implicit none
  
  ! Variables definition (velocity components and scalar field)
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3))             :: ux1,uy1,uz1
  real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),1:numscalar) :: phi1
                                                                                                              
  ! IBM parameter
  real(mytype) :: lind
    
  ! Vorticity (average vorticity components, y-pencils)
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: vortxmean2,vortymean2,vortzmean2   
  
  ! Mean gradients (mean gradients, y-pencils) (x:streamwise, z: spanwise)
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: mean_gradientx2, mean_gradientz2
  
  ! Mean scalar gradient (y-pencils)
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: mean_gradphi2
  
  ! Total average dissipation (y-pencils)
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: epsmean2  
        
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
  call derz (tc3,tf3,di3,sz,ffz, fsz, fwz, zsize(1),zsize(2),zsize(3),0,lind)
  
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
    
  !---- Mean gradients ----!
    
  ! Mean streamwise gradient
  di1 = td1 ! du/dy
  
  ! Transpose array along y and sum
  call transpose_x_to_y(di1,di2)
  mean_gradientx2 = mean_gradientx2 + di2/den
  
  ! Mean spanwise gradient
  di1 = tf1 ! dw/dy
  
  ! Transpose array along y and sum
  call transpose_x_to_y(di1,di2)
  mean_gradientz2 = mean_gradientz2 + di2/den
  
   
  !---Vorticity average---!
  
  ! Vorticity along x 
  di1 = tf1 - th1  !dw/dy - dv/dz
  
  ! Transpose array along y and sum
  call transpose_x_to_y(di1,di2)
  vortxmean2 = vortxmean2 + di2/den
    
  ! Vorticity along y
  di1 = tg1 - tc1  !du/dz - dw/dx
  
  ! Transpose array along y and sum
  call transpose_x_to_y(di1,di2)
  vortymean2 = vortymean2 + di2/den
     
  ! Vorticity along z
  di1 = tb1 - td1  !dv/dx - du/dy
  
  ! Transpose array along y and sum
  call transpose_x_to_y(di1,di2)
  vortzmean2 = vortzmean2 + di2/den
  
  !--- Scalar gradient ---!
  if (numscalar == 1) then
  
      ! Transpose to y-pencils
      call transpose_x_to_y(phi1(:,:,:,1),td2)
  
      ! y-derivative
      call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,lind)
  
      ! Summation
      mean_gradphi2 = mean_gradphi2 + ta2/den
   
  end if

  !---- Instantaneous total dissipation rate ----!
    
  di1 = (2*ta1**2 + 2*te1**2 + 2*ti1**2 + &  !(2*du/dx**2 + 2*dv/dy**2 + 2*dw/dz**2 +
         (tb1+td1)**2 +                   &  !(dv/dx+du/dy)**2 +
         (tf1+th1)**2 +                   &  !(dw/dy+dv/dz)**2 +
         (tg1+tc1)**2)*xnu                   !(du/dz+dw/dx)**2)*nu
  
  ! Transpose array along y and sum
  call transpose_x_to_y(di1,di2)
  epsmean2 = epsmean2 + di2/den
  
end subroutine stat_gradients

!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculate the (auto) correlation functions Rii and Ruv in 
!              z-direction for the velocity components.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
!              Roberto Corsini <roberto.corsini@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine stat_correlation_z(ux2,uy2,uz2,phi2,RuuzH1,RvvzH1,RwwzH1,RuvzH1,RppzH1)

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use variables,       only : numscalar,nz
  use post_processing, only : read_phi,nr,nt
  use var,             only : phi3
  
  use post_processing, only : den

  implicit none
   
  ! Variables definition (velocity and scalar field fluctuations, y-pencils)
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3))           :: ux2,uy2,uz2          
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3),numscalar) :: phi2 
    
  ! Local work arrays
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3,ta3
  
  ! Correlation functions (first index: j (rows); second index: r (columns))
  real(mytype),intent(inout),dimension(zsize(2),zsize(3)) :: RuuzH1, RvvzH1, RwwzH1, RuvzH1, RppzH1
    
  ! Indexes for cycles          
  integer :: i,j,k,rr,kpr 

  ! Transpose arrays along z
  call transpose_y_to_z(ux2,ux3)
  call transpose_y_to_z(uy2,uy3)
  call transpose_y_to_z(uz2,uz3)
  
  if(read_phi) call transpose_y_to_z(phi2(:,:,:,1),phi3(:,:,:,1))

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
                  
                  !--- Vertical fluctuations correlation ---!
                  
                  ! Product of fluctuations at distance 'r'
                  ta3(i,j,k) = uy3(i,j,k)*uy3(i,j,kpr)
                  
                  ! Accumulation inside the correlation function variable (at each subdomain)
                  RvvzH1(j,rr) = RvvzH1(j,rr) + ta3(i,j,k)/den
                  
                  !--- Spanwise fluctuations correlation ---!
                  
                  ! Product of fluctuations at distance 'r'
                  ta3(i,j,k) = uz3(i,j,k)*uz3(i,j,kpr)
                  
                  ! Accumulation inside the correlation function variable (at each subdomain)
                  RwwzH1(j,rr) = RwwzH1(j,rr) + ta3(i,j,k)/den
                  
                  !--- Mixed fluctuations correlation (u'v') ---!
                  
                  ! Product of fluctuations at distance 'r'
                  ta3(i,j,k) = ux3(i,j,k)*uy3(i,j,kpr)
                  
                  ! Accumulation inside the correlation function variable (at each subdomain)
                  RuvzH1(j,rr) = RuvzH1(j,rr) + ta3(i,j,k)/den
                  
                  !--- Scalar fluctuations correlation (phi'phi') ---!
                  if (numscalar == 1) then
                  
                      ! Product of fluctuations at distance 'r'
                      ta3(i,j,k) = phi3(i,j,k,1)*phi3(i,j,kpr,1)
                  
                      ! Accumulation inside the correlation function variable (at each subdomain)
                      RppzH1(j,rr) = RppzH1(j,rr) + ta3(i,j,k)/den
                  
                  end if
                  
              enddo
          enddo
      enddo
  enddo

end subroutine stat_correlation_z

!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculation of Turbulent Kinetic Energy (TKE) terms that cannot 
!              be obtained from mean flow statistics data.
!              Used for for Channels and TTBLs (statistical homogeneity 
!              in x and z directions). All results obtained here are 
!              in y-pencils. Results are later derived in y-direction.
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
!-----------------------------------------------------------------------------!
subroutine extra_terms_tke(ux2,uy2,uz2,kvprime_mean,pseudo_eps_tke_mean)   

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use param
  use variables
  
  use var,  only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  use var,  only : ta2,tb2,tc2,td2,te2,tf2,di2
  use var,  only : ta3,tb3,tc3,td3,te3,tf3,di3
  
  use post_processing, only : nr,nt,den
  
  implicit none
  
  ! Fluctuations
  real(mytype),intent(in),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2
                                                                                                   
  ! Arrays for derivative calculations
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    
  ! IBM parameter
  real(mytype) :: lind
  
  ! Turbulent transport of TKE by v'
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: kvprime_mean
    
  ! Pseudo-dissipation for TKE
  real(mytype),intent(inout),dimension(ysize(1),ysize(2),ysize(3)) :: pseudo_eps_tke_mean   
            
  ! Turbulent transport of TKE by v' [ 0.5 * (u'^2 + v'^2 + w'^2) * v']
  kvprime_mean = kvprime_mean + zpfive*(ux2**2 + uy2**2 + uz2**2)*uy2 / den
  
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
  call derz (tc3,tf3,di3,sz,ffz, fsz, fwz, zsize(1),zsize(2),zsize(3),0,lind)
  
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
  pseudo_eps_tke_mean = pseudo_eps_tke_mean + di2/den
    
end subroutine extra_terms_tke
!-----------------------------------------------------------------------------!


