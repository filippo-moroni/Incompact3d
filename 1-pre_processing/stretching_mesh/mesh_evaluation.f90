!----------------------------------------------------------!
!     This program is used to estimate the dimension       !
!    of the grid elements for a temporal BL simulation.    !
!      Adapted from original Incompact3d file (master)     !
!            stretching_parameter_channel.f90              !
!----------------------------------------------------------!
program mesh_evaluation

  implicit none

  ! Inputs
  real(8)            :: yly = 2.0    ! domain dimension in y direction
  integer, parameter :: ny = 65      ! number of points in y direction
  integer, parameter :: nym = 64     ! if periodic BC is imposed, nym = ny, otherwise nym = ny - 1
  integer, parameter :: istret = 3   ! y mesh refinement (0:no, 1:center, 2:both sides, 3:bottom)
  real(8)            :: beta = 0.25  ! beta parameter for mesh stretching
  real(8)            :: cf = 0.007   ! maximum cf estimated 
  real(8)            :: nu = 0.002   ! kinematic viscosity
  real(8)            :: uwall = 1.0  ! velocity of the wall
  
  ! Declare local constant Pi
  real, parameter    :: pi = 3.1415927
  
  ! Work variables
  real(8) :: yinf, den, xnum, alpha, xcx, xnum1, cst
  real(8) :: den1, den3, den4
  real(8) :: sh_vel, delta_nu
  real(8) :: delta_y1, delta_yn, height
  integer :: j, npvis
    
  ! Variables 
  real(8), dimension(ny)  ::  yeta, yp, yetai, ypi
    
  ! Start of calculations as the original stretching subroutine in Incompact3d
  yinf=-yly/2.0
  den=2.0*beta*yinf
  xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
  alpha=abs(xnum/den)
  xcx=1.0/beta/alpha
  
  if (alpha.ne.0.0) then
     if (istret.eq.1) yp(1)=0.0
     if (istret.eq.2) yp(1)=0.0
     if (istret.eq.1) yeta(1)=0.0
     if (istret.eq.2) yeta(1)=-0.5
     if (istret.eq.3) yp(1)=0.0
     if (istret.eq.3) yeta(1)=-0.5
     
     do j=2,ny
        if (istret==1) yeta(j)=(j - 1.0)*(1.0/nym)
        if (istret==2) yeta(j)=(j - 1.0)*(1.0/nym) - 0.5
        if (istret==3) yeta(j)=(j - 1.0)*(0.5/nym) - 0.5
        
        den1=sqrt(alpha*beta + 1.0)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
        den=2.0*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
        den4=2.0*alpha*beta-cos(2.0*pi*yeta(j)) + 1.0
        xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(2.0*sqrt(alpha)*sqrt(alpha*beta + 1.0))
        
        if (istret==1) then
        
           if (yeta(j).lt. 0.5) then          
           	yp(j)=xnum1-cst-yinf 
           endif
           
           if (yeta(j).eq. 0.5) then 
           	yp(j)=0.0-yinf
           endif
           
           if (yeta(j).gt. 0.5) then
           	yp(j)=xnum1+cst-yinf
           endif
        endif
        
        if (istret==2) then
           if (yeta(j).lt. 0.5) then 
           	yp(j)=xnum1-cst+yly
           endif
           
           if (yeta(j).eq. 0.5) then 
           	yp(j)=0.0+yly
           endif
           
           if (yeta(j).gt. 0.5) then
           	yp(j)=xnum1+cst+yly
           endif
        endif
        
        if (istret==3) then
           if (yeta(j).lt. 0.5) then
           	yp(j)=(xnum1-cst+yly)*2.0
           endif
           
           if (yeta(j).eq. 0.5) then
           	yp(j)=(0.0+yly)*2.0
           endif
           if (yeta(j).gt. 0.5) then
           	yp(j)=(xnum1+cst+yly)*2.0
           endif
        endif
        
     enddo
  endif
  
  if (alpha.eq.0.) then
     yp(1)=-1.e10
     do j=2,ny
        yeta(j)=(j - 1.0)*(1.0/ny)
        yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
     enddo
  endif
  
  if (alpha.ne.0.) then
     do j=1,ny
        if (istret==1) then
        	yetai(j)=(j - 0.5)*(1.0/nym)
        endif
        
        if (istret==2) then
        	yetai(j)=(j - 0.5)*(1.0/nym) - 0.5
        endif
        
        if (istret==3) then
        	yetai(j)=(j - 0.5)*(0.5/nym) - 0.5
        endif
        
        den1=sqrt(alpha*beta + 1.0)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
        den=2.0*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
        den4=2.0*alpha*beta-cos(2.0*pi*yetai(j)) + 1.0
        xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(2.0*sqrt(alpha)*sqrt(alpha*beta + 1.0))
        
        if (istret==1) then
           if (yetai(j).lt. 0.5) then
           	ypi(j)=xnum1-cst-yinf
           endif
           
           if (yetai(j).eq. 0.5) then
           	ypi(j)=0.0-yinf
           endif
           
           if (yetai(j).gt. 0.5) then
           	ypi(j)=xnum1+cst-yinf
           endif          
        endif
        
        if (istret==2) then
           if (yetai(j).lt. 0.5) then
           	ypi(j)=xnum1-cst+yly
           endif
           
           if (yetai(j).eq. 0.5) then
           	ypi(j)=0.0+yly
           endif
           
           if (yetai(j).gt. 0.5) then
           	ypi(j)=xnum1+cst+yly
           endif
        endif
        
        if (istret==3) then
           if (yetai(j).lt. 0.5) then
           	ypi(j)=(xnum1-cst+yly)*2.0 
           endif
           
           if (yetai(j).eq. 0.5) then
           	ypi(j)=(0.0+yly)*2.0
           endif
           
           if (yetai(j).gt. 0.5) then 
           	ypi(j)=(xnum1+cst+yly)*2.0
           endif
        endif
        
     enddo
  endif
  
  if (alpha.eq.0.0) then
     ypi(1)=-1.e10
     do j=2,ny
        yetai(j)=(j - 1.0)*(1.0/ny)
        ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
     enddo
  endif
  
  ! This part is valid for meshes with refinement at the bottom boundary only
  if (istret .eq. 3) then
   
  ! Calculate shear velocity
  sh_vel = sqrt((cf/2.0))*uwall

  ! Calculate viscous unit
  delta_nu = nu/sh_vel
  
  ! Rescaling the coordinates
  yp = yp/delta_nu
  ypi = ypi/delta_nu
 
  ! First and last elements' dimension
  delta_y1 = yp(2) - yp(1)
  delta_yn = yp(ny) - yp(ny-1)
 
  ! Calculation of the number of mesh nodes in the viscous sublayer
  npvis=0     ! number of points viscous sublayer
  height=0.0  ! cumulative height in viscous unit (y+)
  
  do j=2,ny    
     if (height + yp(j) - yp(j-1) .le. 5) npvis = npvis + 1 
     height = height + yp(j) - yp(j-1)
 enddo
  
  ! Printing the useful informations to the screen
  print *,'Number of mesh nodes in wall normal direction : ',ny
  print *,'Skin friction coefficient employed = ',cf
  print *,'Beta parameter = ',beta
  print *
  print *,'Mesh size at the first element near the wall: delta_y1+ = ',delta_y1
  print *,'Mesh size at the last element away from the wall: delta_yn+ = ',delta_yn
  print *,'Number of mesh nodes in viscous sublayer : ', npvis
  print *
  
  endif
  
  ! Writing the y coordinates of the trial mesh' nodes and centers
  !open(10,file='yp_trial.dat', form='formatted')
  !   do j=1,ny
  !      write(10,*)yp(j)
  !   enddo
  !close(10)
  !open(10,file='ypi_trial.dat', form='formatted')
  !   do j=1,nym
  !      write(10,*)ypi(j)
  !   enddo
  !close(10)

end program mesh_evaluation


