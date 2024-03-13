!##################################################################
program mesh_evaluation

  implicit none

  ! Inputs
  real(8)            :: yly = 2.0    ! domain dimension in y direction
  integer, parameter :: ny = 65      ! number of points in y direction
  integer, parameter :: istret = 3   ! y mesh refinement (0:no, 1:center, 2:both sides, 3:bottom)
  real(8)            :: beta = 0.25  ! beta parameter for mesh stretching
  
  
  ! Declare local constant Pi
  real, parameter    :: pi = 3.1415927
  
  ! Work variables
  real(8) :: yinf, den, xnum, alpha, xnum1, cst
  integer :: j, nym
  real(8) :: den1, den3, den4 
  
  ! If periodic BC is imposed, comment the following line
  nym=ny-1
  
  ! Variables 
  real(8), dimension(ny)  ::  yeta, yp
  real(8), dimension(nym) ::  ypi
  
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
           if (yeta(j).lt. 0.5) yp(j)=xnum1-cst-yinf
           if (yeta(j).eq. 0.5) yp(j)=0.0-yinf
           if (yeta(j).gt. 0.5) yp(j)=xnum1+cst-yinf
        endif
        if (istret==2) then
           if (yeta(j).lt. 0.5) yp(j)=xnum1-cst+yly
           if (yeta(j).eq. 0.5) yp(j)=0.0+yly
           if (yeta(j).gt. 0.5) yp(j)=xnum1+cst+yly
        endif
        if (istret==3) then
           if (yeta(j).lt. 0.5) yp(j)=(xnum1-cst+yly)*2.0
           if (yeta(j).eq. 0.5) yp(j)=(0.0+yly)*2.0
           if (yeta(j).gt. 0.5) yp(j)=(xnum1+cst+yly)*2.0
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
        if (istret==1) yetai(j)=(j - 0.5)*(1.0/nym)
        if (istret==2) yetai(j)=(j - 0.5)*(1.0/nym) - 0.5
        if (istret==3) yetai(j)=(j - 0.5)*(0.5/nym) - 0.5
        
        den1=sqrt(alpha*beta + 1.0)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
        den=2.0*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
        den4=2.0*alpha*beta-cos(2.0*pi*yetai(j)) + 1.0
        xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(2.0*sqrt(alpha)*sqrt(alpha*beta + 1.0))
        
        if (istret==1) then
           if (yetai(j).lt. 0.5) ypi(j)=xnum1-cst-yinf
           if (yetai(j).eq. 0.5) ypi(j)=0.0-yinf
           if (yetai(j).gt. 0.5) ypi(j)=xnum1+cst-yinf
        endif
        if (istret==2) then
           if (yetai(j).lt. 0.5) ypi(j)=xnum1-cst+yly
           if (yetai(j).eq. 0.5) ypi(j)=0.0+yly
           if (yetai(j).gt. 0.5) ypi(j)=xnum1+cst+yly
        endif
        if (istret==3) then
           if (yetai(j).lt. 0.5) ypi(j)=(xnum1-cst+yly)*2.0
           if (yetai(j).eq. 0.5) ypi(j)=(0.0+yly)*2.0
           if (yetai(j).gt. 0.5) ypi(j)=(xnum1+cst+yly)*2.0
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

  ! Writing the y coordinates of the trial mesh' nodes and centers
  if (nrank == 0) then
     open(10,file='yp_trial.dat', form='formatted')
     do j=1,ny
        write(10,*)yp(j)
     enddo
     close(10)
     open(10,file='ypi_trial.dat', form='formatted')
     do j=1,nym
        write(10,*)ypi(j)
     enddo
     close(10)
  endif

end subroutine mesh_evaluation
