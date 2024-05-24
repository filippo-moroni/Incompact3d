!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################
module tools

  implicit none

  logical, save :: adios2_restart_initialised = .false.

  character(len=*), parameter :: io_restart = "restart-io"
  character(len=*), parameter :: resfile = "checkpoint"
  character(len=*), parameter :: io_ioflow = "in-outflow-io"
  
  private

  public :: reading_input_file,   &
            test_scalar_min_max,  & 
            test_speed_min_max,   &
            simu_stats,           &
            apply_spatial_filter, &
            mean_plane_x,         &
            mean_plane_y,         &
            mean_plane_z,         &
            avg3d

contains

  !##################################################################
  subroutine reading_input_file()

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE var
  USE MPI
  
  implicit none
  
  ! Variables to read the input.i3d file
  integer :: nargin, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase
  
  ! Reading of the input file as Xcompact3d does
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     if (nrank==0) write(*,*) 'PostIncompact3d is run with the default file -->', trim(InputFN)
  elseif (nargin >= 1) then
     call get_command_argument(1,InputFN,FNLength,status)
     back=.true.
     FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
     DecInd=index(FNBase,'.',back)
     if (DecInd >1) then
        FNBase=FNBase(1:(DecInd-1))
     end if
     if (nrank==0) write(*,*) 'PostIncompact3d is run with the provided file -->', trim(InputFN)
  endif
  
  ! Reading the input file for geometry and numerics
  call parameter(InputFN)
  
  end subroutine reading_input_file
  !##################################################################
  subroutine test_scalar_min_max(phi)

    use decomp_2d
    use variables
    use param
    use var
    use mpi
    use dbg_schemes, only: abs_prec

    implicit none

    integer :: code,ierror,i,j,k,is,jglob
    real(mytype) :: phimax,phimin,phimax1,phimin1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(2,numscalar) :: phimaxin,phimaxout

    do is=1, numscalar

      ta1(:,:,:) = phi(:,:,:,is)
      ! ibm
      if (iibm > 0) then
        ta1(:,:,:) = (one - ep1(:,:,:)) * ta1(:,:,:)
      endif

      phimax=-1609._mytype
      phimin=1609._mytype
      phimax = maxval(ta1(:,:,:))
      phimin =-minval(ta1(:,:,:))
      phimaxin(:,is) =  (/phimin, phimax /)
    enddo

    call MPI_REDUCE(phimaxin,phimaxout,numscalar*2,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)

    do is=1,numscalar
      if (nrank == 0) then
        phimin1 = -phimaxout(1,is)
        phimax1 =  phimaxout(2,is)

        write(*,*) 'Phi'//char(48+is)//' min max=', real(phimin1,4), real(phimax1,4)

        if (abs_prec(phimax1) > 100._mytype) then !if phi control turned off
           write(*,*) 'Scalar diverged! SIMULATION IS STOPPED!'
           call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
        endif
      endif
    enddo

    return
  end subroutine test_scalar_min_max
  !##################################################################
  subroutine test_speed_min_max(ux,uy,uz)

    use decomp_2d
    use variables
    use param
    use var
    use mpi
    use dbg_schemes, only: abs_prec

    implicit none

    integer :: code,ierror,i,j,k
    real(mytype) :: uxmax,uymax,uzmax,uxmin,uymin,uzmin
    real(mytype) :: uxmax1,uymax1,uzmax1,uxmin1,uymin1,uzmin1
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux,uy,uz
    real(mytype),dimension(6) :: umaxin, umaxout

    if (iibm > 0) then
       ux(:,:,:) = (one - ep1(:,:,:)) * ux(:,:,:)
       uy(:,:,:) = (one - ep1(:,:,:)) * uy(:,:,:)
       uz(:,:,:) = (one - ep1(:,:,:)) * uz(:,:,:)
    endif

    uxmax=-1609.;uymax=-1609.;uzmax=-1609.;uxmin=1609.;uymin=1609.;uzmin=1609.
    !
    ! More efficient version
    uxmax=maxval(ux)
    uymax=maxval(uy)
    uzmax=maxval(uz)
    uxmin=-minval(ux)
    uymin=-minval(uy)
    uzmin=-minval(uz)

    umaxin = (/uxmax, uymax, uzmax, uxmin, uymin, uzmin/)
    call MPI_REDUCE(umaxin,umaxout,6,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)

    uxmax1= umaxout(1)
    uymax1= umaxout(2)
    uzmax1= umaxout(3)
    uxmin1=-umaxout(4)
    uymin1=-umaxout(5)
    uzmin1=-umaxout(6)

    if (nrank == 0) then

       write(*,*) 'U,V,W min=',real(uxmin1,4),real(uymin1,4),real(uzmin1,4)
       write(*,*) 'U,V,W max=',real(uxmax1,4),real(uymax1,4),real(uzmax1,4)
       !print *,'CFL=',real(abs(max(uxmax1,uymax1,uzmax1)*dt)/min(dx,dy,dz),4)

       if((abs_prec(uxmax1)>=onehundred).or.(abs_prec(uymax1)>=onehundred).OR.(abs_prec(uzmax1)>=onehundred)) then
         write(*,*) 'Velocity diverged! SIMULATION IS STOPPED!'
         call MPI_ABORT(MPI_COMM_WORLD,code,ierror)
         stop
       endif

    endif

    return
  end subroutine test_speed_min_max
  !##################################################################
  subroutine simu_stats(iwhen)

    use decomp_2d
    use simulation_stats
    use var
    use MPI

    implicit none

    integer :: iwhen

    if (iwhen == 1) then !AT THE START OF THE SIMULATION
       tstart=zero
       time1=zero
       trank=zero
       tranksum=zero
       ttotal=zero
       call cpu_time(tstart)
    else if (iwhen == 2) then !AT THE START OF A TIME STEP
       if (nrank == 0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime==ilast)) then
          call cpu_time(time1)
          write(*,*) '==========================================================='
          write(*,"(' Time step =',i7,'/',i7,', Time unit =',F9.4)") itime,ilast,t
       endif
    else if ((iwhen == 3).and.(itime > ifirst)) then !AT THE END OF A TIME STEP
       if (nrank == 0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime==ilast)) then
          call cpu_time(trank)
          if (nrank==0) write(*,*) 'Time for this time step (s):',real(trank-time1)
          telapsed = (trank-tstart)/threethousandsixhundred
          tremaining  = telapsed*(ilast-itime)/(itime-ifirst)
          write(*,"(' Remaining time:',I8,' h ',I2,' min')") int(tremaining), int((tremaining-int(tremaining))*sixty)
          write(*,"(' Elapsed time:  ',I8,' h ',I2,' min')") int(telapsed), int((telapsed-int(telapsed))*sixty)
       endif
    else if (iwhen == 4) then !AT THE END OF THE SIMULATION
       call cpu_time(trank)
       ttotal=trank-tstart
       if (nrank == 0) then
          write(*,*) '==========================================================='
          write(*,*) '                                                           '
          write(*,*) 'Good job! Xcompact3d finished successfully!                '
          write(*,*) '                                                           '
          write(*,*) '2DECOMP with p_row*p_col=',p_row,p_col
          write(*,*) '                                                           '
          write(*,*) 'nx*ny*nz=',nx*ny*nz
          write(*,*) 'nx,ny,nz=',nx,ny,nz
          write(*,*) 'dx,dy,dz=',dx,dy,dz
          write(*,*) '                                                           '
          write(*,*) 'Averaged time per step (s):',real(ttotal/(ilast-(ifirst-1)),4)
          write(*,*) 'Total wallclock (s):',real(ttotal,4)
          write(*,*) 'Total wallclock (m):',real(ttotal/sixty,4)
          write(*,*) 'Total wallclock (h):',real(ttotal/threethousandsixhundred,4)
          write(*,*) '                                                           '
       endif
    endif

  end subroutine simu_stats
  !##############################################################################
  subroutine init_restart_adios2()

    use decomp_2d, only : mytype, phG
    use decomp_2d_io, only : decomp_2d_register_variable, decomp_2d_init_io
    use variables, only : numscalar
    use param, only : ilmn, nrhotime, ntime
    use var, only : itimescheme, iibm
    
    implicit none

    integer :: ierror
    
    integer :: is
    character(len=80) :: varname
    
    call decomp_2d_init_io(io_restart)
    
    call decomp_2d_register_variable(io_restart, "ux", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart, "uy", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart, "uz", 1, 0, 0, mytype)

    call decomp_2d_register_variable(io_restart, "pp", 3, 0, 0, mytype, phG) !! XXX: need some way to handle the different grid here...

    do is = 1, numscalar
       write(varname,*) "phi-", is
       call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
    end do

    if ((itimescheme.eq.2) .or. (itimescheme.eq.3)) then
       call decomp_2d_register_variable(io_restart, "dux-2", 1, 0, 0, mytype)
       call decomp_2d_register_variable(io_restart, "duy-2", 1, 0, 0, mytype)
       call decomp_2d_register_variable(io_restart, "duz-2", 1, 0, 0, mytype)

       do is = 1, numscalar
          write(varname,*) "dphi-", is, "-2"
          call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
       end do

       if (itimescheme.eq.3) then
          call decomp_2d_register_variable(io_restart, "dux-3", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "duy-3", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "duz-3", 1, 0, 0, mytype)

          do is = 1, numscalar
             write(varname,*) "dphi-", is, "-3"
             call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
          end do
       endif
    endif

    if (iibm .ne. 0) then
       call decomp_2d_register_variable(io_restart, "ep", 1, 0, 0, mytype)
    endif

    if (ilmn) then
       do is = 1, nrhotime
          write(varname, *) "rho-", is
          call decomp_2d_register_variable(io_restart, varname, 1, 0, 0, mytype)
       end do
       do is = 1, ntime
          write(varname, *) "drho-", is
          call decomp_2d_register_variable(io_restart, varname, 1, 0, 0, mytype)
       end do
    end if
    
  end subroutine init_restart_adios2
  !############################################################################
  !!  SUBROUTINE: apply_spatial_filter
  !############################################################################
  subroutine apply_spatial_filter(ux1,uy1,uz1,phi1)

    use decomp_2d
    use param
    use var, only: uxf1,uyf1,uzf1,uxf2,uyf2,uzf2,uxf3,uyf3,uzf3,di1,di2,di3,phif1,phif2,phif3
    use variables
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(inout) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3), numscalar), intent(inout) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi11
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2, phi2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3, phi3

    integer :: i,j,k,npaire

    !if (iscalar == 1) phi11=phi1(:,:,:,1) !currently only first scalar
    if (ifilter==1.or.ifilter==2) then
      call filx(uxf1,ux1,di1,fisx,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0,ubcx)
      call filx(uyf1,uy1,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
      call filx(uzf1,uz1,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    else
      uxf1=ux1
      uyf1=uy1
      uzf1=uz1
      !if (iscalar == 1) phif1=phi11
    end if

    call transpose_x_to_y(uxf1,ux2)
    call transpose_x_to_y(uyf1,uy2)
    call transpose_x_to_y(uzf1,uz2)
    !if (iscalar == 1) call transpose_x_to_y(phif1,phi2)

    if (ifilter==1.or.ifilter==3) then ! all filter or y filter
      call fily(uxf2,ux2,di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcx)
      call fily(uyf2,uy2,di2,fisy,fiffy,fifsy,fifwy,ysize(1),ysize(2),ysize(3),0,ubcy)
      call fily(uzf2,uz2,di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcz)
      !if (iscalar.eq.1) call fily(phif2,phi2,di2,fisy,fiffy,fifsy,fifwy,ysize(1),ysize(2),ysize(3),0)
    else
      uxf2=ux2
      uyf2=uy2
      uzf2=uz2
      !if (iscalar == 1) phif2=phi2
    end if

    call transpose_y_to_z(uxf2,ux3)
    call transpose_y_to_z(uyf2,uy3)
    call transpose_y_to_z(uzf2,uz3)
    !if (iscalar == 1) call transpose_y_to_z(phif2,phi3)

    if (ifilter==1.or.ifilter==2) then
      call filz(uxf3,ux3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
      call filz(uyf3,uy3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
      call filz(uzf3,uz3,di3,fisz,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),0,ubcz)
      !if (iscalar.eq.1) call filz(phif3,phi3,di3,fisz,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),0)
    else
      uxf3=ux3
      uyf3=uy3
      uzf3=uz3
      !if (iscalar == 1) phif3=phi3
    end if

    call transpose_z_to_y(uxf3,ux2)
    call transpose_z_to_y(uyf3,uy2)
    call transpose_z_to_y(uzf3,uz2)
    !if (iscalar == 1) call transpose_z_to_y(phif3,phi2)

    call transpose_y_to_x(ux2,ux1)
    call transpose_y_to_x(uy2,uy1)
    call transpose_y_to_x(uz2,uz1)
    !if (iscalar == 1) call transpose_y_to_x(phi2,phi11)

    !if (iscalar == 1) phi1(:,:,:,1)=phi11

  end subroutine apply_spatial_filter
  !############################################################################
  subroutine mean_plane_x (f1,nx,ny,nz,fm1)

    use param, only : mytype, zero

    implicit none

    integer,intent(in) :: nx, ny, nz
    real(mytype),intent(in),dimension(nx,ny,nz) :: f1
    real(mytype),intent(out),dimension(ny,nz) :: fm1
    integer :: i,j,k

    fm1 = sum(f1, DIM=1) / real(nx, mytype)
    return

  end subroutine mean_plane_x
  !##################################################################
  subroutine mean_plane_y (f2,nx,ny,nz,fm2)

    use param, only : mytype, zero

    implicit none

    integer,intent(in) :: nx, ny, nz
    real(mytype),intent(in),dimension(nx,ny,nz) :: f2
    real(mytype),intent(out),dimension(nx,nz) :: fm2
    integer :: i,j,k

    fm2 = sum(f2, DIM=2) / real(ny, mytype)
    return

  end subroutine mean_plane_y
  !##################################################################
  subroutine mean_plane_z (f3,nx,ny,nz,fm3)

    use param, only : mytype, zero

    implicit none

    integer,intent(in) :: nx, ny, nz
    real(mytype),intent(in),dimension(nx,ny,nz) :: f3
    real(mytype),intent(out),dimension(nx,ny) :: fm3
    integer :: i,j,k

    fm3 = sum(f3, DIM=3) / real(nz,mytype)
    return

  end subroutine mean_plane_z
  !############################################################################
  !!
  !!  SUBROUTINE: avg3d
  !!      AUTHOR: Stefano Rolfo
  !! DESCRIPTION: Compute the total sum of a a 3d field
  !!
  !############################################################################
  subroutine avg3d (var, avg)

    use decomp_2d, only: real_type, xsize, xend
    use param
    use variables, only: nx,ny,nz,nxm,nym,nzm
    use mpi

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: var
    real(mytype), intent(out) :: avg
    real(mytype)              :: dep

    integer :: i,j,k, code
    integer :: nxc, nyc, nzc, xsize1, xsize2, xsize3

    if (nclx1==1.and.xend(1)==nx) then
       xsize1=xsize(1)-1
    else
       xsize1=xsize(1)
    endif
    if (ncly1==1.and.xend(2)==ny) then
       xsize2=xsize(2)-1
    else
       xsize2=xsize(2)
    endif
    if (nclz1==1.and.xend(3)==nz) then
       xsize3=xsize(3)-1
    else
       xsize3=xsize(3)
    endif
    if (nclx1==1) then
       nxc=nxm
    else
       nxc=nx
    endif
    if (ncly1==1) then
       nyc=nym
    else
       nyc=ny
    endif
    if (nclz1==1) then
       nzc=nzm
    else
       nzc=nz
    endif

    dep=zero
    do k=1,xsize3
       do j=1,xsize2
          do i=1,xsize1
             !dep=dep+var(i,j,k)**2
             dep=dep+var(i,j,k)
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(dep,avg,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    avg=avg/(nxc*nyc*nzc)

    return

  end subroutine avg3d
end module tools
!##################################################################
subroutine stretching()

  use decomp_2d
  !use decomp_2d_poisson
  use variables
  use param
  use var
  use mpi
  use dbg_schemes, only: abs_prec, sqrt_prec, sin_prec, cos_prec, tan_prec, atan_prec

  implicit none

  real(mytype) :: yinf,den,xnum,xcx,den1,den2,den3,den4,xnum1,cst
  integer :: j

  yinf=-yly/two
  den=two*beta*yinf
  xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
  alpha=abs(xnum/den)
  xcx=one/beta/alpha
  if (alpha.ne.0.) then
     if (istret.eq.1) yp(1)=zero
     if (istret.eq.2) yp(1)=zero
     if (istret.eq.1) yeta(1)=zero
     if (istret.eq.2) yeta(1)=-half
     if (istret.eq.3) yp(1)=zero
     if (istret.eq.3) yeta(1)=-half
     do j=2,ny
        if (istret==1) yeta(j)=real(j-1,mytype)*(one/nym)
        if (istret==2) yeta(j)=real(j-1,mytype)*(one/nym)-half
        if (istret==3) yeta(j)=real(j-1,mytype)*(half/nym)-half
        den1=sqrt(alpha*beta+one)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
        den=two*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
        den4=two*alpha*beta-cos(two*pi*yeta(j))+one
        xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(two*sqrt(alpha)*sqrt(alpha*beta+one))
        if (istret==1) then
           if (yeta(j).lt.half) yp(j)=xnum1-cst-yinf
           if (yeta(j).eq.half) yp(j)=zero-yinf
           if (yeta(j).gt.half) yp(j)=xnum1+cst-yinf
        endif
        if (istret==2) then
           if (yeta(j).lt.half) yp(j)=xnum1-cst+yly
           if (yeta(j).eq.half) yp(j)=zero+yly
           if (yeta(j).gt.half) yp(j)=xnum1+cst+yly
        endif
        if (istret==3) then
           if (yeta(j).lt.half) yp(j)=(xnum1-cst+yly)*two
           if (yeta(j).eq.half) yp(j)=(zero+yly)*two
           if (yeta(j).gt.half) yp(j)=(xnum1+cst+yly)*two
        endif
     enddo
  endif
  if (alpha.eq.0.) then
     yp(1)=-1.e10
     do j=2,ny
        yeta(j)=real(j-1,mytype)*(one/ny)
        yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
     enddo
  endif
  if (alpha.ne.0.) then
     do j=1,ny
        if (istret==1) yetai(j)=(real(j,mytype)-half)*(one/nym)
        if (istret==2) yetai(j)=(real(j,mytype)-half)*(one/nym)-half
        if (istret==3) yetai(j)=(real(j,mytype)-half)*(half/nym)-half
        den1=sqrt(alpha*beta+one)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
        den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
        den4=two*alpha*beta-cos(two*pi*yetai(j))+one
        xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(two*sqrt(alpha)*sqrt(alpha*beta+one))
        if (istret==1) then
           if (yetai(j).lt.half) ypi(j)=xnum1-cst-yinf
           if (yetai(j).eq.half) ypi(j)=zero-yinf
           if (yetai(j).gt.half) ypi(j)=xnum1+cst-yinf
        endif
        if (istret==2) then
           if (yetai(j).lt.half) ypi(j)=xnum1-cst+yly
           if (yetai(j).eq.half) ypi(j)=zero+yly
           if (yetai(j).gt.half) ypi(j)=xnum1+cst+yly
        endif
        if (istret==3) then
           if (yetai(j).lt.half) ypi(j)=(xnum1-cst+yly)*two
           if (yetai(j).eq.half) ypi(j)=(zero+yly)*two
           if (yetai(j).gt.half) ypi(j)=(xnum1+cst+yly)*two
        endif
     enddo
  endif
  if (alpha.eq.0.) then
     ypi(1)=-1.e10
     do j=2,ny
        yetai(j)=real(j-1,mytype)*(one/ny)
        ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
     enddo
  endif

  !Mapping!!, metric terms
  if (istret .ne. 3) then
     do j=1,ny
        ppy(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yeta(j))*sin(pi*yeta(j)))
        pp2y(j)=ppy(j)*ppy(j)
        pp4y(j)=(-two/beta*cos(pi*yeta(j))*sin(pi*yeta(j)))
     enddo
     do j=1,ny
        ppyi(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yetai(j))*sin(pi*yetai(j)))
        pp2yi(j)=ppyi(j)*ppyi(j)
        pp4yi(j)=(-two/beta*cos(pi*yetai(j))*sin(pi*yetai(j)))
     enddo
  endif

  if (istret .eq. 3) then
     do j=1,ny
        ppy(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yeta(j))*sin(pi*yeta(j)))
        pp2y(j)=ppy(j)*ppy(j)
        pp4y(j)=(-two/beta*cos(pi*yeta(j))*sin(pi*yeta(j)))/two
     enddo
     do j=1,ny
        ppyi(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yetai(j))*sin(pi*yetai(j)))
        pp2yi(j)=ppyi(j)*ppyi(j)
        pp4yi(j)=(-two/beta*cos(pi*yetai(j))*sin(pi*yetai(j)))/two
     enddo
  endif

  !   yp(1) = 0.0
  !   yp(2) = 0.01
  !   coeff0= 1.1
  !   blender1 = 0.0
  !   blender2 = 0.0
  !   do j=3,ny
  !!      yeta(j)=(j-1.)*(1./ny)
  !!      yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
  !
  !     if (yp(j-1).LE.3.5*1.0) then
  !       dy_plus_target = 8.0
  !       !Calculate re_tau guess somewhere
  !      dy_plus_current= (yp(j-1)-yp(j-2))*85.0
  !       !dy_plus_coeff is from 1 to 0
  !       dy_plus_coeff = (dy_plus_target-dy_plus_current)/dy_plus_target
  !       coeff = coeff0**dy_plus_coeff
  !
  !       dy_plus_coeff_old1 = dy_plus_coeff   !will be required for blenders
  !     else if (yp(j-1).GE.39.0*1.0) then
  !       dy_plus_target = 10.0
  !       !Calculate re_tau guess somewhere
  !       dy_plus_current= (yp(j-1)-yp(j-2))*85.0
  !       !dy_plus_coeff is from 1 to 0
  !       dy_plus_coeff = (dy_plus_target-dy_plus_current)/dy_plus_target
  !
  !       if (blender2.LT.1.0) blender2 = blender2 + 0.1   !carry the coeff smoothly
  !       coeff = coeff0**((1.0-blender2)*dy_plus_coeff_old2+blender2*dy_plus_coeff)
  !     else
  !       dy_plus_target = 80.0
  !       !Calculate re_tau guess somewhere
  !       dy_plus_current= (yp(j-1)-yp(j-2))*85.0
  !       !dy_plus_coeff is from 1 to 0
  !       dy_plus_coeff = (dy_plus_target-dy_plus_current)/dy_plus_target
  !
  !       if (blender1.LT.1.0) blender1 = blender1 + 0.1   !carry the coeff smoothly
  !       coeff = coeff0**((1.0-blender1)*dy_plus_coeff_old1+blender1*dy_plus_coeff)
  !
  !       dy_plus_coeff_old2 = dy_plus_coeff   !will be required for blenders
  !     endif
  !     yp(j) = yp(j-1)+(yp(j-1)-yp(j-2))*coeff
  !   enddo
  !
  !   !Normalize to yly
  !   ypmax = yp(ny)
  !   yp = yp/ypmax*yly

  if (nrank == 0) then
     open(10,file='yp.dat', form='formatted')
     do j=1,ny
        write(10,*)yp(j)
     enddo
     close(10)
     open(10,file='ypi.dat', form='formatted')
     do j=1,nym
        write(10,*)ypi(j)
     enddo
     close(10)
  endif

end subroutine stretching
!##################################################################
subroutine inversion5_v1(aaa_in,eee,spI)

  use decomp_2d
  !use decomp_2d_poisson
  use variables
  use param
  use var
  use mpi
  use dbg_schemes, only: abs_prec

  implicit none

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

  complex(mytype),dimension(spI%yst(1):spI%yen(1),ny/2,spI%yst(3):spI%yen(3),5) :: aaa, aaa_in
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(2):spI%yen(2),spI%yst(3):spI%yen(3)) :: eee
  integer :: i,j,k,m,mi,jc
  integer,dimension(2) :: ja,jb
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

  real(mytype) :: tmp1,tmp2,tmp3,tmp4

  complex(mytype) :: cx
  real(mytype) :: rl, iy
  external cx, rl, iy

  aaa = aaa_in

  do i = 1, 2
     ja(i) = 4 - i
     jb(i) = 5 - i
  enddo
  do m = 1, ny/2 - 2
     do i = 1, 2
        mi = m + i
        do k = spI%yst(3), spI%yen(3)
           do j = spI%yst(1), spI%yen(1)
              if (rl(aaa(j,m,k,3)) /= zero) tmp1 = rl(aaa(j,mi,k,3-i)) / rl(aaa(j,m,k,3))
              if (iy(aaa(j,m,k,3)) /= zero) tmp2 = iy(aaa(j,mi,k,3-i)) / iy(aaa(j,m,k,3))
              sr(j,k)=cx(tmp1,tmp2)
              eee(j,mi,k)=cx(rl(eee(j,mi,k)) - tmp1 * rl(eee(j,m,k)),&
                             iy(eee(j,mi,k)) - tmp2 * iy(eee(j,m,k)))
           enddo
        enddo
        do jc = ja(i), jb(i)
           do k = spI%yst(3), spI%yen(3)
              do j = spI%yst(1), spI%yen(1)
                 aaa(j,mi,k,jc) = cx(rl(aaa(j,mi,k,jc)) - rl(sr(j,k)) * rl(aaa(j,m,k,jc+i)),&
                                     iy(aaa(j,mi,k,jc)) - iy(sr(j,k)) * iy(aaa(j,m,k,jc+i)))
              enddo
           enddo
        enddo
     enddo
  enddo

  do k = spI%yst(3), spI%yen(3)
     do j = spI%yst(1), spI%yen(1)
        if (abs_prec(rl(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp1 = rl(aaa(j,ny/2,k,2)) / rl(aaa(j,ny/2-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp2 = iy(aaa(j,ny/2,k,2)) / iy(aaa(j,ny/2-1,k,3))
        else
           tmp2 = zero
        endif
        sr(j,k) = cx(tmp1,tmp2)
        b1(j,k) = cx(rl(aaa(j,ny/2,k,3)) - tmp1 * rl(aaa(j,ny/2-1,k,4)),&
                     iy(aaa(j,ny/2,k,3)) - tmp2 * iy(aaa(j,ny/2-1,k,4)))

        if (abs_prec(rl(b1(j,k))) > epsilon) then
           tmp1 = rl(sr(j,k)) / rl(b1(j,k))
           tmp3 = rl(eee(j,ny/2,k)) / rl(b1(j,k)) - tmp1 * rl(eee(j,ny/2-1,k))
        else
           tmp1 = zero
           tmp3 = zero
        endif
        if (abs_prec(iy(b1(j,k))) > epsilon) then
           tmp2 = iy(sr(j,k)) / iy(b1(j,k))
           tmp4 = iy(eee(j,ny/2,k)) / iy(b1(j,k)) - tmp2 * iy(eee(j,ny/2-1,k))
        else
           tmp2 = zero
           tmp4 = zero
        endif
        a1(j,k) = cx(tmp1,tmp2)
        eee(j,ny/2,k) = cx(tmp3,tmp4)

        if (abs_prec(rl(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp1 = one / rl(aaa(j,ny/2-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp2 = one / iy(aaa(j,ny/2-1,k,3))
        else
           tmp2 = zero
        endif
        b1(j,k) = cx(tmp1, tmp2)
        a1(j,k) = cx(rl(aaa(j,ny/2-1,k,4)) * rl(b1(j,k)),&
                     iy(aaa(j,ny/2-1,k,4)) * iy(b1(j,k)))
        eee(j,ny/2-1,k) = cx(rl(eee(j,ny/2-1,k)) * rl(b1(j,k)) - rl(a1(j,k)) * rl(eee(j,ny/2,k)),&
                             iy(eee(j,ny/2-1,k)) * iy(b1(j,k)) - iy(a1(j,k)) * iy(eee(j,ny/2,k)))
     enddo
  enddo

  do i = ny/2 - 2, 1, -1
     do k = spI%yst(3), spI%yen(3)
        do j = spI%yst(1), spI%yen(1)
           if (abs_prec(rl(aaa(j,i,k,3))) > epsilon) then
              tmp1 = one / rl(aaa(j,i,k,3))
           else
              tmp1 = zero
           endif
           if (abs_prec(iy(aaa(j,i,k,3))) > epsilon) then
              tmp2 = one/iy(aaa(j,i,k,3))
           else
              tmp2 = zero
           endif
           sr(j,k) = cx(tmp1,tmp2)
           a1(j,k) = cx(rl(aaa(j,i,k,4)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,4)) * iy(sr(j,k)))
           b1(j,k) = cx(rl(aaa(j,i,k,5)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,5)) * iy(sr(j,k)))
           eee(j,i,k) = cx(rl(eee(j,i,k)) * rl(sr(j,k)) - rl(a1(j,k)) * rl(eee(j,i+1,k)) - rl(b1(j,k)) * rl(eee(j,i+2,k)),&
                           iy(eee(j,i,k)) * iy(sr(j,k)) - iy(a1(j,k)) * iy(eee(j,i+1,k)) - iy(b1(j,k)) * iy(eee(j,i+2,k)))
        enddo
     enddo
  enddo

  return

end subroutine inversion5_v1
!##################################################################
subroutine inversion5_v2(aaa,eee,spI)

  use decomp_2d
  !use decomp_2d_poisson
  use variables
  use param
  use var
  use MPI
  use dbg_schemes, only: abs_prec

  implicit none

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

  complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3),5) :: aaa
  complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3)) :: eee
  integer :: i,j,k,m,mi,jc
  integer,dimension(2) :: ja,jb
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

  real(mytype) :: tmp1,tmp2,tmp3,tmp4

  complex(mytype) :: cx
  real(mytype) :: rl, iy
  external cx, rl, iy

  do i = 1, 2
     ja(i) = 4 - i
     jb(i) = 5 - i
  enddo
  do m = 1, nym - 2
     do i = 1, 2
        mi = m + i
        do k = spI%yst(3), spI%yen(3)
           do j = spI%yst(1), spI%yen(1)
              if (rl(aaa(j,m,k,3)) /= zero) tmp1 = rl(aaa(j,mi,k,3-i)) / rl(aaa(j,m,k,3))
              if (iy(aaa(j,m,k,3)) /= zero) tmp2 = iy(aaa(j,mi,k,3-i)) / iy(aaa(j,m,k,3))
              sr(j,k) = cx(tmp1, tmp2)
              eee(j,mi,k) = cx(rl(eee(j,mi,k)) - tmp1 * rl(eee(j,m,k)),&
                               iy(eee(j,mi,k)) - tmp2 * iy(eee(j,m,k)))
           enddo
        enddo
        do jc = ja(i), jb(i)
           do k = spI%yst(3), spI%yen(3)
              do j = spI%yst(1), spI%yen(1)
                 aaa(j,mi,k,jc) = cx(rl(aaa(j,mi,k,jc)) - rl(sr(j,k)) * rl(aaa(j,m,k,jc+i)),&
                                     iy(aaa(j,mi,k,jc)) - iy(sr(j,k)) * iy(aaa(j,m,k,jc+i)))
              enddo
           enddo
        enddo
     enddo
  enddo
  do k = spI%yst(3), spI%yen(3)
     do j = spI%yst(1), spI%yen(1)
        if (abs_prec(rl(aaa(j,nym-1,k,3))) > epsilon) then
           tmp1 = rl(aaa(j,nym,k,2)) / rl(aaa(j,nym-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,nym-1,k,3))) > epsilon) then
           tmp2 = iy(aaa(j,nym,k,2)) / iy(aaa(j,nym-1,k,3))
        else
           tmp2 = zero
        endif
        sr(j,k) = cx(tmp1,tmp2)
        b1(j,k) = cx(rl(aaa(j,nym,k,3)) - tmp1 * rl(aaa(j,nym-1,k,4)),&
                     iy(aaa(j,nym,k,3)) - tmp2 * iy(aaa(j,nym-1,k,4)))
        if (abs_prec(rl(b1(j,k))) > epsilon) then
           tmp1 = rl(sr(j,k)) / rl(b1(j,k))
           tmp3 = rl(eee(j,nym,k)) / rl(b1(j,k)) - tmp1 * rl(eee(j,nym-1,k))
        else
           tmp1 = zero
           tmp3 = zero
        endif
        if (abs_prec(iy(b1(j,k))) > epsilon) then
           tmp2 = iy(sr(j,k)) / iy(b1(j,k))
           tmp4 = iy(eee(j,nym,k)) / iy(b1(j,k)) - tmp2 * iy(eee(j,nym-1,k))
        else
           tmp2 = zero
           tmp4 = zero
        endif
        a1(j,k) = cx(tmp1, tmp2)
        eee(j,nym,k) = cx(tmp3, tmp4)

        if (abs_prec(rl(aaa(j,nym-1,k,3))) > epsilon) then
           tmp1 = one / rl(aaa(j,nym-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs_prec(iy(aaa(j,nym-1,k,3))) > epsilon) then
           tmp2 = one / iy(aaa(j,nym-1,k,3))
        else
           tmp2 = zero
        endif
        b1(j,k) = cx(tmp1,tmp2)
        a1(j,k) = cx(rl(aaa(j,nym-1,k,4)) * rl(b1(j,k)),&
                     iy(aaa(j,nym-1,k,4)) * iy(b1(j,k)))
        eee(j,nym-1,k) = cx(rl(eee(j,nym-1,k)) * rl(b1(j,k)) - rl(a1(j,k)) * rl(eee(j,nym,k)),&
                            iy(eee(j,nym-1,k)) * iy(b1(j,k)) - iy(a1(j,k)) * iy(eee(j,nym,k)))
     enddo
  enddo

  do i = nym - 2, 1, -1
     do k = spI%yst(3), spI%yen(3)
        do j = spI%yst(1), spI%yen(1)
           if (abs_prec(rl(aaa(j,i,k,3))) > epsilon) then
              tmp1 = one / rl(aaa(j,i,k,3))
           else
              tmp1 = zero
           endif
           if (abs_prec(iy(aaa(j,i,k,3))) > epsilon) then
              tmp2 = one / iy(aaa(j,i,k,3))
           else
              tmp2 = zero
           endif
           sr(j,k) = cx(tmp1,tmp2)
           a1(j,k) = cx(rl(aaa(j,i,k,4)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,4)) * iy(sr(j,k)))
           b1(j,k) = cx(rl(aaa(j,i,k,5)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,5)) * iy(sr(j,k)))
           eee(j,i,k) = cx(rl(eee(j,i,k)) * rl(sr(j,k)) - rl(a1(j,k)) * rl(eee(j,i+1,k)) -rl(b1(j,k)) * rl(eee(j,i+2,k)),&
                           iy(eee(j,i,k)) * iy(sr(j,k)) - iy(a1(j,k)) * iy(eee(j,i+1,k)) -iy(b1(j,k)) * iy(eee(j,i+2,k)))
        enddo
     enddo
  enddo

  return

end subroutine inversion5_v2
!##################################################################
function rl(complexnumber)

  use param

  implicit none

  real(mytype) :: rl
  complex(mytype) :: complexnumber

  rl = real(complexnumber, kind=mytype)

end function rl
!##################################################################
function iy(complexnumber)

  use param

  implicit none

  real(mytype) :: iy
  complex(mytype) :: complexnumber

  iy = aimag(complexnumber)

end function iy
!##################################################################
function cx(realpart,imaginarypart)

  use param

  implicit none

  complex(mytype) :: cx
  real(mytype) :: realpart, imaginarypart

  cx = cmplx(realpart, imaginarypart, kind=mytype)

end function cx
!##################################################################
subroutine calc_mweight(mweight, phi, xlen, ylen, zlen)

  use decomp_2d
  use param, only : zero, one
  use param, only : massfrac, mol_weight
  use var, only : numscalar

  implicit none

  integer, intent(in) :: xlen, ylen, zlen
  real(mytype), intent(in), dimension(xlen, ylen, zlen, numscalar) :: phi

  !! LOCALS
  real(mytype), dimension(xlen, ylen, zlen) :: mweight
  integer :: is

  mweight(:,:,:) = zero
  do is = 1, numscalar
     if (massfrac(is)) then
        mweight(:,:,:) = mweight(:,:,:) + phi(:,:,:,is) / mol_weight(is)
     endif
  enddo
  mweight(:,:,:) = one / mweight(:,:,:)

endsubroutine calc_mweight
!##################################################################
function r8_random ( s1, s2, s3 )

!*****************************************************************************80
!
!! R8_RANDOM returns a pseudorandom number between 0 and 1.
!
!  Discussion:
!
!    This function returns a pseudo-random number rectangularly distributed
!    between 0 and 1.   The cycle length is 6.95E+12.  (See page 123
!    of Applied Statistics (1984) volume 33), not as claimed in the
!    original article.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2008
!
!  Author:
!
!    FORTRAN77 original version by Brian Wichman, David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Brian Wichman, David Hill,
!    Algorithm AS 183: An Efficient and Portable Pseudo-Random
!    Number Generator,
!    Applied Statistics,
!    Volume 31, Number 2, 1982, pages 188-190.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) S1, S2, S3, three values used as the
!    seed for the sequence.  These values should be positive
!    integers between 1 and 30,000.
!
!    Output, real ( kind = 8 ) R8_RANDOM, the next value in the sequence.
!
  implicit none

  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) s3
  real ( kind = 8 ) r8_random

  s1 = mod ( 171 * s1, 30269 )
  s2 = mod ( 172 * s2, 30307 )
  s3 = mod ( 170 * s3, 30323 )

  r8_random = mod ( real ( s1, kind = 8 ) / 30269.0D+00 &
                  + real ( s2, kind = 8 ) / 30307.0D+00 &
                  + real ( s3, kind = 8 ) / 30323.0D+00, 1.0D+00 )

  return
end
!##################################################################
function return_30k(x) result(y)

  integer ( kind = 4 ), intent(in) :: x
  integer ( kind = 4 )             :: y
  integer ( kind = 4 ), parameter  :: xmax = 30000

  y = iabs(x) - int(iabs(x)/xmax)*xmax
end function return_30k
!##################################################################
function r8_uni ( s1, s2 )

!*****************************************************************************80
!
!! R8_UNI returns a pseudorandom number between 0 and 1.
!
!  Discussion:
!
!    This function generates uniformly distributed pseudorandom numbers
!    between 0 and 1, using the 32-bit generator from figure 3 of
!    the article by L'Ecuyer.
!
!    The cycle length is claimed to be 2.30584E+18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2008
!
!  Author:
!
!    Original Pascal original version by Pierre L'Ecuyer
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Pierre LEcuyer,
!    Efficient and Portable Combined Random Number Generators,
!    Communications of the ACM,
!    Volume 31, Number 6, June 1988, pages 742-751.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) S1, S2, two values used as the
!    seed for the sequence.  On first call, the user should initialize
!    S1 to a value between 1 and 2147483562;  S2 should be initialized
!    to a value between 1 and 2147483398.
!
!    Output, real ( kind = 8 ) R8_UNI, the next value in the sequence.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uni
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) z

  k = s1 / 53668
  s1 = 40014 * ( s1 - k * 53668 ) - k * 12211
  if ( s1 < 0 ) then
    s1 = s1 + 2147483563
  end if

  k = s2 / 52774
  s2 = 40692 * ( s2 - k * 52774 ) - k * 3791
  if ( s2 < 0 ) then
    s2 = s2 + 2147483399
  end if

  z = s1 - s2
  if ( z < 1 ) then
    z = z + 2147483562
  end if

  r8_uni = real ( z, kind = 8 ) / 2147483563.0D+00

  return
end
!##################################################################
subroutine test_min_max(name,text,array_tmp,i_size_array_tmp)

  use param
  use variables
  use decomp_2d
  use MPI

  implicit none

  integer :: ierror, i, i_size_array_tmp
  real(mytype) :: max_tmp, min_tmp, tot_tmp, max_tot, min_tot, tot_tot
  real(mytype), dimension(i_size_array_tmp) :: array_tmp
  character(len=5) :: name
  character(len=15) :: text

  max_tmp=-0.000000000000000001_mytype
  tot_tmp=0._mytype
  min_tmp=+1000000000000000000._mytype
  do i=1,size(array_tmp)
    max_tmp=max(max_tmp,array_tmp(i))
    tot_tmp=tot_tmp + array_tmp(i)
    min_tmp=min(min_tmp,array_tmp(i))
  enddo
  call MPI_ALLREDUCE(max_tmp,max_tot,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
  call MPI_ALLREDUCE(min_tmp,min_tot,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierror)
  call MPI_ALLREDUCE(tot_tmp,tot_tot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
  if (nrank == 0) then
     write(*,*) " "
     write(*,*) trim(text)//' Max ',name,max_tot
     write(*,*) trim(text)//' Tot ',name,tot_tot
     write(*,*) trim(text)//' Min ',name,min_tot
     write(*,*) " "
     call flush(6)
  endif

  return
end subroutine test_min_max
