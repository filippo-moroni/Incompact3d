!----------------------------------------------------------!
!         This module is used to run temporal BLs          !
!                       simulations.                       !
!      Adapted from original Incompact3d file (v4.0)       !
!               for channel flow simulations.              !
!----------------------------------------------------------!

module temporal_bl

  use decomp_2d
  use variables
  use param

  implicit none

  integer :: FS
  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE ! All functions/subroutines private by default
  PUBLIC :: 

contains
  !############################################################################
  subroutine init_temporal_bl (ux1,uy1,uz1,ep1,phi1)

    use decomp_2d
    use decomp_2d_io
    use variables
    use param
    use MPI
    use dbg_schemes, only: exp_prec, abs_prec, sqrt_prec
    

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y,r,um,r3,x,z,h,ct
    real(mytype) :: cx0,cy0,cz0,hg,lg
    real(mytype) :: ftent
    integer :: k,j,i,fh,ierror,ii,is,it,code, jj

    !integer, dimension (:), allocatable :: seed
    !real(mytype), dimension(:,:,:), allocatable :: urand

    integer :: xshift, yshift, zshift

    integer ( kind = 4 ) :: seed1, seed2, seed3, seed11, seed22, seed33
    integer ( kind = 4 ) :: return_30k
    integer ( kind = 4 ), parameter :: nsemini = 1000 ! For the moment we fix it but after this can go in the input file
    real(mytype), dimension(3,nsemini) :: eddy, posvor
    real(mytype)     :: volsemini, rrand, ddx, ddy, ddz, lsem, upr, vpr, wpr, rrand1
    real(mytype), dimension(3) :: dim_min, dim_max
    real( kind = 8 ) :: r8_random
    external r8_random, return_30k

    ! to do: add tanh profile + noise
    
    ! to be modified
    if (iscalar==1) then
       if (nrank==0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime == ilast)) then
          write(*,*) 'Imposing linear temperature profile'
       end if
       do k=1,xsize(3)
          do j=1,xsize(2)
             if (istret==0) y=real(j+xstart(2)-2,mytype)*dy
             if (istret/=0) y=yp(j+xstart(2)-1)
             do i=1,xsize(1)
                phi1(i,j,k,:) = one - y/yly
             enddo
          enddo
       enddo

       phi1(:,:,:,:) = zero !change as much as you want
       if ((nclyS1 == 2).and.(xstart(2) == 1)) then
         !! Generate a hot patch on bottom boundary
         phi1(:,1,:,:) = one
       endif
       if ((nclySn == 2).and.(xend(2) == ny)) then
         phi1(:,xsize(2),:,:) = zero
       endif
    endif


    ux1=zero
    uy1=zero
    uz1=zero
    
    byx1=zero;
    byy1=zero;
    byz1=zero

    ! Traditional init to turbulent flows using random numbers + lam profile
    
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)
       !modulation of the random noise + initial velocity profile
       do k=1,xsize(3)
          do j=1,xsize(2)
             if (istret==0) y=real(j+xstart(2)-1-1,mytype)*dy-yly*half
             if (istret/=0) y=yp(j+xstart(2)-1)-yly*half
             um=exp_prec(-zptwo*y*y)
             do i=1,xsize(1)
   
                   ux1(i,j,k)=init_noise*um*(two*ux1(i,j,k)-one)+one-y*y
                   uy1(i,j,k)=init_noise*um*(two*uy1(i,j,k)-one)
                   uz1(i,j,k)=init_noise*um*(two*uz1(i,j,k)-one)

             enddo
          enddo
       enddo
       
    return
  end subroutine init_temporal_bl
  !############################################################################
  !############################################################################
  subroutine boundary_conditions_ttbl (ux,uy,uz,phi)

    use param
    use variables
    use decomp_2d

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    
    ! check if BCs needs to be explicitly declared
    
  end subroutine boundary_conditions_ttbl
  
end module temporal_bl



