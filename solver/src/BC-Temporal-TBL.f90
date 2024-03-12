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

  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE   ! All functions/subroutines private by default
  PUBLIC :: init_temporal_tbl,boundary_conditions_ttbl

contains
  !############################################################################
  subroutine init_temporal_tbl (ux1,uy1,uz1)

    use decomp_2d
    use decomp_2d_io
    use variables
    use param
    use MPI
    use dbg_schemes, only: tanh_prec
    
    implicit none
       
    real(mytype),intent(inout),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1

    real(mytype) :: y
    real(mytype) :: theta_sl  ! momentum thickness of the initial shear layer
    real(mytype) :: um        ! mean streamwise initial velocity profile
    real(mytype) :: diff      ! difference between wall and mean velocities
    
    
    ! For random numbers generation (check better what 'seed' means)
    integer :: ii,code  
    integer :: i,j,k
    
    ! Momentum thickness calculation
    theta_sl = 54.0*xnu/uwall

    ! Initialize velocity fields
    ux1=zero
    uy1=zero
    uz1=zero
 
    ! Initialization as Kozul et al. (JFM, 2016) (tanh + noise)
    
       ! Noise (random numbers from 0 to 1)
       call system_clock(count=code)
       if (iin.eq.2) code=0
       call random_seed(size = ii)
       call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

       call random_number(ux1)
       call random_number(uy1)
       call random_number(uz1)
             
       ! Noise superimposed to the tanh velocity profile
       do k=1,xsize(3)
          do j=1,xsize(2)
          
             if (istret==0) y=real(j+xstart(2)-1-1,mytype)*dy
             if (istret/=0) y=yp(j+xstart(2)-1)
             
             ! Initial streamwise velocity profile
             um = uwall*(0.5 + 0.5*(tanh_prec(twd/2*theta_sl)*(1.0 - y/twd)))
             
             ! Difference between wall and mean velocities
             diff = uwall - um
             
             ! Area near the wall, we add noise to all velocity components
             if (diff < noise_loc*uwall) then
             
             do i=1,xsize(1)
                  ! Rescaling the noise with a percentage of the wall velocity
                  ux1 = ux1*init_noise*uwall
                  uy1 = uy1*init_noise*uwall
                  uz1 = uz1*init_noise*uwall
                 
                  ux1(i,j,k)= ux1(i,j,k) + um 
             enddo
             
             ! Area away from the wall, no noise, only mean velocity profile
             else 
             do i=1,xsize(1)
   
                  ux1 = zero
                  uy1 = zero
                  uz1 = zero
                 
                  ux1(i,j,k)= ux1(i,j,k) + um 
             enddo  
          enddo
       enddo
       
    return
  end subroutine init_temporal_tbl
  !############################################################################
  !############################################################################
  subroutine boundary_conditions_ttbl()

    use param
    use variables
    use decomp_2d

    implicit none
    
    ! Bottom boundary (other BCs should not be defined explicitly) (to check)
    if (ncly1 == 2) then
      do k = 1, xsize(3)
        do i = 1, xsize(1)
          byx1(i,k) = uwall  
          byy1(i,k) = zero
          byz1(i,k) = zero
        enddo
      enddo
    endif
    
  end subroutine boundary_conditions_ttbl
  
end module temporal_bl



