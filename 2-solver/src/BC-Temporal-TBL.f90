!----------------------------------------------------------!
!         This module is used to run temporal TBLs         !
!                       simulations.                       !
!      Adapted from original Incompact3d file (v4.0)       !
!               for channel flow simulations.              !
!----------------------------------------------------------!

module temporal_tbl

  use decomp_2d
  use variables
  use param

  implicit none

  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10)  ! new line character

  private   ! All functions/subroutines private by default
  public :: init_temporal_tbl,        &
            boundary_conditions_ttbl, &
            postprocess_ttbl,         &
            visu_ttbl_init,           &
            visu_ttbl

contains
  !---------------------------------------------------------------------------!
  ! Initial condition for a temporal turbulent boundary layer (TTBL) as
  ! Kozul et al. (2016).
  !---------------------------------------------------------------------------!
  subroutine init_temporal_tbl (ux1,uy1,uz1,phi1)

    use tools,       only : apply_spatial_filter 
    use ibm_param
         
    implicit none
       
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1

    real(mytype) :: y           ! y coordinate of grid points
    real(mytype) :: theta_sl    ! momentum thickness of the initial shear layer
    real(mytype) :: um          ! initial mean streamwise velocity profile 
    real(mytype) :: phim        ! initial mean scalar profile
    real(mytype) :: diff        ! difference between wall and mean velocities
    real(mytype) :: mg          ! mean gradient of the initial velocity profile
    real(mytype) :: sh_vel_ic   ! shear velocity of the initial velocity profile
    real(mytype) :: delta_nu_ic ! viscous length of the initial velocity profile    
            
    integer      :: ii,code     ! for random numbers generation 
    integer      :: i,j,k
        
    ! Momentum thickness calculation 
    theta_sl = 54.0_mytype*xnu/uwall
    
    ! Initial mean gradient
    mg = - (uwall / (four * theta_sl)) * (one / cosh(twd / two / theta_sl))**2
    
    ! Initial shear velocity
    sh_vel_ic = sqrt(xnu * abs(mg))
    
    ! Initial viscous length
    delta_nu_ic = xnu / sh_vel_ic

    ! Initialize velocity fields
    ux1 = zero 
    uy1 = zero 
    uz1 = zero   
    
    ! Initialize BC of the bottom wall
    byx1 = uwall 
    byy1 = zero
    byz1 = zero
    
    ! Initialize the spanwise velocity at the bottom wall
    span_vel = zero
 
    !--- Initialization as Kozul et al. (JFM, 2016) (tanh + noise) ---!
    
    ! Noise (random numbers from 0 to 1)
    call system_clock(count=code)
    
    ! Always the same noise if selected from the user, otherwise pseudo-random
    if (iin.eq.0) code=0 
    
    call random_seed(size = ii)
    call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

    call random_number(ux1)
    call random_number(uy1)
    call random_number(uz1)
        
    ! Initialize the filter 
    call filter(C_filter)  ! the argument is alpha £ [-0.5, 0.5] (0.5: no filtering, -0.5: maximum filtering)

    ! Filtering (no IBM is considered, we employ the subroutine from 'tools' module)
    call apply_spatial_filter(ux1,uy1,uz1)
                          
    ! Noise superimposed to the tanh velocity profile
    do k=1,xsize(3)
       do j=1,xsize(2)
       
       ! y-coordinate calculation
       if (istret==0) y=real(j+xstart(2)-1-1,mytype)*dy
       if (istret/=0) y=yp(j+xstart(2)-1)
                    
       ! Initial streamwise velocity profile
       um = uwall*(half + half*(tanh((twd/two/theta_sl)*(one - y/twd))))
           
       ! Difference between wall and mean velocities
       diff = uwall - um
                    
             ! Add noise near the wall, excluding first grid points if requested 
             if (diff < uln*uwall .and. y/delta_nu_ic > lln) then            
                               
                do i=1,xsize(1)
                                
                ! Rescaling the noise with a percentage of the wall velocity and center it with respect to zero
                ux1(i,j,k) = (ux1(i,j,k)*two - one)*init_noise*uwall
                uy1(i,j,k) = (uy1(i,j,k)*two - one)*init_noise*uwall
                uz1(i,j,k) = (uz1(i,j,k)*two - one)*init_noise*uwall
                 
                ux1(i,j,k) = ux1(i,j,k) + um 
                enddo
             
             ! Area with no noise, only mean velocity profile
             else 
                do i=1,xsize(1)
                  
                  ! No noise
                  ux1(i,j,k) = zero
                  uy1(i,j,k) = zero
                  uz1(i,j,k) = zero
                 
                  ux1(i,j,k)= ux1(i,j,k) + um 
                enddo
             end if 
       enddo
    enddo
    
    !--- Initialization of scalar fields ---!  
    
    ! Initialize only if scalar fields are present
    if(iscalar==1) then
    
    phi1=zero
    
       do k=1,xsize(3)
          do j=1,xsize(2)
             if (istret==0) y=real(j+xstart(2)-2,mytype)*dy
             if (istret/=0) y=yp(j+xstart(2)-1)
             
             ! Initial scalar profile (same thickness of the velocity BL, change if Pr =/ 1)
             phim = phiwall*(half + half*(tanh((twd/two/theta_sl)*(one - y/twd))))
            
             ! Only mean profile, no noise            
             do i=1,xsize(1)
                phi1(i,j,k,:) = phim
             enddo
             
          enddo
       enddo
    end if   
       
    return
  end subroutine init_temporal_tbl
  !---------------------------------------------------------------------------!
  ! Boundary conditions for a temporal turbulent boundary layer (TTBL) with
  ! possibility of scalar field and spanwise wall oscillations.
  !---------------------------------------------------------------------------!
  subroutine boundary_conditions_ttbl(phi)

    implicit none
    
    integer :: i,j,k
    
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    
    ! Top boundary Neumann BCs (e.g. free-slip for velocity) 
    ! do not need to be explicitly re-defined (both for velocity and scalar fields).
     
    ! If spanwise motion is present, span_vel is calculated at each sub-time step,
    ! otherwise it is zero.
     
    ! Bottom boundary (Dirichlet, imposed velocity of the wall)
    if (ncly1 == 2) then
      do k = 1, xsize(3)
        do i = 1, xsize(1)
          byx1(i,k) = uwall 
          byy1(i,k) = zero
          byz1(i,k) = span_vel
        enddo
      enddo
    endif
        
    ! Top boundary (Dirichlet, no-slip condition) 
    if (nclyn == 2) then
      do k = 1, xsize(3)
        do i = 1, xsize(1)
          byxn(i,k) = zero  
          byyn(i,k) = zero
          byzn(i,k) = zero
        enddo
      enddo
    endif
    
    ! Scalar fields
    if(iscalar==1) then
    
       ! Bottom boundary (Dirichlet, imposed scalar value at the wall)   
       if ((nclyS1 == 2).and.(xstart(2) == 1)) then
         phi(:,1,:,:) = phiwall
       endif
       
       ! Top boundary (Dirichlet, imposed scalar value at the freestream)
       if ((nclySn == 2).and.(xend(2) == ny)) then
         phi(:,xsize(2),:,:) = zero
       endif
       
    endif
      
  end subroutine boundary_conditions_ttbl
  !---------------------------------------------------------------------------!
  subroutine postprocess_ttbl(ux1,uy1,uz1,pp3,phi1,ep1)

    use var, only : nzmsize

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

  end subroutine postprocess_ttbl
  !---------------------------------------------------------------------------!
  subroutine visu_ttbl_init(visu_initialised)

    use decomp_2d_io, only : decomp_2d_register_variable
    use visu,         only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    ! Q-criterion
    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)
    
    visu_initialised = .true.
    
  end subroutine visu_ttbl_init
  !---------------------------------------------------------------------------!
  subroutine visu_ttbl(ux1, uy1, uz1, pp3, phi1, ep1, num)

    use var,  only : ux2, uy2, uz2, ux3, uy3, uz3
    use var,  only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    use var,  only : ta2,tb2,tc2,td2,te2,tf2,di2,ta3,tb3,tc3,td3,te3,tf3,di3
    use var,  only : nzmsize
    use visu, only : write_field
    
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    character(len=32), intent(in) :: num

    ! Write Q-criterion as an example of post processing

    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uy1,uy2)
      call transpose_x_to_y(uz1,uz2)
      call transpose_y_to_z(ux2,ux3)
      call transpose_y_to_z(uy2,uy3)
      call transpose_y_to_z(uz2,uz3)
      sync_vel_needed = .false.
    endif

    ! x-derivatives
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,ubcx)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    ! y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tb2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,ubcy)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    ! z-derivatives
    call derz (ta3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
    call derz (tb3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
    call derz (tc3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,ubcz)
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

    !Q=-0.5*(ta1**2+te1**2+di1**2)-td1*tb1-tg1*tc1-th1*tf1
    di1 = zero
    di1(:,:,:) = - half*(ta1(:,:,:)**2 + te1(:,:,:)**2 + ti1(:,:,:)**2) &
                 - td1(:,:,:) * tb1(:,:,:) &
                 - tg1(:,:,:) * tc1(:,:,:) &
                 - th1(:,:,:) * tf1(:,:,:)
                     
    call write_field(di1, ".", "critq", trim(num), flush = .true.)  ! Reusing temporary array, force flush

  end subroutine visu_ttbl
  
end module temporal_tbl



