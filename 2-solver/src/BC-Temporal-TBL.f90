!----------------------------------------------------------!
!         This module is used to run temporal BLs          !
!                       simulations.                       !
!      Adapted from original Incompact3d files (v4.0)      !
!               for channel flow simulations.              !
!----------------------------------------------------------!

module temporal_tbl

  use decomp_2d
  use variables
  use param

  implicit none

  character(len=100) :: fileformat
  character(len=1),parameter :: NL=char(10) !new line character

  PRIVATE   ! All functions/subroutines private by default
  PUBLIC :: init_temporal_tbl,boundary_conditions_ttbl, postprocess_ttbl, &
            visu_ttbl_init, visu_ttbl

contains
  !############################################################################
  subroutine init_temporal_tbl (ux1,uy1,uz1)

    use decomp_2d
    use decomp_2d_io
    use variables
    use param
    use MPI
    use dbg_schemes, only : tanh_prec,cosh_prec,sqrt_prec,abs_prec
    use var,         only : ta1,tb1,tc1,di1
    use var,         only : ta2,tb2,tc2,di2
    use var,         only : ta3,tb3,tc3,di3
    use ibm_param    
        
    implicit none
       
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1

    real(mytype) :: y         ! y-coordinate of grid points
    real(mytype) :: theta_sl  ! momentum thickness of the initial shear layer
    real(mytype) :: um        ! mean streamwise initial velocity profile
    real(mytype) :: diff      ! difference between wall and mean velocities
    real(mytype) :: mg        ! mean gradient of the initial velocity profile
    real(mytype) :: sh_vel    ! shear velocity of the initial velocity profile
    real(mytype) :: delta_nu  ! viscous length of the initial velocity profile
    
    real(mytype) :: lind 
    
    ! Filtered velocity fields, with pencils for trasposition for filtering
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1f, uy1f, uz1f
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)) :: ux2f, uy2f, uz2f
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)) :: ux3f, uy3f, uz3f
        
    ! For random numbers generation
    integer :: ii,code  
    integer :: i,j,k
    
    ! Momentum thickness calculation
    theta_sl = 54.0*xnu/uwall
    
    ! Initial mean gradient
    mg = - (uwall / (4.0 * theta_sl)) * (1.0 / cosh_prec(twd / 2.0 / theta_sl))**2
    
    ! Initial shear velocity
    sh_vel = sqrt_prec(xnu * abs_prec(mg))
    
    ! Initial viscous length
    delta_nu = xnu / sh_vel

    ! Initialize velocity fields and boundaries
    ux1=zero 
    uy1=zero 
    uz1=zero   
 
    ! Initialization as Kozul et al. (JFM, 2016) (tanh + noise)
    
    ! Noise (random numbers from 0 to 1)
    call system_clock(count=code)
    
    ! Uncomment the following line to obtain always the same noise
    ! if (iin.eq.2) code=0 
    
    call random_seed(size = ii)
    call random_seed(put = code+63946*(nrank+1)*(/ (i - 1, i = 1, ii) /))

    call random_number(ux1)
    call random_number(uy1)
    call random_number(uz1)
    
    ! No ibm considered (see dynsmag subroutine in les_models.f90)
    ta1 = ux1
    tb1 = uy1
    tc1 = uz1
    
    ! Initialize the filter 
    call filter(zero)  ! the argument is alpha £ [-0.5, 0.5] (0.5: no filtering, -0.5: maximum filtering)

    ! Low-pass filtering of the white noise, x-direction
    call filx(ux1f,ta1,di1,fisx,fiffx ,fifsx ,fifwx ,xsize(1),xsize(2),xsize(3),0,ubcx) !ux1
    call filx(uy1f,tb1,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcy) !uy1
    call filx(uz1f,tc1,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcz) !uz1
    
    ! Transposition along y
    call transpose_x_to_y(ux1f,ta2)
    call transpose_x_to_y(uy1f,tb2)
    call transpose_x_to_y(uz1f,tc2)
    
    ! Filtering along y   
    call fily(ux2f,ta2,di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcx) 
    call fily(uy2f,tb2,di2,fisy,fiffy ,fifsy ,fifwy ,ysize(1),ysize(2),ysize(3),0,ubcy) 
    call fily(uz2f,tc2,di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcz) 
    
    ! Transposition along z
    call transpose_y_to_z(ux2f,ta3)
    call transpose_y_to_z(uy2f,tb3)
    call transpose_y_to_z(uz2f,tc3)
    
    ! Filtering along z
    call filz(ux3f,ta3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcx) 
    call filz(uy3f,tb3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcy) 
    call filz(uz3f,tc3,di3,fisz,fiffz ,fifsz ,fifwz ,zsize(1),zsize(2),zsize(3),0,ubcz) 
    
    ! Back to x pencils
    call transpose_z_to_y(ux3f,ux2f)
    call transpose_z_to_y(uy3f,uy2f)
    call transpose_z_to_y(uz3f,uz2f)
    
    call transpose_y_to_x(ux2f,ux1f)
    call transpose_y_to_x(uy2f,uy1f)
    call transpose_y_to_x(uz2f,uz1f)
                          
    ! Noise superimposed to the tanh velocity profile
    do k=1,xsize(3)
       do j=1,xsize(2)
       
       ! y-coordinate calculation
       if (istret==0) y=real(j+xstart(2)-1-1,mytype)*dy
       if (istret/=0) y=yp(j+xstart(2)-1)
                    
       ! Initial streamwise velocity profile
       um = uwall*(half + half*(tanh_prec((twd/two/theta_sl)*(one - y/twd))))
             
       ! Difference between wall and mean velocities
       diff = uwall - um
                    
             ! Add noise in an intermediate region near the wall but excluding first grid points
             if (diff < uln*uwall .and. y/delta_nu > lln) then            
                do i=1,xsize(1)
                
                ! Rescaling the noise with a percentage of the wall velocity and center it with respect to zero
                ux1(i,j,k) = (ux1f(i,j,k)*two - one)*init_noise*uwall
                uy1(i,j,k) = (uy1f(i,j,k)*two - one)*init_noise*uwall
                uz1(i,j,k) = (uz1f(i,j,k)*two - one)*init_noise*uwall
                 
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
       
    return
  end subroutine init_temporal_tbl
  !############################################################################
  subroutine boundary_conditions_ttbl()

    use param
    use variables
    use decomp_2d

    implicit none
    
    integer :: i,j,k
       
    ! Bottom boundary (Dirichlet, imposed velocity of the wall)
    if (ncly1 == 2) then
      do k = 1, xsize(3)
        do i = 1, xsize(1)
          byx1(i,k) = uwall 
          byy1(i,k) = zero
          byz1(i,k) = zero
        enddo
      enddo
    endif
    
    ! Top Neumann BC (free-slip) (not needed to be explicitly re-defined)
    
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
       
  end subroutine boundary_conditions_ttbl
  !############################################################################
  subroutine postprocess_ttbl(ux1,uy1,uz1,pp3,phi1,ep1)

    use var, only : nzmsize

    implicit none

    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1, ep1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize,npress) :: pp3

  end subroutine postprocess_ttbl
  !############################################################################
  subroutine visu_ttbl_init(visu_initialised)

    use decomp_2d,    only : mytype
    use decomp_2d_io, only : decomp_2d_register_variable
    use visu,         only : io_name, output2D
    
    implicit none

    logical, intent(out) :: visu_initialised

    call decomp_2d_register_variable(io_name, "critq", 1, 0, output2D, mytype)

    visu_initialised = .true.
    
  end subroutine visu_ttbl_init
  !############################################################################
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
  !############################################################################
   
end module temporal_tbl



