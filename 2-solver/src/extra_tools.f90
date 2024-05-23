!----------------------------------------------------------!
!         This module is used to store useful              !
!           subroutines for general purpose.               !
!----------------------------------------------------------!

module extra_tools

  implicit none
  
  private 
  
  public :: print_cf,                   &
            calculate_shear_velocity,   &
            spanwise_wall_oscillations, &
            calculate_ubulk,            &
            update_time_int_coeff

contains
  
  !---------------------------------------------------------------------------!
  ! Write skin friction coefficients, shear velocity,
  ! viscous time unit and time unit and stores
  ! them in a .txt file (used for TTBL and Channel).
  !---------------------------------------------------------------------------!
  subroutine print_cf(ux,uz)
  
  use param
  use decomp_2d
  use dbg_schemes, only : sqrt_prec
      
  implicit none
  
  integer :: iunit
  logical :: exists
  
  ! Inputs 
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux, uz
  
  ! TTBL
  if(itype .eq. itype_ttbl) then
      ! Shear velocity bottom wall
      call calculate_shear_velocity(ux,uz,sh_vel,sh_velx,sh_velz)
      
      ! Create or open a file to store sh_vel, cf coefficients, viscous time and time unit
      if(nrank .eq. 0) then
          ! Calculate friction coefficients
          fric_coeff  = two * ((sh_vel  / uwall)**2)
          fric_coeffx = two * ((sh_velx / uwall)**2)
          fric_coeffz = two * ((sh_velz / uwall)**2)
          
          ! Calculate viscous time unit
          if(iswitch_wo .eq. 1) then
              ! Moving walls
              t_viscous = xnu / (sh_velx**2)
          else
              ! Fixed walls
              t_viscous = xnu / (sh_vel**2)
          end if
          
          inquire(file="cf_history.txt", exist=exists)
          if (exists) then
              open(newunit=iunit, file="cf_history.txt", status="old", position="append", action="write")
              
              write(iunit, '(F12.6,A,F16.10,A,F16.10,A,F16.10,A,F12.6,A,F12.4,A,I12)') &
              sh_vel, ',', fric_coeff, ',', fric_coeffx, ',', fric_coeffz, ',', t_viscous, ',', t, ',', itime
          else
              open(newunit=iunit, file="cf_history.txt", status="new", action="write")
              ! Header
              write(iunit, '(A12,  A,A16,   A,A16,   A,A16,   A,A12,  A,A12,  A,A12)') &
              'sh_vel', ',', 'cf,tot', ',', 'cf,x', ',', 'cf,z', ',', 't_nu', ',', 'T', ',', 'ts'          
              
              write(iunit, '(F12.6,A,F16.10,A,F16.10,A,F16.10,A,F12.6,A,F12.4,A,I12)') &
              sh_vel, ',', fric_coeff, ',', fric_coeffx, ',', fric_coeffz, ',', t_viscous, ',', t, ',', itime
          end if
              close(iunit)
      end if
 
  ! Channel
  else if(itype .eq. itype_channel) then
      ! Shear velocity bottom wall
      call calculate_shear_velocity(ux,uz,sh_vel,sh_velx,sh_velz)
      ! Bulk velocity
      call calculate_ubulk(ux,ubulk)
    
      ! Create or open a file to store sh_vel, cf coefficients, viscous time, time unit and bulk velocity
      if(nrank .eq. 0) then
          ! Calculate friction coefficients
          fric_coeff  = two * ((sh_vel  / ubulk)**2)
          fric_coeffx = two * ((sh_velx / ubulk)**2)
          fric_coeffz = two * ((sh_velz / ubulk)**2)
          
          ! Calculate viscous time unit
          if(iswitch_wo .eq. 1) then
              ! Moving walls
              t_viscous = xnu / (sh_velx**2)
          else
              ! Fixed walls
              t_viscous = xnu / (sh_vel**2)
          end if
          
          inquire(file="cf_history.txt", exist=exists)
          if (exists) then
              open(newunit=iunit, file="cf_history.txt", status="old", position="append", action="write")
              
              write(iunit, '(F12.6,A,F16.10,A,F16.10,A,F16.10,A,F12.6,A,F12.4,A,F12.4,A,I12)') &
              sh_vel, ',', fric_coeff, ',', fric_coeffx, ',', fric_coeffz, ',', t_viscous, ',', t, ',', ubulk, ',', itime
          else
              open(newunit=iunit, file="cf_history.txt", status="new", action="write")
              ! Header
              write(iunit, '(A12,  A,A16,   A,A16,   A,A16,   A,A12,  A,A12,  A,A12,  A,A12)') &
              'sh_vel', ',', 'cf,tot', ',', 'cf,x', ',', 'cf,z', ',', 't_nu', ',', 'T', ',', 'Ubulk', ',', 'ts'          
              
              write(iunit, '(F12.6,A,F16.10,A,F16.10,A,F16.10,A,F12.6,A,F12.4,A,F12.4,A,I12)') &
              sh_vel, ',', fric_coeff, ',', fric_coeffx, ',', fric_coeffz, ',', t_viscous, ',', t, ',', ubulk, ',', itime
          end if
              close(iunit)
      end if  
  
  end if ! closing of the if condition for different flow cases
             
  end subroutine print_cf
  
  !---------------------------------------------------------------------------!
  ! Calculate total shear velocity and its x and z components
  !
  ! - Used in BC-Temporal-TBL and in BC-Channel-flow
  !   for the spanwise wall oscillations. 
  ! - Used to print cf coefficients and shear velocity to an overall 
  !   .txt file for time evolution check.
  !---------------------------------------------------------------------------!
  subroutine calculate_shear_velocity(ux,uz,sh_vel,sh_velx,sh_velz)
    
    use var,         only : ux2, uz2     
    use ibm_param,   only : ubcx,ubcy,ubcz
    use dbg_schemes, only : sqrt_prec, abs_prec
    
    use var,         only : ta2,tc2,di2
    use ibm_param,   only : ubcx,ubcz
    
    use MPI
    use decomp_2d,   only : mytype, real_type, nrank
    use decomp_2d,   only : xsize, ysize
    use decomp_2d,   only : transpose_x_to_y
    
    use param,       only : zero, xnu
    use variables
    
    implicit none
    
    ! Inputs
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux, uz

    ! Outputs
    real(mytype), intent(out) :: sh_vel  ! Total shear velocity 
    real(mytype), intent(out) :: sh_velx ! Shear velocity along x 
    real(mytype), intent(out) :: sh_velz ! Shear velocity along z
    
    ! Work variables
    real(mytype) :: mean_gw    ! Mean total parallel gradient at each processor
    real(mytype) :: mean_gwx   ! Mean gradient direction x at each processor
    real(mytype) :: mean_gwz   ! Mean gradient direction z at each processor
    real(mytype) :: den        ! Denominator of the divisions
       
    integer      :: ierr         
    integer      :: i,k
        
    ! Set again variables to zero
    mean_gw  = zero
    mean_gwx = zero
    mean_gwz = zero
    sh_vel   = zero
    sh_velx  = zero    
    sh_velz  = zero 
    
    ! Denominator of the divisions
    den = real(nx*nz,mytype)
    
    ! Transpose to y-pencils
    call transpose_x_to_y(ux,ux2)
    call transpose_x_to_y(uz,uz2)
 
    ! y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    
    ! du/dy=ta2   
    ! dw/dy=tc2
    
    ! Index for j is 1, since we are dealing with y-pencils and summation over all points (each processor)
    do k=1,ysize(3)
       do i=1,ysize(1)
              
           ! Total velocity gradient at the wall, sqrt[(du/dy)**2 + (dw/dy)**2] 
           mean_gw = mean_gw + sqrt_prec(ta2(i,1,k)**2 + tc2(i,1,k)**2) / den
           
           ! Mean streamwise gradient dU/dy
           mean_gwx = mean_gwx + ta2(i,1,k) / den
           
           ! Mean spanwise gradient dW/dy
           mean_gwz = mean_gwz + tc2(i,1,k) / den
                         
       enddo
    enddo
         
    ! Summation over all MPI processes and broadcast the result          
    call MPI_ALLREDUCE(mean_gw, sh_vel, 1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(mean_gwx,sh_velx,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(mean_gwz,sh_velz,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    ! Finalize shear velocities calculation
    sh_vel  = sqrt_prec(sh_vel  * xnu)
    sh_velx = sqrt_prec(abs_prec(sh_velx) * xnu)
    sh_velz = sqrt_prec(abs_prec(sh_velz) * xnu)  
                  
  end subroutine calculate_shear_velocity
  
  !---------------------------------------------------------------------------!
  ! Calculate the spanwise velocity at the wall due to the imposed
  ! sinusoidal oscillations.
  ! 
  ! Parameters of the non-dimensional oscillation (A^+ and T^+) are read from
  ! the input file.
  !---------------------------------------------------------------------------!
  subroutine spanwise_wall_oscillations(ux1,uz1)
  
  use dbg_schemes, only : sin_prec, sqrt_prec
  use param,       only : sh_vel, sh_velx, sh_velz, span_vel, t
  use param,       only : a_plus_cap, t_plus_cap
  use param,       only : two, xnu, pi
  use decomp_2d,   only : xsize
  use decomp_2d,   only : mytype
  
  implicit none
  
  real(mytype) :: amplitude, period
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uz1
  
  ! Calculate shear velocity    
  call calculate_shear_velocity(ux1,uz1,sh_vel,sh_velx,sh_velz)
    
  ! Maximum amplitude of spanwise oscillations
  !amplitude = sh_vel * a_plus_cap
  amplitude = sh_velx * a_plus_cap
  
  ! Period of oscillation
  !period = xnu * t_plus_cap / (sh_vel**2)
  period = xnu * t_plus_cap / (sh_velx**2)
  
  ! Calculation of the spanwise wall velocity
  span_vel = amplitude * sin_prec(two*pi*t/period)
   
  end subroutine spanwise_wall_oscillations
  
  !---------------------------------------------------------------------------!
  ! Calculate bulk velocity for a channel.
  ! Adapted from 'channel_cfr' subroutine.
  !---------------------------------------------------------------------------!
  subroutine calculate_ubulk(ux,uball)
  
  use decomp_2d, only : mytype, real_type, nrank
  use decomp_2d, only : xsize, zsize, xstart
  use param
  use MPI
  use variables, only : ppy
  
  implicit none

  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux

  integer                   :: code, i, j, k, jloc
  real(mytype)              :: ub, coeff
  real(mytype), intent(out) :: uball

  ub    = zero
  uball = zero
  coeff = dy / (yly * real(xsize(1) * zsize(3), kind=mytype))

  do k = 1, xsize(3)
     do jloc = 1, xsize(2)
        j = jloc + xstart(2) - 1
        do i = 1, xsize(1)
          ub = ub + ux(i,jloc,k) / ppy(j)
        enddo
     enddo
  enddo

  ub = ub * coeff
  
  call MPI_REDUCE(ub,uball,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    
  end subroutine calculate_ubulk 
  
  !---------------------------------------------------------------------------!
  ! Update coefficients for time integration schemes
  ! if adaptive time-step is used (max CFL condition)
  ! (taken from init_variables in variables module).
  !
  ! - Used in subroutine for CFL calculation in tools module.
  ! - At the moment, valid only for RK3 scheme.
  ! - To do: AB schemes.
  !---------------------------------------------------------------------------!
  
  subroutine update_time_int_coeff()
  
  use param
  
  implicit none
  
    ! Runge-Kutta 3 (RK3)
    adt(1)=(eight/fifteen)*dt
    bdt(1)= zero
    gdt(1)=adt(1)
    adt(2)=(      five/twelve)*dt
    bdt(2)=(-seventeen/ sixty)*dt
    gdt(2)=adt(2)+bdt(2)
    adt(3)=( three/four)*dt
    bdt(3)=(-five/twelve)*dt
    gdt(3)=adt(3)+bdt(3)
           
  end subroutine update_time_int_coeff
  !---------------------------------------------------------------------------! 
  
end module extra_tools




