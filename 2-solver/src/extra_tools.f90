!----------------------------------------------------------!
!         This module is used to store useful              !
!           subroutines for general purpose.               !
!----------------------------------------------------------!

module extra_tools

  implicit none
  
  private 
  
  public :: print_cf,                       &
            calculate_mgw,                  &
            spanwise_wall_oscillations,     &
            update_time_int_coeff

contains
  
  !---------------------------------------------------------------------------!
  ! Write skin friction coefficient, shear velocity and time unit and stores
  ! them in a .txt file for time-evolution (used for TTBL and Channel).
  !---------------------------------------------------------------------------!
  subroutine print_cf(ux,uz)
  
  use param
  use decomp_2d
  use dbg_schemes, only : sqrt_prec
    
  implicit none
  
  integer :: iunit
  logical :: exists
  
  ! Inputs from 'calculate_mgw' subroutine
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux, uz
  
  ! 
  mean_gw_tot = zero
   
  ! Calculate skin friction coefficient and shear velocity at the bottom wall (TTBL and Channel)
  if(itype .eq. itype_ttbl .or. itype .eq. itype_channel) then
      call calculate_mgw(ux,uz)  
  end if   
     
  ! Create or open a file to store sh_vel, cf and time unit
  if(nrank .eq. 0) then
      
      ! Calculate shear velocity and friction coefficient
      sh_vel     = sqrt_prec(xnu * mean_gw_tot)
      fric_coeff = two * ((sh_vel / uwall)**2)
      
      inquire(file="cf_history.txt", exist=exists)
      if (exists) then
          open(newunit=iunit, file="cf_history.txt", status="old", position="append", action="write")
          write(iunit, '(F8.6,A,F8.6,A,F8.6)') sh_vel, ',', fric_coeff, ',', t
      else
          open(newunit=iunit, file="cf_history.txt", status="new", action="write")
          write(iunit, '(A8,A,A8,A,A8)')    'sh_vel', ',', 'cf', ',', 'T'          
          write(iunit, '(F8.6,A,F8.6,A,F8.6)') sh_vel, ',', fric_coeff, ',', t
      end if
      close(iunit)
  end if
           
  end subroutine print_cf
  
  !---------------------------------------------------------------------------!
  ! Calculate mean gradient at the wall
  !
  ! - Used in BC-Temporal-TBL and in BC-Channel-flow
  !   for the spanwise wall oscillations. 
  !---------------------------------------------------------------------------!
  subroutine calculate_mgw(ux,uz)
    
    use var,         only : ux2, uz2     
    use ibm_param,   only : ubcx,ubcy,ubcz
    use dbg_schemes, only : sqrt_prec
    
    use var,         only : ta2,tc2,di2
    use ibm_param,   only : ubcx,ubcz
    
    use MPI
    use decomp_2d
    
    use param,       only : zero, mean_gw_tot
    use variables
    
    implicit none
    
    ! Inputs
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux, uz

    ! Outputs
   
    !real(mytype), intent(out) :: mean_gw_tot  ! mean gradient at the wall total
    
    ! Work variables
    real(mytype) :: mean_gw      ! mean gradient at the wall at each processor
    integer      :: ierr         
    integer      :: i,k
        
    ! Set again variables to zero
    mean_gw     = zero
    mean_gw_tot = zero
    !mean_gw_tot = zero
    
    ! Transpose to y-pencils
    call transpose_x_to_y(ux,ux2)
    call transpose_x_to_y(uz,uz2)
 
    ! y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    
    ! du/dy=ta2   
    ! dw/dy=tc2
        
    ! Mean velocity gradient at the wall, sqrt[(du/dy)**2 + (dw/dy)**2] and summation over all points
    do k=ystart(3),yend(3)
       do i=ystart(1),yend(1)
           
           ! Index for j is 1, since we are in global coordinates (y-pencils)
           mean_gw = mean_gw + sqrt_prec(ta2(i,1,k)**2 + tc2(i,1,k)**2) / real(nx*nz,mytype)                   
       enddo
    enddo
     
    ! Summation over all MPI processes and broadcast the result
    !call MPI_ALLREDUCE(mean_gw,mean_gw_tot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    write(*,*) mean_gw_tot
    
    call MPI_REDUCE(mean_gw,mean_gw_tot,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
    write(*,*) mean_gw_tot
              
  end subroutine calculate_mgw
  
  !---------------------------------------------------------------------------!
  ! Calculate the spanwise velocity at the wall due to the imposed
  ! sinusoidal oscillation.
  ! 
  ! Parameters of the non-dimensional oscillation (A^+ and T^+) are read from
  ! the input file.
  !---------------------------------------------------------------------------!
  subroutine spanwise_wall_oscillations(ux1,uz1)
  
  use dbg_schemes, only : sin_prec, sqrt_prec
  use param,       only : sh_vel, span_vel, fric_coeff, t
  use param,       only : a_plus_cap, t_plus_cap
  use param,       only : two, xnu, pi, uwall, mean_gw_tot
  use decomp_2d,   only : xsize, xstart, xend 
  use decomp_2d,   only : mytype
  
  implicit none
  
  real(mytype) :: amplitude, period
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uz1
    
  call calculate_mgw(ux1,uz1)
  
  ! Calculate shear velocity and friction coefficient
  sh_vel     = sqrt_prec(xnu * mean_gw_tot)
  fric_coeff = two * ((sh_vel / uwall)**2)
  
  ! Maximum amplitude of spanwise oscillations
  amplitude = sh_vel * a_plus_cap
  
  ! Period of oscillation
  period = xnu * t_plus_cap / (sh_vel**2)
  
  ! Calculation of the spanwise wall velocity
  span_vel = amplitude * sin_prec(two*pi*t/period)
   
  end subroutine spanwise_wall_oscillations
  
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
  
end module extra_tools




