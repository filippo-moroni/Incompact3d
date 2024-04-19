!----------------------------------------------------------!
!         This module is used to store useful              !
!      subroutines for temporal TBLs simulations.          !
!----------------------------------------------------------!
module tools_for_ttbl

  use decomp_2d
  use decomp_2d_io
  use variables
  use param 

  implicit none
  
  private 
  
  public :: calculate_friction_coefficient, update_time_int_coeff

contains

  !---------------------------------------------------------------------------!
  ! Calculate skin friction coefficient at the bottom wall and shear velocity
  ! Adapted from visu_ttbl subroutine.
  !
  ! - Used in BC-Temporal-TBL for the spanwise oscillations.
  ! - Used in tools module to calculate cf at the restart.
  !---------------------------------------------------------------------------!
  subroutine calculate_friction_coefficient(ux,uz)
    
    use var     
    use ibm_param,   only : ubcx,ubcy,ubcz
    use dbg_schemes, only : sqrt_prec
    
    use MPI
    
    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux, uz
    
    real(mytype) :: mean_gw      ! mean gradient at the wall at each processor
    real(mytype) :: mean_gw_tot  ! mean gradient at the wall total
    
    integer :: ierr  ! for MPI (initialized in init_xcompact3d subroutine)
    integer :: i,k
    
    ! Set again variables to zero
    mean_gw     = zero
    mean_gw_tot = zero
    fric_coeff  = zero
    sh_vel      = zero
    
    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux,ux2)
      call transpose_x_to_y(uz,uz2)
      sync_vel_needed = .false.
    endif

    ! y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    
    ! du/dy=ta2   
    ! dw/dy=tc2
    
    ! Mean velocity gradient at the wall, sqrt(du/dy**2 + dw/dy**2) and summation over all points
    do k=ystart(3),yend(3)
       do i=ystart(1),yend(1)
           
             ! Index for j is 1, since we are in global coordinates (y-pencils)
             mean_gw = mean_gw + (sqrt_prec(ta2(i,1,k)**2 + tc2(i,1,k)**2)) / real(nx*nz,mytype)                  
       enddo
    enddo
    
    ! Summation over all MPI processes and broadcast the result
    call MPI_ALLREDUCE(mean_gw,mean_gw_tot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    ! Calculate cf and shear velocity from the mean gradient at the wall    
    fric_coeff = mean_gw_tot * two * xnu / (uwall**2)
    sh_vel     = sqrt_prec(xnu * mean_gw_tot)
        
    return 
  
  end subroutine calculate_friction_coefficient
  
  !---------------------------------------------------------------------------!
  ! Update coefficients for time integration schemes
  ! if adaptive time-step is used (max CFL condition)
  ! (taken from init_variables in variables module).
  !
  ! - Used in subroutine for CFL calculation in tools module.
  !---------------------------------------------------------------------------!
  
  ! to be completed: time-steps are different among different steps!
  ! to be done: modify the AB schemes for variable time-step
  subroutine update_time_int_coeff()
  
  use param
  
  implicit none
  
    adt=zero
    bdt=zero
    cdt=zero
    gdt=zero

    if (itimescheme.eq.1) then ! Euler

       iadvance_time=1

       adt(1)=one*dt
       bdt(1)=zero
       gdt(1)=adt(1)+bdt(1)
       gdt(3)=gdt(1)

       ntime = 1
       nrhotime = 2
    elseif (itimescheme.eq.2) then ! AB2
       iadvance_time=1
       adt(1)=onepfive*dt
       bdt(1)=-half*dt
       gdt(1)=adt(1)+bdt(1)
       gdt(3)=gdt(1)

       ntime = 2
       nrhotime = 3
    elseif (itimescheme.eq.3) then ! AB3
       iadvance_time=1

       adt(1)= (twentythree/twelve)*dt
       bdt(1)=-(    sixteen/twelve)*dt
       cdt(1)= (       five/twelve)*dt
       gdt(1)=adt(1)+bdt(1)+cdt(1)
       gdt(3)=gdt(1)

       ntime = 3
       nrhotime = 4
    elseif(itimescheme==4) then  ! AB4
       iadvance_time=1

       adt(1)= (  fiftyfive/twentyfour)*dt
       bdt(1)=-(  fiftynine/twentyfour)*dt
       cdt(1)= (thirtyseven/twentyfour)*dt
       ddt(1)=-(       nine/twentyfour)*dt
       gdt(1)=adt(1)+bdt(1)+cdt(1)+ddt(1)
       gdt(3)=gdt(1)

       ntime    = 4
       nrhotime = 5
    elseif(itimescheme.eq.5) then !RK3
       iadvance_time=3

       adt(1)=(eight/fifteen)*dt
       bdt(1)= zero
       gdt(1)=adt(1)
       adt(2)=(      five/twelve)*dt
       bdt(2)=(-seventeen/ sixty)*dt
       gdt(2)=adt(2)+bdt(2)
       adt(3)=( three/four)*dt
       bdt(3)=(-five/twelve)*dt
       gdt(3)=adt(3)+bdt(3)

       ntime = 2
       nrhotime = 3
    elseif(itimescheme.eq.6) then !RK4 Carpenter and Kennedy
       iadvance_time=5
       adt(1)=zero
       adt(2)=-0.4178904745_mytype
       adt(3)=-1.192151694643_mytype
       adt(4)=-1.697784692471_mytype
       adt(5)=-1.514183444257_mytype
       bdt(1)=0.1496590219993_mytype
       bdt(2)=0.3792103129999_mytype
       bdt(3)=0.8229550293869_mytype
       bdt(4)=0.6994504559488_mytype
       bdt(5)=0.1530572479681_mytype
       gdt(1)=0.1496590219993_mytype*dt
       gdt(2)=0.220741935365_mytype*dt
       gdt(3)=0.25185480577_mytype*dt
       gdt(4)=0.33602636754_mytype*dt
       gdt(5)=0.041717869325_mytype*dt

       ntime = 2
       nrhotime = 5 ! (A guess)
       
    end if   
  
  return
  
  end subroutine update_time_int_coeff  
  
end module tools_for_ttbl




