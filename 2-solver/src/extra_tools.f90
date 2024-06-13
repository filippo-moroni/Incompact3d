!----------------------------------------------------------!
!         This module is used to store useful              !
! subroutines for general purpose, not present in standard !
!                  Incompact3d releases.                   !
!----------------------------------------------------------!

module extra_tools

  implicit none
  
  private 
  
  public :: print_cf,                   &
            calculate_shear_velocity,   &
            spanwise_wall_oscillations, &
            calculate_ubulk,            &
            calculate_bl_thick,         &
            write_scalar_plane_z 

contains
  
  !---------------------------------------------------------------------------!
  ! Write shear velocities, skin friction coefficients,
  ! viscous time unit, time unit, bulk velocity (channel only) 
  ! boundary layer thickness and Re_tau (TTBLs only) and stores
  ! them in a .txt file (used for TTBL and Channel).
  !---------------------------------------------------------------------------!
  subroutine print_cf(ux,uz)
  
  use param
  use decomp_2d
  use dbg_schemes, only : sqrt_prec
      
  implicit none
 
  ! Inputs 
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux, uz
  
  ! Locals
  integer :: iunit
  logical :: exists
  character(len=90) :: filename
  
  ! Write filename
  if (nrank .eq. 0) then
      write(filename,"('monitoring/cf_history.txt')") 
  end if
  
  ! TTBL
  if(itype .eq. itype_ttbl) then
      ! Shear velocity bottom wall
      call calculate_shear_velocity(ux,uz,sh_vel,sh_velx,sh_velz)
      
      ! Boundary layer thickness
      call calculate_bl_thick(ux,delta_99)
      
      ! Create or open a file to store sh_vel, cf coefficients, viscous time and time unit
      if(nrank .eq. 0) then
          ! Calculate friction coefficients
          fric_coeff  = two * ((sh_vel  / uwall)**2)
          fric_coeffx = two * ((sh_velx / uwall)**2)
          fric_coeffz = two * ((sh_velz / uwall)**2)
          
          ! Calculate friction Re number for a TBL
          re_tau_tbl = delta_99 * sh_velx / xnu
          
          ! Calculate viscous time unit
          if(iswitch_wo .eq. 1) then
              ! Moving walls
              t_viscous = xnu / (sh_velx**2)
          else
              ! Fixed walls
              t_viscous = xnu / (sh_vel**2)
          end if
          
          inquire(file=filename, exist=exists)
          if (exists) then
              open(newunit=iunit, file=filename, status="old", position="append", action="write")
              
              write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F16.10,A,F16.10,A,F16.10,A, F12.6,A,F12.4,A,I12,A, F12.6,A,F12.6)') &
                             sh_vel,     ',', sh_velx,     ',', sh_velz,     ',',                                         & 
                             fric_coeff, ',', fric_coeffx, ',', fric_coeffz, ',',                                         &
                             t_viscous,  ',', t,           ',', itime,       ',',                                         &
                             delta_99,   ',', re_tau_tbl
          else
              open(newunit=iunit, file=filename, status="new", action="write")
              ! Header
              write(iunit, '(A12,A,A12,A,A12,A, A16,A,A16,A,A16,A, A12,A,A12,A,A12,A, A12,A,A12)') &
                            'sh_vel',    ',', 'sh_velx', ',', 'sh_velz', ',',                      &
                            'cf,tot',    ',', 'cf,x',    ',', 'cf,z',    ',',                      &
                            't_nu',      ',', 'T',       ',', 'ts',      ',',                      &
                            'delta_99',  ',', 'Re_tau'          
              
              write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F16.10,A,F16.10,A,F16.10,A, F12.6,A,F12.4,A,I12,A, F12.6,A,F12.6)') &
                             sh_vel,     ',', sh_velx,     ',', sh_velz,     ',',                                         & 
                             fric_coeff, ',', fric_coeffx, ',', fric_coeffz, ',',                                         &
                             t_viscous,  ',', t,           ',', itime,       ',',                                         &
                             delta_99,   ',', re_tau_tbl
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
          
          inquire(file=filename, exist=exists)
          if (exists) then
              open(newunit=iunit, file=filename, status="old", position="append", action="write")
              
              write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F16.10,A,F16.10,A,F16.10,A, F12.6,A,F12.4,A,I12,A, F12.4)') &
                             sh_vel,     ',', sh_velx,     ',', sh_velz,     ',',                                 & 
                             fric_coeff, ',', fric_coeffx, ',', fric_coeffz, ',',                                 &
                             t_viscous,  ',', t,           ',', itime,       ',',                                 &
                             ubulk
          else
              open(newunit=iunit, file=filename, status="new", action="write")
              ! Header
              write(iunit, '(A12,A,A12,A,A12,A, A16,A,A16,A,A16,A, A12,A,A12,A,A12,A, A12)') &
                            'sh_vel', ',', 'sh_velx', ',', 'sh_velz', ',',                   &
                            'cf,tot', ',', 'cf,x',    ',', 'cf,z',    ',',                   &
                            't_nu',   ',', 'T',       ',', 'ts',      ',',                   &
                            'Ubulk'          
              
              write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F16.10,A,F16.10,A,F16.10,A, F12.6,A,F12.4,A,I12,A, F12.4)') &
                             sh_vel,     ',', sh_velx,     ',', sh_velz,     ',',                                 & 
                             fric_coeff, ',', fric_coeffx, ',', fric_coeffz, ',',                                 &
                             t_viscous,  ',', t,           ',', itime,       ',',                                 &
                             ubulk
          end if
              close(iunit)
      end if  
  
  end if ! closing of the if condition for different flow cases
             
  end subroutine print_cf
  
  !---------------------------------------------------------------------------!
  ! Calculate total shear velocity and its x and z components at bottom wall
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
    
    use param,       only : zero, xnu, itype, itype_channel, itype_ttbl
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
           
           ! TTBL, only bottom wall
           if (itype .eq. itype_ttbl) then
              
               ! Total velocity gradient at the wall, sqrt[(du/dy)**2 + (dw/dy)**2] 
               mean_gw = mean_gw + sqrt_prec(ta2(i,1,k)**2 + tc2(i,1,k)**2) / den
           
               ! Mean streamwise gradient dU/dy
               mean_gwx = mean_gwx + ta2(i,1,k) / den
           
               ! Mean spanwise gradient dW/dy
               mean_gwz = mean_gwz + tc2(i,1,k) / den
           
           ! Channel, upper wall summation too, with opposite sign
           else if (itype .eq. itype_channel) then
           
               ! Total velocity gradient at the wall, sqrt[(du/dy)**2 + (dw/dy)**2] 
               mean_gw = mean_gw + (sqrt_prec(ta2(i,1,k)**2 + tc2(i,1,k)**2) + sqrt_prec(ta2(i,ysize(2),k)**2 + tc2(i,ysize(2),k)**2)) / den
           
               ! Mean streamwise gradient dU/dy
               mean_gwx = mean_gwx + (ta2(i,1,k) - ta2(i,ysize(2),k)) / den
           
               ! Mean spanwise gradient dW/dy
               mean_gwz = mean_gwz + (tc2(i,1,k) - tc2(i,ysize(2),k)) / den
           
           end if              
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
  use param,       only : a_wo, t_wo, ifeedback_control
  use param,       only : two, xnu, pi
  use decomp_2d,   only : xsize
  use decomp_2d,   only : mytype
  
  implicit none
  
  real(mytype) :: amplitude, period
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uz1
  
  ! Rescale amplitude and period in friction units if feedback control is enabled (closed loop)
  if (ifeedback_control .eq. 1) then 
  
      ! Calculate shear velocity    
      call calculate_shear_velocity(ux1,uz1,sh_vel,sh_velx,sh_velz)
    
      ! Maximum amplitude of spanwise oscillations, based on longitudinal shear velocity
      amplitude = sh_velx * a_wo
  
      ! Period of oscillation, based on longitudinal shear velocity
      period = xnu * t_wo / (sh_velx**2)
  
  end if
  
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
  ! Calculate the boundary layer thickness for a TTBL.
  !---------------------------------------------------------------------------!
 
  subroutine calculate_bl_thick(ux,delta_99)
  
  use var,         only : ux2    
  use MPI
  use decomp_2d,   only : mytype, real_type, nrank
  use decomp_2d,   only : xsize, ysize, ystart, yend
  use decomp_2d,   only : transpose_x_to_y
    
  use param,       only : zpzeroone 
  use variables
    
  implicit none
  
  ! Inputs
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux
  
  ! Outputs
  real(mytype), intent(out) :: delta_99  ! BL thickness 
   
  ! Local variables 
  real(mytype), dimension(ysize(2)) :: u1meanH1,u1meanHT
  real(mytype)                      :: temp
  integer                           :: code,i,j,k
    
  ! Transpose data to y-pencils 
  call transpose_x_to_y(ux,ux2)
  
  ! Summation over x and z directions
  do k=ystart(3),yend(3)
      do i=ystart(1),yend(1)
          do j=ystart(2),yend(2)          
              u1meanH1(j)=u1meanH1(j)+ux2(i,j,k)/real(nx*nz,mytype)                                
          enddo          
      enddo
  enddo

  ! Summation over all MPI processes
  call MPI_REDUCE(u1meanH1,u1meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  
  ! Calculate temporal BL thickness (delta_99)
  if(nrank .eq. 0) then
      
      do j = ystart(2),yend(2)
     
          ! BL thickness equivalent to the y-coordinate
          delta_99 = yp(j)  
        
          ! %1 of the velocity at the wall (valid only for TTBLs with translating wall)
          temp = zpzeroone*u1meanHT(ystart(2))
                 
          if(u1meanHT(j) < temp) exit  
               
      end do
    
  end if
  
  end subroutine calculate_bl_thick
  
  !---------------------------------------------------------------------------! 
  ! Write a plane with z-dir. normal of the scalar field for visualization.
  ! 
  ! .xdmf header generation is present.
  !---------------------------------------------------------------------------!
   
  subroutine write_scalar_plane_z(phi1,itime)
 
  use visu
  
  use decomp_2d,    only : mytype, xsize, nrank
  use decomp_2d_io, only : decomp_2d_start_io
  use param,        only : ioutput_plane
  use variables,    only : numscalar
  
  implicit none
  
  ! Inputs
  real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
  integer,                                                          intent(in) :: itime
  
  ! Locals
  character(len=32) :: num  ! taken from write_snapshot in visu module
      
  ! Switch to output2D with z-normal plane
  output2D = 3

!--- Write snapshot part ---!
  
#ifdef ADIOS2
  call decomp_2d_start_io(io_name, "data")
#endif
    
  ! Snapshot number
#ifndef ADIOS2
  if (filenamedigits) then
     ! New enumeration system, it works integrated with xcompact3d_toolbox
     write(num, ifilenameformat) itime
  else
     ! Classic enumeration system
     write(num, ifilenameformat) itime/ioutput_plane
  endif
#else
  ! ADIOS2 is zero-indexed
  write(num, '(I0)') itime/ioutput_plane - 1
#endif
    
  ! Write XDMF header
  call write_xdmf_header("planes", "phiplane", trim(num))
   
  ! Write first scalar field
  call write_field(phi1(:,:,:,1), "planes", "phiplane", trim(num))
  
!--- End snapshot part ---!
  
  ! Write XDMF footer
  call write_xdmf_footer()
  
#ifdef ADIOS2
  call decomp_2d_end_io(io_name, "data")
#endif

!-------------------------!
  
  ! Switch back to 3D output for default snapshots
  output2D = 0
  
  ! call the creation of the .xdmf file
  
  ! Write the plane
  !call decomp_2d_write_plane(1,phi1(:,:,:,1),3,1,"data",gen_filename("planes", "phi01-p", num, 'bin'),io_name)
  
  end subroutine write_scalar_plane_z
 
end module extra_tools




