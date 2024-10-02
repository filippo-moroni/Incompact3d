
!This file is not part of standard Xcompact3d releases (xcompact3d.com).

!-----------------------------------------------------------------------------!
! DESCRIPTION: This file is used to store useful subroutines for general 
!              purpose, not present in standard Incompact3d releases:
!              - 'print_cf' 
!              - 'calculate_shear_velocity' 
!              - 'calculate_scalar_grad_wall'
!              - 'spanwise_wall_oscillations' 
!              - 'calculate_ubulk' 
!              - 'calculate_bl_thick'
!              - 'print_mean_stats'
!              - 'print_mean_stats_scalar'
!              - 'write_scalar_plane_z'
!              - 'write_vortx_plane_x'.    
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
! DESCRIPTION: Write shear velocities, skin friction coefficient, viscous 
!              time unit, time unit, bulk velocity (channel only), boundary 
!              layer thickness (delta_99) and Re_tau (TTBL only) and stores
!              them in a .txt file (used for TTBLs and Channels). 
!              We are also calling subroutines to print runtime mean 
!              statistics (TTBL only).    
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine print_cf(ux,uy,uz,phi)
  
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use param
  use variables, only : nx,ny,nz,yp,numscalar,sc
      
  implicit none
 
  ! Inputs 
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)),           intent(in) :: ux,uy,uz
  real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi
  
  ! Locals
  integer :: iunit, j
  integer :: nyh         ! half - 1 points in y direction for a channel (h: half)
  logical :: exists
  character(len=90) :: filename
  
  ! Local value of amplitude and period of wall oscillations in wall units
  real(mytype) :: a_wo_loc,t_wo_loc

  ! Total shear velocity
  real(mytype) :: sh_vel
      
  ! Write filename
  if (nrank .eq. 0) write(filename,"('data/monitoring/cf_history.txt')")
  
  ! Calculate shear velocities
  call calculate_shear_velocity(ux,uz,sh_velx,sh_velz)

  ! Calculate total shear velocity (based on the total wall shear stress)
  sh_vel = (sh_velx**4 + sh_velz**4)**0.25

  ! Calculate viscous time unit (based on the total shear velocity)            
  t_viscous = xnu / (sh_vel**2)

  ! Calculate mean scalar gradient
  if(iscalar .eq. 1) call calculate_scalar_grad_wall(phi,mean_phigwtot)

  ! TTBL
  if(itype .eq. itype_ttbl) then
  
      ! Write mean statistics runtime in case of a TTBL
      call print_mean_stats(ux,uy,uz)
      
      ! Write mean statistics runtime for scalar field in case of a TTBL 
      if(iscalar .eq. 1) call print_mean_stats_scalar(ux,uy,uz,phi) 
             
      ! Boundary layer thickness
      call calculate_bl_thick(ux,delta_99,counter)
      
      ! Create or open a file to store shear velocities, friction coefficient, viscous time, time unit
      ! TTBL thickness delta_99, friction Reynolds number and power input 
      if(nrank .eq. 0) then

          ! Calculate (streamwise) friction coefficient
          fric_coeff  = two * ((sh_velx  / uwall)**2)

          ! Calculate friction Re number for a TTBL
          re_tau_tbl = delta_99 * sh_velx / xnu

          ! Calculate power input (assuming unitary density)
          if(iswitch_wo .eq. 0) then
              
              ! Fixed walls
              powerin = (sh_velx**2)*uwall
          
          else if(iswitch_wo .eq. 1) then
          
              ! Oscillating walls 
              powerin = (sh_velx**2)*uwall + (sh_velz**2)*span_vel
          
          end if
        
          ! Dissimilar control section
          if (iscalar .eq. 1) then
              
              ! Analogy factor
              A_fact = (xnu / sc(1)) * mean_phigwtot / sh_vel**2 
              
          end if
          
          inquire(file=filename, exist=exists)
          if (exists) then
              open(newunit=iunit, file=filename, status="old", position="append", action="write")
              
              write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F16.10,A,F12.6,A,F12.6,A, F12.6,A,F12.4,A,I12,A, F12.6,A,F12.6,A,F12.6)') &
                             sh_vel,     ',', sh_velx,       ',', sh_velz, ',',                                                 & 
                             fric_coeff, ',', mean_phigwtot, ',', A_fact,  ',',                                                 &
                             t_viscous,  ',', t,             ',', itime,   ',',                                                 &
                             delta_99,   ',', re_tau_tbl,    ',', powerin                                         
                             
          else
              open(newunit=iunit, file=filename, status="new", action="write")
              ! Header
              write(iunit, '(A12,A,A12,A,A12,A, A16,A,A12,A,A12,A, A12,A,A12,A,A12,A, A12,A,A12,A,A12)') &
                            'sh_vel',    ',', 'sh_velx',    ',', 'sh_velz', ',',                         &
                            'cfx',       ',', '(dPhi/dy)w', ',', 'A_fact',  ',',                         &
                            't_nu',      ',', 't',          ',', 'ts',      ',',                         &
                            'delta_99',  ',', 'Re_tau',     ',', 'P_in'                                      
                            
              
              write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F16.10,A,F12.6,A,F12.6,A, F12.6,A,F12.4,A,I12,A, F12.6,A,F12.6,A,F12.6)') &
                             sh_vel,     ',', sh_velx,       ',', sh_velz, ',',                                                 & 
                             fric_coeff, ',', mean_phigwtot, ',', A_fact,  ',',                                                 &
                             t_viscous,  ',', t,             ',', itime,   ',',                                                 &
                             delta_99,   ',', re_tau_tbl,    ',', powerin 
          end if
              close(iunit)
      
      end if
 
  ! Channel
  else if(itype .eq. itype_channel) then
            
      ! Bulk velocity
      call calculate_ubulk(ux,ubulk)
    
      ! Create or open a file to store shear velocities, friction coefficient, bulk velocity, viscous time and time unit 
      if(nrank .eq. 0) then

          ! Calculate (streamwise) friction coefficient
          fric_coeff  = two * ((sh_velx  / ubulk)**2)
                   
          inquire(file=filename, exist=exists)
          if (exists) then
              open(newunit=iunit, file=filename, status="old", position="append", action="write")
              
              write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F16.10,A,F12.4,A, F12.6,A,F12.4,A,I12,A)') &
                             sh_vel,     ',', sh_velx,     ',', sh_velz,     ',',                & 
                             fric_coeff, ',', ubulk,       ',',                                  &
                             t_viscous,  ',', t,           ',', itime
                             
          else
              open(newunit=iunit, file=filename, status="new", action="write")
              ! Header
              write(iunit, '(A12,A,A12,A,A12,A, A16,A,A12,A, A12,A,A12,A,A12,A)') &
                            'sh_vel',    ',', 'sh_velx',   ',', 'sh_velz', ',',   &
                            'cfx',       ',', 'Ubulk',     ',',                   &
                            't_nu',      ',', 't',         ',', 'ts'
                                      
              
              write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F16.10,A,F12.4,A, F12.6,A,F12.4,A,I12,A)') &
                             sh_vel,     ',', sh_velx,     ',', sh_velz,   ',',                  & 
                             fric_coeff, ',', ubulk,       ',',                                  &
                             t_viscous,  ',', t,           ',', itime
          end if
              close(iunit)
      end if  
  
  end if ! closing of the if condition for different flow cases
  
  ! Create or open a file to store non-dimensional mesh spacings and domain dimensions
  if(nrank .eq. 0) then
  
      ! Open and read yp coordinates
      open(newunit=iunit, file="yp.dat", status="old", form="formatted")
      do j=1,ny
          read(iunit,*) yp(j)
      enddo
      close(iunit)
  
      ! Write filename
      write(filename,"('data/monitoring/grid_spacings.txt')") 
   
      ! Calculate viscous unit with total shear velocity
      delta_nu = xnu / sh_vel 
   
      ! Calculate non-dimensional grid spacings
      deltaxplus  = xlx / nx / delta_nu
      deltayplusw = yp(2) / delta_nu
      deltazplus  = zlz / nz / delta_nu
      
      ! Calculate non-dimensional domain dimensions
      xlxplus = xlx / delta_nu
      ylyplus = yly / delta_nu
      zlzplus = zlz / delta_nu
      
      ! Calculate grid spacing at interface (TTBL) or at centerline (Channel)
      if(itype .eq. itype_ttbl) then
       
          deltayplusd = (yp(counter + 1) - yp(counter)) / delta_nu
      
      else if(itype .eq. itype_channel) then 
          
          ! Half - 1 points in y-dir. for a channel
          nyh = (ny - 1)/2
          
          deltayplusd = (yp(nyh) - yp(nyh - 1)) / delta_nu
          
      end if
      
      ! Half the height for a channel
      if(itype .eq. itype_channel) ylyplus = ylyplus / two    
          
      inquire(file=filename, exist=exists)
      if (exists) then
          open(newunit=iunit, file=filename, status="old", position="append", action="write")
              
          write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F12.6,A,F12.6,A,F12.6,A, F12.6,A, F12.4,A,I12)') &
                         deltaxplus,  ',', deltayplusw,  ',', deltazplus,  ',',                    & 
                         xlxplus,     ',', ylyplus,      ',', zlzplus,     ',',                    &
                         deltayplusd, ',',                                                         &
                         t,           ',', itime
      else
          open(newunit=iunit, file=filename, status="new", action="write")
          ! Header
          write(iunit, '(A12,A,A12,A,A12,A, A12,A,A12,A,A12,A, A12,A, A12,A,A12)') &
                        'delta_x^+',  ',', 'delta_yw^+', ',', 'delta_z^+', ',',    &
                        'Lx^+',       ',', 'Ly^+',       ',', 'Lz^+',      ',',    &
                        'delta_yd^+', ',',                                         &
                        't',          ',', 'ts'          
              
          write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F12.6,A,F12.6,A,F12.6,A, F12.6,A, F12.4,A,I12)') &
                         deltaxplus,  ',', deltayplusw,  ',', deltazplus,  ',',                    & 
                         xlxplus,     ',', ylyplus,      ',', zlzplus,     ',',                    &
                         deltayplusd, ',',                                                         &
                         t,           ',', itime
      end if
          close(iunit)
  end if
  
  ! Create or open a file to store wall oscillation parameters in case of a TTBL without feedback control
  if(nrank .eq. 0 .and. itype .eq. itype_ttbl .and. iswitch_wo .eq. 1 .and. ifeedback_control .eq. 0) then
   
      ! Write filename
      write(filename,"('data/monitoring/oscill_param.txt')") 
   
      ! Calculate amplitude and period of oscillations in wall units (loc: local value)
      a_wo_loc = a_wo / sh_velx
      t_wo_loc = t_wo / t_viscous
      
      ! Write to file   
      inquire(file=filename, exist=exists)
      if (exists) then
          open(newunit=iunit, file=filename, status="old", position="append", action="write")
              
          write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F12.6,A,F12.4,A,I12)') &
                         a_wo_loc, ',', t_wo_loc, ',', re_tau_tbl, ',',  &  
                         span_vel, ',', t,        ',', itime
                             
      else
          open(newunit=iunit, file=filename, status="new", action="write")
          ! Header
          write(iunit, '(A12,A,A12,A,A12,A, A12,A,A12,A,A12)') &
                        'A^+', ',', 'T^+', ',', 'Re_tau', ',', &                                                        
                        'Ww',  ',', 't',   ',', 'ts'
          
          write(iunit, '(F12.6,A,F12.6,A,F12.6,A, F12.6,A,F12.4,A,I12)') &
                         a_wo_loc, ',', t_wo_loc, ',', re_tau_tbl, ',',  & 
                         span_vel, ',', t,        ',', itime                              
      end if
          close(iunit)       
  end if
                 
end subroutine print_cf

!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculate streamwise and spanwise shear velocities
!              (if squared, they correspond to the mean wall-shear stress
!              along a specific direction). 
!              - Used in 'BC-Temporal-TBL.f90' and in 'BC-Channel-flow.f90'
!                for the spanwise wall oscillations with feedback control 
!                enabled. 
!              - Used to print cf coefficients and shear velocities to an 
!                overall .txt file for time evolution check.    
! ANNOTATIONS: These components cannot be used to calculate
!              the total shear velocity as a vector magnitude. 
!              We must obtain first wall-shear stress components (x and z).
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine calculate_shear_velocity(ux,uz,sh_velx,sh_velz)

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use MPI

  use var,       only : ux2,uz2,ta2,tc2,di2
  use ibm_param, only : ubcx,ubcz
  use param,     only : zero, two, xnu, itype, itype_channel, itype_ttbl
  use variables

  implicit none

  ! Inputs
  real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux, uz

  ! Outputs
  real(mytype), intent(out) :: sh_velx ! Shear velocity along x 
  real(mytype), intent(out) :: sh_velz ! Shear velocity along z

  ! Work variables
  real(mytype) :: mean_gwx   ! Mean gradient direction x at each processor
  real(mytype) :: mean_gwz   ! Mean gradient direction z at each processor
  real(mytype) :: den        ! Denominator of the divisions

  integer      :: ierr
  integer      :: i,k

  ! Set variables to zero
  mean_gwx = zero
  mean_gwz = zero
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

         ! TTBL, only bottom wall (j = 1)
         if (itype .eq. itype_ttbl) then

             ! Mean streamwise gradient dU/dy
             mean_gwx = mean_gwx + ta2(i,1,k) / den

             ! Mean spanwise gradient dW/dy
             mean_gwz = mean_gwz + tc2(i,1,k) / den

         ! Channel, upper wall summation too, with opposite sign (j = ysize(2))
         else if (itype .eq. itype_channel) then

             ! Mean streamwise gradient dU/dy
             mean_gwx = mean_gwx + (ta2(i,1,k) - ta2(i,ysize(2),k)) / den / two

             ! Mean spanwise gradient dW/dy
             mean_gwz = mean_gwz + (tc2(i,1,k) - tc2(i,ysize(2),k)) / den / two

         end if
     enddo
  enddo

  ! Summation over all MPI processes and broadcast the result          
  call MPI_ALLREDUCE(mean_gwx,sh_velx,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(mean_gwz,sh_velz,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)

  ! Finalize shear velocities calculation
  sh_velx = sqrt(abs(sh_velx) * xnu)
  sh_velz = sqrt(abs(sh_velz) * xnu)

end subroutine calculate_shear_velocity

!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculate mean scalar gradient at the wall.    
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!  
subroutine calculate_scalar_grad_wall(phi,mean_phigwtot)
         
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use MPI 
              
  use var,       only : td2,di2,phi2
  use param,     only : zero,two,itype,itype_ttbl,itype_channel
  use variables
    
  implicit none
    
  ! Inputs
  real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi

  ! Outputs
  real(mytype), intent(out) :: mean_phigwtot  ! Mean scalar gradient at the wall (all processors)
      
  ! Work variables
  real(mytype) :: mean_phigw  ! Mean scalar gradient at each processor
  real(mytype) :: den         ! Denominator of the divisions 
  integer      :: ierr         
  integer      :: i,k
                
  ! Set variables to zero
  mean_phigw = zero  
  mean_phigwtot = zero  
  
  ! Denominator of the divisions
  den = real(nx*nz,mytype)
      
  ! Transpose to y-pencils
  call transpose_x_to_y(phi(:,:,:,1), phi2(:,:,:,1))
     
  ! y-derivative, scalar field (IBM parameter is zero, similarly to what is done for 2nd derivative in navier module for mass fraction)
  call dery (td2,phi2(:,:,:,1),di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,zero)
     
  ! dphi/dy = td2
    
  ! Index for j is 1, since we are dealing with y-pencils and summation over all points (each processor)
  do k=1,ysize(3)
     do i=1,ysize(1)
         
        ! TTBL, only bottom wall (j = 1)
        if (itype .eq. itype_ttbl) then
                                 
            ! Mean scalar gradient dPhi/dy at the wall
            mean_phigw = mean_phigw + td2(i,1,k) / den
         
        ! Channel, upper wall summation too, with opposite sign (j = ysize(2))
        else if (itype .eq. itype_channel) then

            ! Mean scalar gradient dPhi/dy at the wall
            mean_phigw = mean_phigw + (td2(i,1,k) - td2(i,ysize(2),k)) / den / two
        
        end if     
     enddo
  enddo
         
  ! Summation over all MPI processes and broadcast the result          
  call MPI_ALLREDUCE(mean_phigw,mean_phigwtot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
  
  ! Finalize calculation of mean scalar gradient at the wall
  mean_phigwtot = abs(mean_phigwtot)
      
end subroutine calculate_scalar_grad_wall

!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculate the spanwise velocity at the wall due to the 
!              imposed sinusoidal oscillations. Oscillation parameters 
!              (A and T) are read from the input file. 
!              With no feedback control, A and T are in outer units.
!              With feedback control,    A and T are in viscous units, 
!              (A^+, T^+).    
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine spanwise_wall_oscillations(ux1,uz1)

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d

  use param,       only : sh_velx, sh_velz, span_vel, t
  use param,       only : a_wo, t_wo, ifeedback_control, in_phase
  use param,       only : two, xnu, pi

  implicit none

  real(mytype) :: amplitude, period
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uz1

  ! Amplitude and period in outer units if feedback control is disabled (open loop)
  if (ifeedback_control .eq. 0) then

      amplitude = a_wo

      period = t_wo

  ! Rescale amplitude and period in friction units if feedback control is enabled (closed loop)
  else if (ifeedback_control .eq. 1) then

      ! Calculate shear velocity    
      call calculate_shear_velocity(ux1,uz1,sh_velx,sh_velz)

      ! Maximum amplitude of spanwise oscillations, based on longitudinal shear velocity
      amplitude = sh_velx * a_wo

      ! Period of oscillation, based on longitudinal shear velocity
      period = xnu * t_wo / (sh_velx**2)

  end if

  ! Calculation of the spanwise wall velocity
  span_vel = amplitude * sin(two*pi*t/period + in_phase*pi)

end subroutine spanwise_wall_oscillations

!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculate bulk velocity for a channel.
!              Adapted from 'channel_cfr' subroutine.    
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------! 
subroutine calculate_ubulk(ux,uball)
  
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
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

!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculate the boundary layer thickness (delta_99) for a TTBL.    
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!       
subroutine calculate_bl_thick(ux,delta_99,counter)
  
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use var,         only : ux2   
  use MPI
    
  use param,       only : zpzeroone, zero
  use variables,   only : yp, nx, ny, nz
    
  implicit none
  
  ! Inputs
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux
  
  ! Outputs
  real(mytype), intent(out) :: delta_99  ! BL thickness
  integer,      intent(out) :: counter   ! counter for the index of the BL mean interface for yp coordinates 
   
  ! Local variables 
  real(mytype), dimension(ysize(2)) :: u1meanH1,u1meanHT
  real(mytype)                      :: temp, den
  integer                           :: code,i,j,k,iunit
  
  ! Set variables to zero
  u1meanH1 = zero
  u1meanHT = zero
    
  ! Denominator of the divisions
  den = real(nx*nz,mytype)
      
  ! Transpose data to y-pencils 
  call transpose_x_to_y(ux,ux2)
  
  ! Summation over x and z directions
  do k=1,ysize(3)
      do i=1,ysize(1)
          do j=1,ysize(2)          
              u1meanH1(j)=u1meanH1(j)+ux2(i,j,k)/den                               
          enddo          
      enddo
  enddo

  ! Summation over all MPI processes
  call MPI_REDUCE(u1meanH1,u1meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  
  ! Calculate temporal BL thickness (delta_99)
  if(nrank .eq. 0) then
      
      ! Set again delta_99 and counter to 0 
      delta_99 = zero
      counter  = 0
      
      ! Open and read yp coordinates
      open(newunit=iunit, file="yp.dat", status="old", form="formatted")
      do j=1,ny
          read(iunit,*) yp(j)
      enddo
      close(iunit)
      
      ! Condition: %1 of the velocity at the wall (valid only for TTBLs with translating wall)
      temp = zpzeroone*u1meanHT(1)
      
      ! Initialize the index
      j = 1
      
      ! Cycle to check the BL thickness
      do while(u1meanHT(j) > temp)
      
          delta_99 = yp(j)
      
          j = j + 1
      
          counter = j
      
      end do
          
  end if
  
end subroutine calculate_bl_thick

!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculate runtime mean statistics (umean, wmean, variances and
!              Reynolds stress) for a TTBL.
!              This subroutine is called inside 'print_cf'. 
!              The streamwise mean velocity profile is used to calculate 
!              high order integrals of BL thickness parameters (delta*, theta) 
!              during post_processing.
!              These statistics must be finalized with different flow 
!              realizations. var[u] must be finalized further, as:
!                         var[u] = var[u] - mean[u]**2 
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!       
subroutine print_mean_stats(ux,uy,uz)
  
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use MPI
  
  use var,       only : ux2,uy2,uz2     
  use param,     only : zero,itime
  use variables, only : nx,ny,nz
    
  implicit none
  
  ! Inputs
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
     
  ! Local variables
  real(mytype)                      :: den
  integer                           :: code,i,j,k,iunit
  character(99)                     :: filename 
  
  ! Mean velocity profile
  real(mytype), dimension(ysize(2)) :: u1meanH1,u1meanHT  ! mean streamwise  velocity
  real(mytype), dimension(ysize(2)) :: v1meanH1,v1meanHT  ! mean wall-normal velocity
  real(mytype), dimension(ysize(2)) :: w1meanH1,w1meanHT  ! mean spanwise    velocity
  
  ! Variances
  real(mytype), dimension(ysize(2)) :: u2meanH1,u2meanHT  ! streamwise  velocity variance
  real(mytype), dimension(ysize(2)) :: v2meanH1,v2meanHT  ! wall-normal velocity variance
  real(mytype), dimension(ysize(2)) :: w2meanH1,w2meanHT  ! spanwise    velocity variance

  ! Reynolds stress
  real(mytype), dimension(ysize(2)) :: uvmeanH1,uvmeanHT  ! Reynolds stress <u'v'>
  real(mytype), dimension(ysize(2)) :: uwmeanH1,uwmeanHT  ! Reynolds stress <u'w'>
  real(mytype), dimension(ysize(2)) :: vwmeanH1,vwmeanHT  ! Reynolds stress <v'w'>
  
  ! Set variables to zero
  u1meanH1 = zero; u1meanHT = zero
  v1meanH1 = zero; v1meanHT = zero  
  w1meanH1 = zero; w1meanHT = zero

  u2meanH1 = zero; u2meanHT = zero
  v2meanH1 = zero; v2meanHT = zero  
  w2meanH1 = zero; w2meanHT = zero
  
  uvmeanH1 = zero; uvmeanHT = zero
  uwmeanH1 = zero; uwmeanHT = zero  
  vwmeanH1 = zero; vwmeanHT = zero

  ! Denominator of the divisions
  den = real(nx*nz,mytype)
      
  ! Transpose data to y-pencils 
  call transpose_x_to_y(ux,ux2)
  call transpose_x_to_y(uy,uy2)
  call transpose_x_to_y(uz,uz2)
    
  ! Summation over x and z directions
  do k=1,ysize(3)
      do i=1,ysize(1)
          do j=1,ysize(2)
          
              ! Mean velocity field          
              u1meanH1(j)=u1meanH1(j)+ux2(i,j,k)/den
              v1meanH1(j)=v1meanH1(j)+uy2(i,j,k)/den
              w1meanH1(j)=w1meanH1(j)+uz2(i,j,k)/den
              
              ! Variances
              u2meanH1(j)=u2meanH1(j)+ux2(i,j,k)*ux2(i,j,k)/den
              v2meanH1(j)=v2meanH1(j)+uy2(i,j,k)*uy2(i,j,k)/den
              w2meanH1(j)=w2meanH1(j)+uz2(i,j,k)*uz2(i,j,k)/den                                                              
          
              ! Reynolds stresses
              uvmeanH1(j)=uvmeanH1(j)+ux2(i,j,k)*uy2(i,j,k)/den
              uwmeanH1(j)=uwmeanH1(j)+ux2(i,j,k)*uz2(i,j,k)/den
              vwmeanH1(j)=vwmeanH1(j)+uy2(i,j,k)*uz2(i,j,k)/den
          
          enddo          
      enddo
  enddo

  ! Summation over all MPI processes
  call MPI_REDUCE(u1meanH1,u1meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(v1meanH1,v1meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(w1meanH1,w1meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

  call MPI_REDUCE(u2meanH1,u2meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(v2meanH1,v2meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(w2meanH1,w2meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

  call MPI_REDUCE(uvmeanH1,uvmeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(uwmeanH1,uwmeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(vwmeanH1,vwmeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  
  ! Print mean statistics for velocity field 
  if(nrank .eq. 0) then

      ! Write filename 
      write(filename, "('data/mean_stats_runtime/velocity/mean_stats_runtime-ts',I7.7,'.txt')") itime
            
      ! Open and write the mean statistics for velocity field
      open(newunit=iunit, file=filename, status="unknown", form="formatted")

      ! Header
      write(iunit, *) 'Mean statistics calculated runtime.'
      write(iunit, *) ' '
      write(iunit, *) 'Pay attention that these statistics need to be averaged later'
      write(iunit, *) 'with different flow realizations.'
      write(iunit, *) ' '
      write(iunit, *) 'Only var[u] needs to be finalized later by subtracting (mean[u])**2,'
      write(iunit, *) 'due to the symmetries of the TTBL.'
      write(iunit, *) ' '     
      write(iunit, *) 'The finalization of all 2nd order statistics can be performed anyway in order to check'
      write(iunit, *) 'the correctness of the calculations.' 
      write(iunit, *) ' '
      write(iunit, '(9(A14, A1, 1X))') 'umean(y,t,nr)' , ',', &
                                       'vmean(y,t,nr)' , ',', &
                                       'wmean(y,t,nr)' , ',', &
                                       'var[u](y,t,nr)', ',', &
                                       'var[v](y,t,nr)', ',', &
                                       'var[w](y,t,nr)', ',', &
                                       "<u'v'>(y,t,nr)", ',', &
                                       "<u'w'>(y,t,nr)", ',', &
                                       "<v'w'>(y,t,nr)"

      ! Write mean statistics, function of y-direction, time and specific realization
      do j=1,ny
          write(iunit, '(9(F14.9, A1, 1X))') u1meanHT(j), ',', &  ! Mean streamwise  velocity
                                             v1meanHT(j), ',', &  ! Mean wall-normal velocity
                                             w1meanHT(j), ',', &  ! Mean spanwise    velocity
                                             u2meanHT(j), ',', &  ! Streamwise  velocity variance
                                             v2meanHT(j), ',', &  ! Wall-normal velocity variance
                                             w2meanHT(j), ',', &  ! Spanwise    velocity variance
                                             uvmeanHT(j), ',', &  ! Reynolds stress <u'v'>
                                             uwmeanHT(j), ',', &  ! Reynolds stress <u'w'>
                                             vwmeanHT(j)          ! Reynolds stress <v'w'>
      enddo
      
      close(iunit)
                      
  end if
  
end subroutine print_mean_stats

!-----------------------------------------------------------------------------!
! DESCRIPTION: Calculate runtime mean statistics for the scalar field
!              (mean[phi], var[phi] and mixed fluctuations) for a TTBL.
!              This subroutine is called inside 'print_cf'. 
!              These statistics must be averaged with different flow 
!              realizations. 2nd order statistics must be finalized after 
!              the averaging procedure (subtraction of averages).  
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!       
subroutine print_mean_stats_scalar(ux,uy,uz,phi)
  
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use MPI
  
  use var,       only : ux2,uy2,uz2,phi2
  use param,     only : zero,itime
  use variables, only : nx,ny,nz
    
  implicit none
  
  ! Inputs
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)          ) :: ux,uy,uz
  real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
     
  ! Local variables
  real(mytype)                      :: den
  integer                           :: code,i,j,k,iunit
  character(99)                     :: filename 
  
  ! Mean statistics for scalar field
  real(mytype), dimension(ysize(2)) :: phi1meanH1,phi1meanHT  ! mean scalar field
  real(mytype), dimension(ysize(2)) :: phi2meanH1,phi2meanHT  ! scalar field variance
  
  ! Mixed fluctuations (velocity/scalar)
  real(mytype), dimension(ysize(2)) :: uphimeanH1,uphimeanHT  ! Mixed fluctuations <u'phi'>
  real(mytype), dimension(ysize(2)) :: vphimeanH1,vphimeanHT  ! Mixed fluctuations <v'phi'>
  real(mytype), dimension(ysize(2)) :: wphimeanH1,wphimeanHT  ! Mixed fluctuations <w'phi'>
  
  ! Set variables to zero
  phi1meanH1 = zero; phi1meanHT = zero
  phi2meanH1 = zero; phi2meanHT = zero  

  uphimeanH1 = zero; uphimeanHT = zero
  vphimeanH1 = zero; vphimeanHT = zero  
  wphimeanH1 = zero; wphimeanHT = zero
  
  ! Denominator of the divisions
  den = real(nx*nz,mytype)
      
  ! Transpose data to y-pencils 
  call transpose_x_to_y(ux,ux2)
  call transpose_x_to_y(uy,uy2)
  call transpose_x_to_y(uz,uz2)
  call transpose_x_to_y(phi(:,:,:,1), phi2(:,:,:,1))
    
  ! Summation over x and z directions
  do k=1,ysize(3)
      do i=1,ysize(1)
          do j=1,ysize(2)
          
              ! Mean scalar field          
              phi1meanH1(j)=phi1meanH1(j)+phi2(i,j,k,1)/den
              
              ! Variance
              phi2meanH1(j)=phi2meanH1(j)+phi2(i,j,k,1)*phi2(i,j,k,1)/den                                                             
          
              ! Mixed fluctuations
              uphimeanH1(j)=uphimeanH1(j)+ux2(i,j,k)*phi2(i,j,k,1)/den
              vphimeanH1(j)=vphimeanH1(j)+uy2(i,j,k)*phi2(i,j,k,1)/den
              wphimeanH1(j)=wphimeanH1(j)+uz2(i,j,k)*phi2(i,j,k,1)/den
                                      
          enddo          
      enddo
  enddo

  ! Summation over all MPI processes
  call MPI_REDUCE(phi1meanH1,phi1meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(phi2meanH1,phi2meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)

  call MPI_REDUCE(uphimeanH1,uphimeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(vphimeanH1,vphimeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  call MPI_REDUCE(wphimeanH1,wphimeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  
  ! Print mean statistics for scalar field
  if(nrank .eq. 0) then

      ! Write filename 
      write(filename, "('data/mean_stats_runtime/scalar/mean_stats_scalar_runtime-ts',I7.7,'.txt')") itime
            
      ! Open and write mean statistics for scalar field
      open(newunit=iunit, file=filename, status="unknown", form="formatted")

      ! Header
      write(iunit, *) 'Mean statistics for scalar field calculated runtime.'
      write(iunit, *) ' '
      write(iunit, *) 'Pay attention that these statistics need to be averaged later'
      write(iunit, *) 'with different flow realizations.'
      write(iunit, *) ' '    
      write(iunit, *) 'The finalization of all 2nd order statistics can be performed anyway in order to check'
      write(iunit, *) 'the correctness of the calculations.' 
      write(iunit, *) ' '
      write(iunit, '(5(A17, A1, 1X))') 'mean[phi](y,t,nr)', ',', &
                                       'var[phi](y,t,nr)' , ',', &
                                       "<u'phi'>(y,t,nr)" , ',', &
                                       "<v'phi'>(y,t,nr)" , ',', &
                                       "<w'phi'>(y,t,nr)"

      ! Write mean scalar field statistics, function of y-direction, time and specific realization
      do j=1,ny
          write(iunit, '(5(F17.9, A1, 1X))') phi1meanHT(j), ',', &  ! Mean scalar field
                                             phi2meanHT(j), ',', &  ! Scalar field variance
                                             uphimeanHT(j), ',', &  ! Mixed fluctuations <u'phi'>
                                             vphimeanHT(j), ',', &  ! Mixed fluctuations <v'phi'>
                                             wphimeanHT(j)          ! Mixed fluctuations <w'phi'>
      enddo
      
      close(iunit)
                      
  end if
  
end subroutine print_mean_stats_scalar

!-----------------------------------------------------------------------------!
! DESCRIPTION: Write an instantaneous plane with z-dir. normal of the scalar 
!              field for visualization.
!              Adapted from 'write_snapshot' and 'end_snapshot' subroutines.    
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------! 
subroutine write_scalar_plane_z(phi1,ux1,uz1,itime)
 
  use visu

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use decomp_2d_io, only : decomp_2d_start_io
  use param,        only : ioutput_plane
  use variables,    only : numscalar
  
  implicit none
  
  ! Inputs
  real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)),            intent(in) :: ux1,uz1
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
  call write_xdmf_header("planes", "phiplanez", trim(num))
   
  ! Write first scalar field
  call write_field(phi1(:,:,:,1), "planes", "phiplanez", trim(num), flush = .true.)
  
!--- End snapshot part ---!
  
  ! Write XDMF footer
  call write_xdmf_footer(ux1,uz1)
  
#ifdef ADIOS2
  call decomp_2d_end_io(io_name, "data")
#endif

!-------------------------!
  
  ! Switch back to 3D output for default snapshots
  output2D = 0
  
end subroutine write_scalar_plane_z

!-----------------------------------------------------------------------------!
! DESCRIPTION: Write an instantaneous plane with x-dir. normal of  
!              streamwise vorticity for visualization.
!              Adapted from 'write_snapshot' and 'end_snapshot' subroutines.    
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!   
subroutine write_vortx_plane_x(ux1,uy1,uz1,itime)
 
  use visu
  
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  
  use decomp_2d_io, only : decomp_2d_start_io
  use param,        only : ioutput_plane, zero
  use variables
  
  implicit none
  
  ! Inputs
  real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1,uy1,uz1
  integer,                                               intent(in) :: itime
  
  ! Locals
  character(len=32) :: num  ! taken from write_snapshot in visu module
  
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,di2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,di3
  real(mytype) :: lind
  
  ! x-derivatives
  call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0,lind)
  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,lind)
  call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1,lind)
  ! y-derivatives
  call transpose_x_to_y(ux1,td2)
  call transpose_x_to_y(uy1,te2)
  call transpose_x_to_y(uz1,tf2)
  call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,lind)
  call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0,lind)
  call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,lind)
  ! z-derivatives
  call transpose_y_to_z(td2,td3)
  call transpose_y_to_z(te2,te3)
  call transpose_y_to_z(tf2,tf3)
  call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,lind)
  call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1,lind)
  call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0,lind)
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
        
  ! Vorticity along x
  di1 = zero
     
  di1 = tf1 - th1  !dw/dy - dv/dz
        
  ! Switch to output2D with x-normal plane
  output2D = 1

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
  call write_xdmf_header("planes", "vortxplanex", trim(num))
       
  ! Write streamwise vorticity
  call write_field(di1(:,:,:), "planes", "vortxplanex", trim(num), flush = .true.)
  
!--- End snapshot part ---!
  
  ! Write XDMF footer
  call write_xdmf_footer(ux1,uz1)
  
#ifdef ADIOS2
  call decomp_2d_end_io(io_name, "data")
#endif

!-------------------------!
  
  ! Switch back to 3D output for default snapshots
  output2D = 0
  
end subroutine write_vortx_plane_x
!-----------------------------------------------------------------------------!


