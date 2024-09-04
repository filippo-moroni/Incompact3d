!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

program xcompact3d

  use var
  use case

  use transeq,          only : calculate_transeq_rhs
  use time_integrators, only : int_time
  use navier,           only : velocity_to_momentum, momentum_to_velocity, pre_correc, &
                               calc_divu_constraint, solve_poisson, cor_vel
  use tools,            only : restart, simu_stats, apply_spatial_filter
  use ibm_param
  use ibm,              only : body
  use genepsi,          only : genepsi3d
    
  implicit none
  
  call init_xcompact3d()

  t = t0

  do itime=ifirst,ilast
       
     !t=itime*dt
     !t=t0 + (itime0 + itime + 1 - ifirst)*dt
      
     t = t + dt
     
     call simu_stats(2)
 
     ! Sub-time steps cycle (for RK schemes)
     do itr=1,iadvance_time

        call set_fluid_properties(rho1,mu1)
        
        ! BCs are updated at each sub-time step        
        call boundary_conditions(rho1,ux1,uy1,uz1,phi1,ep1)

        !if (imove.eq.1) then ! update epsi for moving objects
        !  if ((iibm.eq.2).or.(iibm.eq.3)) then
        !     call genepsi3d(ep1)
        !  else if (iibm.eq.1) then
        !     call body(ux1,uy1,uz1,ep1)
        !  endif
        !endif
        
        call calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)
        
#ifdef DEBG
        call check_transients()
#endif
        
        !if (ilmn) then
        !   !! XXX N.B. from this point, X-pencil velocity arrays contain momentum (LMN only).
        !   call velocity_to_momentum(rho1,ux1,uy1,uz1)
        !endif

        call int_time(rho1,ux1,uy1,uz1,phi1,drho1,dux1,duy1,duz1,dphi1)
        call pre_correc(ux1,uy1,uz1,ep1)

        call calc_divu_constraint(divu3,rho1,phi1)
        call solve_poisson(pp3,px1,py1,pz1,rho1,ux1,uy1,uz1,ep1,drho1,divu3)
        call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

        !if (ilmn) then
        !   call momentum_to_velocity(rho1,ux1,uy1,uz1)
        !   !! XXX N.B. from this point, X-pencil velocity arrays contain velocity (LMN only).
        !   !! Note - all other solvers work on velocity always
        !endif
                
     enddo ! End of sub-time steps cycle
     
     ! Calculation of numerics-related parameters (CFL, Pé, S) and if the simulation diverged
     call test_flow(rho1,ux1,uy1,uz1,phi1,ep1,drho1,divu3)

     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,1)

     call simu_stats(3)

     call postprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
     
     ! Print of cf for monitoring
     if ((mod(itime, ioutput_cf) .eq. 0) .or. (itime.eq.ifirst) .and. (itime .ge. start_output)) then
         call print_cf(ux1,uz1,phi1)
     end if
     
     ! Save planes of relevant quantities for quick visualization and low memory requirements
     if (mod(itime, ioutput_plane) .eq. 0 .or. (itime.eq.ifirst) .and. (itime .ge. start_output)) then
     
     ! Save a streamwise vorticity plane with x-normal
     call write_vortx_plane_x(ux1,uy1,uz1,itime)
         
         ! Save a scalar plane with z-normal 
         if (iscalar .eq. 1) then
             call write_scalar_plane_z(phi1,ux1,uz1,itime)       
         end if
     end if
     
  enddo ! End time loop

  call finalise_xcompact3d()

end program xcompact3d

!-----------------------------------------------------------------------------!
subroutine init_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_io,      only : decomp_2d_io_init
  USE decomp_2d_poisson, only : decomp_2d_poisson_init
  use case

  use var

  use navier, only : calc_divu_constraint
  use tools,  only : test_speed_min_max, test_scalar_min_max, &
                     restart, simu_stats, compute_cfldiff

  use param, only : ilesmod, jles,itype
  use param, only : irestart, nvisu, ilist

  use variables, only : nx, ny, nz, nxm, nym, nzm
  use variables, only : p_row, p_col
  use variables, only : nstat, nprobe

  use les, only: init_explicit_les

  use visu, only : visu_init, visu_ready

  use genepsi, only : genepsi3d, epsi_init
  use ibm, only : body

  implicit none

  integer :: ierr

  integer :: nargin, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase
    
  ! Initialise MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! Reading input file 
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     if (nrank==0) write(*,*) 'Xcompact3d is run with the default file -->', trim(InputFN)
  elseif (nargin >= 1) then
     call get_command_argument(1,InputFN,FNLength,status)
     back=.true.
     FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
     DecInd=index(FNBase,'.',back)
     if (DecInd >1) then
        FNBase=FNBase(1:(DecInd-1))
     end if
     if (nrank==0) write(*,*) 'Xcompact3d is run with the provided file -->', trim(InputFN)
  endif

#ifdef ADIOS2
  if (nrank .eq. 0) then
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
     print *, " WARNING: Running Xcompact3d with ADIOS2"
     print *, "          this is currently experimental"
     print *, "          for safety of results it is recommended"
     print *, "          to run the default build as this feature"
     print *, "          is developed. Thank you for trying it."
     print *, " WARNING === WARNING === WARNING === WARNING === WARNING"
  endif
#endif
  
  call parameter(InputFN)

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_io_init()
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)    !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true
  
  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  call decomp_info_init(nxm, ny, nz, ph4)
  
  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)
  call decomp_info_init(nxm, nym, nz, ph3)

  call init_variables()

  call schemes()

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

  if (ilesmod.ne.0) then
     if (jles.gt.0)  call init_explicit_les()
  endif

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if (iibm.eq.1) then
     call epsi_init(ep1)
     call body(ux1,uy1,uz1,ep1)
  endif

  ! initialise visu
  if (ivisu.ne.0) then
     call visu_init()
     call visu_case_init() !! XXX: If you get error about uninitialised IO, look here.
                           !! Ensures additional case-specific variables declared for IO
     call visu_ready()
  end if    
  
  if (irestart==0) then
     call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
     itime = 0
     call preprocessing(rho1,ux1,uy1,uz1,pp3,phi1,ep1)
  else
     itr=1
     call restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3(:,:,:,1),phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,0)
  endif

  if ((iibm.eq.2).or.(iibm.eq.3)) then
     call genepsi3d(ep1)
  else if ((iibm.eq.1).or.(iibm.eq.3)) then
     call body(ux1,uy1,uz1,ep1)
  endif

  if (mod(itime, ilist) == 0 .or. itime == ifirst) then
     call test_speed_min_max(ux1,uy1,uz1)
     if (iscalar==1) call test_scalar_min_max(phi1)
  endif

  call simu_stats(1)

  call calc_divu_constraint(divu3, rho1, phi1)

endsubroutine init_xcompact3d

!-----------------------------------------------------------------------------!
subroutine finalise_xcompact3d()

  use MPI
  use decomp_2d
  use decomp_2d_io, only : decomp_2d_io_finalise

  use tools, only : simu_stats
  use param, only : itype, dt, ifirst, ilast, t
  use visu,  only : visu_finalise

  implicit none

  integer       :: ierr, iunit
  logical       :: exists
  character(99) :: filename, ts_index
    
  !--- Create or open a file to store the dt and ts of stop ---! 
  if(nrank.eq.0) then
 
      ! Write filename
      write(filename,"('data/monitoring/dt_history.txt')") 
      
      inquire(file=filename, exist=exists)
      if (exists) then
          open(newunit=iunit, file=filename, status="old", position="append", action="write")
          write(iunit, '(F8.6,A,I12,A,I12,A,F10.4)')  dt,  ',',  ifirst,  ',',  ilast, ',', t
      else
          open(newunit=iunit, file=filename, status="new", action="write")
          write(iunit, '(A8,A,A12,A,A12,A,A10)') 'dt', ',', 'ifirst', ',', 'ilast', ',', 'time_unit'
          write(iunit, '(F8.6,A,I12,A,I12,A,F10.4)')  dt,  ',',  ifirst,  ',',  ilast, ',', t
      end if
      close(iunit)
  end if
    
  call simu_stats(4)
  call visu_finalise()
  call decomp_2d_io_finalise()
  call decomp_2d_finalize
  CALL MPI_FINALIZE(ierr)

endsubroutine finalise_xcompact3d

!-----------------------------------------------------------------------------!
! Subroutine used in debugging mode in the main file.
!-----------------------------------------------------------------------------!
subroutine check_transients()

  use decomp_2d, only : mytype

  use var
  use tools, only : avg3d
  
  implicit none

  real(mytype) avg_param
  
  avg_param = zero
  call avg3d (dux1, avg_param)
  if (nrank == 0) write(*,*)'## Main dux1 ', avg_param
  avg_param = zero
  call avg3d (duy1, avg_param)
  if (nrank == 0) write(*,*)'## Main duy1 ', avg_param
  avg_param = zero
  call avg3d (duz1, avg_param)
  if (nrank == 0) write(*,*)'## Main duz1 ', avg_param
  
end subroutine check_transients
