!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!-----------------------------------------------------------------------------!
! DESCRIPTION: This file contains pointers to each of the different flow 
!              configurations (the BC files) for the initialisation, 
!              imposition of the boundary conditions, post-processing and 
!              visualisation. It means that the bulk of the code does not 
!              need to be changed for each flow configuration.
!-----------------------------------------------------------------------------!

module case

  use decomp_2d_constants
  use decomp_2d_mpi  
  use decomp_2d
     
  use param
  use variables
  use var, only : nzmsize
  
  use channel
  use temporal_tbl

  implicit none

  logical :: case_visu_init = .false.
  
  private ! All functions/subroutines private by default
  
  public :: init,                 &
            boundary_conditions,  &
            preprocessing,        &
            postprocessing,       &
            visu_case_init,       &
            visu_case,            &
            momentum_forcing,     &
            scalar_forcing,       &
            set_fluid_properties, &
            test_flow 

contains
  !-----------------------------------------------------------------------------!
  ! DESCRIPTION: Initialization of flow cases 
  !              (imposition of initial conditions, ICs).
  !-----------------------------------------------------------------------------!
  subroutine init (rho1, ux1, uy1, uz1, ep1, phi1, drho1, dux1, duy1, duz1, dphi1, &
                   pp3,  px1, py1, pz1)

    implicit none
    
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1, py1, pz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1,drho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3

    integer :: it, is

    ! Zero out the pressure field
    pp3(:,:,:,1) = zero
    px1(:,:,:)   = zero
    py1(:,:,:)   = zero
    pz1(:,:,:)   = zero

    ! Default density and pressure0 to one
    pressure0 = one
    rho1(:,:,:,:) = one

    ! Initialize Channel flow
    if (itype.eq.itype_channel) then
    
       call init_channel (ux1, uy1, uz1, ep1, phi1)
    
    ! Initialize TTBL
    elseif (itype.eq.itype_ttbl) then

       call init_temporal_tbl (ux1, uy1, uz1, phi1)

    else
       
       ! Exit if an unknown flow type has been selected 
       if (nrank.eq.0) then
          print *, "ERROR: Unknown itype: ", itype
          stop
       endif

    endif

    ! Setup old arrays
    do it = 1, ntime
       drho1(:,:,:,it) = rho1(:,:,:,1)
       dux1 (:,:,:,it) = ux1 (:,:,:)
       duy1 (:,:,:,it) = uy1 (:,:,:)
       duz1 (:,:,:,it) = uz1 (:,:,:)
    enddo

    do it = 2, nrhotime
       rho1(:,:,:,it) = rho1(:,:,:,1)
    enddo

    do is = 1, numscalar
       do it = 1, ntime
          dphi1(:,:,:,it,is) = phi1(:,:,:,is)
       enddo
    enddo

  end subroutine init
  !-----------------------------------------------------------------------------!
  ! DESCRIPTION: Imposition of boundary conditions, BCs.
  !-----------------------------------------------------------------------------!
  subroutine boundary_conditions (rho,ux,uy,uz,phi,ep)

    implicit none
    
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho
    
    ! Channel flow
    if (itype.eq.itype_channel) then

       ! Calculate the spanwise wall oscillations
       if (iswitch_wo .eq. 1) then
           call spanwise_wall_oscillations (ux,uz)
       end if

       call boundary_conditions_channel (ux, uy, uz, phi)
    
    ! TTBL
    elseif (itype.eq.itype_ttbl) then
       
       ! Calculate the spanwise wall oscillations
       if (iswitch_wo .eq. 1) then
           call spanwise_wall_oscillations (ux,uz)
       end if
       
       call boundary_conditions_ttbl (phi)   
       
    endif

  end subroutine boundary_conditions
  !-----------------------------------------------------------------------------!
  subroutine preprocessing(rho1, ux1, uy1, uz1, pp3, phi1, ep1)

    use decomp_2d, only : xsize, ph1
    use visu,      only : write_snapshot
    use var,       only : itime, numscalar, nrhotime, npress

    implicit none
    
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime), intent(in) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ep1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3

  end subroutine preprocessing
  !-----------------------------------------------------------------------------!
  ! DESCRIPTION: Saving of snapshots (velocity, pressure, scalar fields) and,
  !              if requested, save flowcase-specific fields for visualization
  !              (e.g. Q-criterion).
  !-----------------------------------------------------------------------------!
  subroutine postprocessing(rho1, ux1, uy1, uz1, pp3, phi1, ep1)

    use decomp_2d,   only : xsize, ph1
    use visu,        only : write_snapshot, end_snapshot
    use var,         only : itime, numscalar, nrhotime, npress
    
    implicit none
    
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),           intent(in) :: ux1, uy1, uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar), intent(in) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime),  intent(in) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),           intent(in) :: ep1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp3

    integer :: j
    character(len=32) :: num
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: T ! FIXME This can be huge
    
    T = zero

    ! Recover temperature when decomposed (pressure to be recovered externally)
    T = phi1

    ! Writing snapshots only if some conditions are met
    if ( ( isnap .ne. 0 )                                        .and. &  ! Snapshots enabled by the switcher 'isnap'
         ((itime .eq. ifirst) .or. (mod(itime, ioutput) .eq. 0)) .and. &  ! Save at the first time step or at every 'ioutput' time-steps  
         ( itime .ge. start_output)                              .and. &  ! Save snapshots if time-step is greater or equal to 'start_output' 
         ( itime .le. end_output  )                                    &  ! Save snapshots if time-step is less or equal to 'end_output' 
       ) then
       
       ! Write snapshot
       call write_snapshot(rho1, ux1, uy1, uz1, pp3, T, ep1, itime, num)

       ! XXX: Ultimate goal for ADIOS2 is to pass do all postproc online - do we need this?
       !      Currently, needs some way to "register" variables for IO
       
       ! Save case-specific field for visualization (e.g. Q-criterion)
       if (ivisu .eq. 1) call visu_case(rho1, ux1, uy1, uz1, pp3, T, ep1, num)

       call end_snapshot(ux1, uz1, itime, num)
                     
    end if

    call postprocess_case(rho1, ux1, uy1, uz1, pp3, T, ep1)

  end subroutine postprocessing
  !-----------------------------------------------------------------------------!
  ! DESCRIPTION: Flowcase-specific post-processing if needed.
  !-----------------------------------------------------------------------------!
  subroutine postprocess_case(rho,ux,uy,uz,pp,phi,ep)

    use param, only : npress

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ep
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress), intent(in) :: pp

    if (itype.eq.itype_channel) then

       call postprocess_channel (ux, uy, uz, pp, phi, ep)
        
    elseif (itype.eq.itype_ttbl) then

       call postprocess_ttbl (ux, uy, uz, pp, phi, ep)
    
    endif

  end subroutine postprocess_case
  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: visu_case_init
  !      AUTHOR: PB
  ! DESCRIPTION: Initialise case-specific visualization.
  !-----------------------------------------------------------------------------!
  subroutine visu_case_init

    implicit none
    
    if (itype .eq. itype_channel) then

       call visu_channel_init(case_visu_init)
    
    else if (itype .eq. itype_ttbl) then

       call visu_ttbl_init(case_visu_init)

    end if
    
  end subroutine visu_case_init
  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: visu_case
  !      AUTHOR: CF
  ! DESCRIPTION: Call case-specific visualization.
  !-----------------------------------------------------------------------------!
  subroutine visu_case(rho1,ux1,uy1,uz1,pp3,phi1,ep1,num)

    use param, only : npress

    implicit none
    
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype), intent(in), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), intent(in), dimension(xsize(1),xsize(2),xsize(3)) :: ep1
    
    character(len=32), intent(in) :: num

    logical :: called_visu = .false.
      
    if (itype.eq.itype_channel) then

       call visu_channel(ux1, uy1, uz1, pp3, phi1, ep1, num)
       called_visu = .true.
       
    elseif (itype.eq.itype_ttbl) then

       call visu_ttbl(ux1, uy1, uz1, pp3, phi1, ep1, num)
       called_visu = .true.
       
    endif

    if (called_visu .and. (.not. case_visu_init)) then

       print *, "ERROR: tried to run case-specific visu without initialisation!"
       print *, "       See the TGV case initialisation for example."
       stop
       
    endif

  end subroutine visu_case
  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: momentum_forcing
  !      AUTHOR: Paul Bartholomew
  ! DESCRIPTION: Calls case-specific forcing functions for the
  !              momentum equations.
  !-----------------------------------------------------------------------------!
  subroutine momentum_forcing(dux1, duy1, duz1, rho1, ux1, uy1, uz1, phi1)

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), numscalar) :: phi1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    if (itype.eq.itype_channel) then

       call momentum_forcing_channel(dux1, duy1, duz1, ux1, uy1, uz1)

    endif

  end subroutine momentum_forcing
  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: scalar_forcing
  !      AUTHOR: Kay Schäfer
  ! DESCRIPTION: Calls case-specific forcing functions for the
  !              scalar transport equations.
  !-----------------------------------------------------------------------------!
  subroutine scalar_forcing(dphi1, rho1, ux1, uy1, uz1, phi1)

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, phi1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime) :: dphi1

  end subroutine scalar_forcing
  !-----------------------------------------------------------------------------!
  subroutine set_fluid_properties(rho1, mu1)

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mu1

  endsubroutine set_fluid_properties
  !-----------------------------------------------------------------------------!
  ! DESCRIPTION: Check if velocity or scalar fields diverged. 
  !              Compute and show:
  !               - CFL (or Co):  Courant number;
  !               - D:            Numerical Fourier number;
  !               - Pé (or Re_c): Numerical Péclet or Reynolds cell number;
  !               - S:            Stability parameter (Thompson et al. (1985)).
  !-----------------------------------------------------------------------------!
  subroutine test_flow(rho1,ux1,uy1,uz1,phi1,ep1,drho1,divu3)

    use decomp_2d
    use param

    use navier, only : divergence

    use var,    only : numscalar, dv3
    use tools,  only : test_speed_min_max, compute_cfl, &
                       test_scalar_min_max, compute_cfldiff, &
                       compute_reynolds_cell, compute_stab_param
    
    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1, ep1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime), intent(in) :: drho1
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: divu3

    if ((mod(itime,ilist)==0 .or. itime == ifirst .or. itime == ilast)) then
       call divergence(dv3,rho1,ux1,uy1,uz1,ep1,drho1,divu3,2)
       call test_speed_min_max(ux1,uy1,uz1)
       call compute_cfl(ux1,uy1,uz1)
       call compute_cfldiff()
       call compute_reynolds_cell(ux1,uy1,uz1)
       call compute_stab_param(ux1,uy1,uz1)
       if (iscalar==1) call test_scalar_min_max(phi1)
    endif

  end subroutine test_flow
  !-----------------------------------------------------------------------------!
  
end module case


