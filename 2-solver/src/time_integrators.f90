!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!---------------------------------------------------------------------------!
! DESCRIPTION: This file contains subroutines for time-advancement of the 
!              Navier-Stokes equations via a fractional step method. 
!              Typically, you would use an explicit time advancement 
!              (Adams Bashforth / Runge Kutta). 
!              For wall-bounded flows, it is possible to use a semi-implicit 
!              approach for the viscous terms.
!---------------------------------------------------------------------------!

module time_integrators

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d

  implicit none

  private 
  
  ! Main subroutine for time-integration
  public  :: int_time

contains

  !-----------------------------------------------------------------------------!
  ! Main subroutine that collects all the other subroutines for
  ! time integration. It is called by the main file xcompact3d.f90.
  !-----------------------------------------------------------------------------!
  subroutine int_time(rho1, ux1, uy1, uz1, phi1, drho1, dux1, duy1, duz1, dphi1)

    use param,     only : zero, one
    use param,     only : ntime, nrhotime, ilmn, iscalar, ilmn_solve_temp,itimescheme
    use param,     only : iimplicit, sc_even
    use param,     only : primary_species, massfrac
    use param,     only : scalar_lbound, scalar_ubound
    use param,     only : nu0nu
    use variables, only : numscalar
    use var,       only : ta1, tb1

#ifdef DEBG 
    use tools, only : avg3d
#endif

    IMPLICIT NONE

    ! INPUT/OUTPUT
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime) :: drho1, dux1, duy1, duz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), ntime, numscalar) :: dphi1

    ! OUTPUT
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3), numscalar) :: phi1

    ! LOCAL
    integer :: is, i, j, k

#ifdef DEBG
    real(mytype) avg_param
    if (nrank .eq. 0) write(*,*)'## Init int_time'
#endif

    call int_time_momentum(ux1, uy1, uz1, dux1, duy1, duz1)

#ifdef DEBG
     avg_param = zero
     call avg3d (dux1, avg_param)
     if (nrank == 0) write(*,*)'## int_time dux1 ', avg_param
     avg_param = zero
     call avg3d (duy1, avg_param)
     if (nrank == 0) write(*,*)'## int_time duy1 ', avg_param
     avg_param = zero
     call avg3d (duz1, avg_param)
     if (nrank == 0) write(*,*)'## int_time duz1 ', avg_param
#endif

    !--- Commenting Low Mach Number (lmn) part to speed up the code ---!
    
    !IF (ilmn) THEN
    !   IF (ilmn_solve_temp) THEN
    !      CALL int_time_temperature(rho1, drho1, dphi1, phi1)
    !   ELSE
    !      CALL int_time_continuity(rho1, drho1)
    !   ENDIF
    !ENDIF

    IF (iscalar.NE.0) THEN
       
       !IF (ilmn.and.ilmn_solve_temp) THEN
       !   ! Compute temperature
       !   call calc_temp_eos(ta1, rho1(:,:,:,1), phi1, tb1, xsize(1), xsize(2), xsize(3))
       !ENDIF

       DO is = 1, numscalar
          IF (is.NE.primary_species) THEN
             IF (iimplicit.ge.1) then
                if (sc_even(is)) then
                   k = 1
                else
                   k = 0
                endif
                CALL intt(phi1(:,:,:,is), dphi1(:,:,:,:,is), npaire=k, isc=is)
             ELSE
                CALL intt(phi1(:,:,:,is), dphi1(:,:,:,:,is))
             ENDIF

             DO k = 1, xsize(3)
                DO j = 1, xsize(2)
                   DO i = 1, xsize(1)
                      phi1(i,j,k,is) = max(phi1(i,j,k,is),scalar_lbound(is))
                      phi1(i,j,k,is) = min(phi1(i,j,k,is),scalar_ubound(is))
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
       ENDDO

       !IF (primary_species.GE.1) THEN
       !   phi1(:,:,:,primary_species) = one
       !   DO is = 1, numscalar
       !      IF ((is.NE.primary_species).AND.massfrac(is)) THEN
       !         phi1(:,:,:,primary_species) = phi1(:,:,:,primary_species) - phi1(:,:,:,is)
       !      ENDIF
       !   ENDDO

       !   DO k = 1, xsize(3)
       !      DO j = 1, xsize(2)
       !         DO i = 1, xsize(1)
       !            phi1(i,j,k,primary_species) = max(phi1(i,j,k,primary_species),zero)
       !            phi1(i,j,k,primary_species) = min(phi1(i,j,k,primary_species),one)
       !         ENDDO
       !      ENDDO
       !   ENDDO
       !ENDIF

       !IF (ilmn.and.ilmn_solve_temp) THEN
       !   ! Compute rho
       !   call calc_rho_eos(rho1(:,:,:,1), ta1, phi1, tb1, xsize(1), xsize(2), xsize(3))
       !ENDIF
       
    ENDIF

#ifdef DEBG
    if (nrank .eq. 0) write(*,*)'## End  int_time'
#endif

  ENDSUBROUTINE int_time

  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: int_time_momentum
  ! DESCRIPTION: Integrates the momentum equations in time by calling time
  !              integrator.
  !      INPUTS: dux1, duy1, duz1 - the RHS(s) of the momentum equations
  !     OUTPUTS: ux1,   uy1,  uz1 - the intermediate momentum state.
  !       NOTES: This is integrating the MOMENTUM in time (!= velocity)
  !      AUTHOR: Paul Bartholomew
  !-----------------------------------------------------------------------------!
  subroutine int_time_momentum(ux1, uy1, uz1, dux1, duy1, duz1)

    use param
    use variables
    use var, only: px1, py1, pz1

    implicit none

    ! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1
    
    ! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1, duy1, duz1

    ! Semi-implicit diffusion
    if (iimplicit.ge.1) then
       
       ! x-dir. momentum
       if (itype .eq. itype_ttbl) then
           call intt(ux1, dux1, npaire=1, isc=0, forcing1=px1, wall_vel=uwall)
       else
           call intt(ux1, dux1, npaire=1, isc=0, forcing1=px1, wall_vel=zero)
       end if
       
       ! y-dir. momentum
       call intt(uy1, duy1, npaire=0, isc=0, forcing1=py1, wall_vel=zero)
       
       ! z-dir. momentum
       if (iswitch_wo .eq. 1) then
           call intt(uz1, duz1, npaire=1, isc=0, forcing1=pz1, wall_vel=span_vel)
       else 
           call intt(uz1, duz1, npaire=1, isc=0, forcing1=pz1, wall_vel=zero)
       end if
    else
       call intt(ux1, dux1)
       call intt(uy1, duy1)
       call intt(uz1, duz1)
    endif

  endsubroutine int_time_momentum

  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: int_time_continuity
  ! DESCRIPTION: Integrates the continuity (aka density transport) equation in
  !              time
  !      INPUTS: drho1 - the RHS(s) of the continuity equation.
  !     OUTPUTS:  rho1 - the density at new time.
  !      AUTHOR: Paul Bartholomew
  !-----------------------------------------------------------------------------!
  subroutine int_time_continuity(rho1, drho1)

    use param
    use variables

    implicit none

    integer :: it, i, j, k
    real(mytype) :: rhomin, rhomax

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: drho1

    !! First, update old density / store old transients depending on scheme
    if (itimescheme.lt.5) then
       !! Euler/AB - Store old density values
       do it = nrhotime, 2, -1
          rho1(:,:,:,it) = rho1(:,:,:,it-1)
       enddo
    elseif (itimescheme.eq.5) then
       !! RK3 - Stores old transients
       if (itr.eq.1) then
          do it = nrhotime, 2, -1
             rho1(:,:,:,it) = rho1(:,:,:,it-1)
          enddo
          rho1(:,:,:,2) = drho1(:,:,:,1)
       endif
    else
       if (nrank  == 0) then
          write(*,*) "int_time_continuity not implemented for itimescheme", itimescheme
          stop
       endif
    endif

    !! Now we can update current density
    call intt(rho1(:,:,:,1), drho1)

    !! Enforce boundedness on density
    if (ilmn_bound) then
       rhomin = min(dens1, dens2)
       rhomax = max(dens1, dens2)
       do k = 1, xsize(3)
          do j = 1, xsize(2)
             do i = 1, xsize(1)
                rho1(i, j, k, 1) = max(rho1(i, j, k, 1), rhomin)
                rho1(i, j, k, 1) = min(rho1(i, j, k, 1), rhomax)
             enddo
          enddo
       enddo
    endif

  endsubroutine int_time_continuity

  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: int_time_temperature
  ! DESCRIPTION: Integrates the temperature equation in time
  !      INPUTS: drho1 - the RHS(s) of the temperature equation.
  !     OUTPUTS:  rho1 - the density at new time.
  !      AUTHOR: Paul Bartholomew
  !-----------------------------------------------------------------------------!
  subroutine int_time_temperature(rho1, drho1, dphi1, phi1)

    use param
    use variables

    use navier, only : lmn_t_to_rho_trans
    use var, only : tc1, tb1

    implicit none

    integer :: it, i, j, k

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: drho1

    !! First, update old density / store old transients depending on scheme
    if (itimescheme.lt.5) then
       !! Euler/AB - Store old density values
       do it = nrhotime, 2, -1
          rho1(:,:,:,it) = rho1(:,:,:,it-1)
       enddo
    elseif (itimescheme.eq.5) then
       !! RK3 - Stores old transients
       if (itr.eq.1) then
          do it = nrhotime, 2, -1
             rho1(:,:,:,it) = rho1(:,:,:,it-1)
          enddo

          !! Convert temperature transient to density transient and store it.
          call lmn_t_to_rho_trans(rho1(:,:,:,2), drho1(:,:,:,1), rho1(:,:,:,1), dphi1, phi1)
       endif
    else
       if (nrank==0) then
          write(*,*) "int_time_continuity not implemented for itimescheme", itimescheme
          stop
       endif
    endif

    !!-------------------------------------------------------------------
    !! XXX We are integrating the temperature equation - get temperature
    !!-------------------------------------------------------------------
    call calc_temp_eos(tc1, rho1(:,:,:,1), phi1, tb1, xsize(1), xsize(2), xsize(3))

    !! Now we can update current temperature
    call intt(tc1, drho1)

    !! Temperature >= 0
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
             tc1(i,j,k) = max(tc1(i,j,k), zero)
          enddo
       enddo
    enddo

    !!-------------------------------------------------------------------
    !! XXX We are integrating the temperature equation - get back to rho
    !!-------------------------------------------------------------------
    call calc_rho_eos(rho1(:,:,:,1), tc1, phi1, tb1, xsize(1), xsize(2), xsize(3))

  endsubroutine int_time_temperature

  !-----------------------------------------------------------------------------!
  ! Perform time-integration of a generic field (velocity, scalars, etc.) with
  ! call to an external subroutine for implicit y-diffusion if necessary.
  !-----------------------------------------------------------------------------!
  subroutine intt(var1,dvar1,npaire,isc,forcing1,wall_vel)

   use MPI
   use param
   use variables
   use ydiff_implicit, only : inttimp

#ifdef DEBG 
   use tools, only : avg3d
#endif

   implicit none

   ! INPUT / OUTPUT
   real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: var1
   real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dvar1

   ! INPUTS
   real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(in), optional :: forcing1
   integer, intent(in), optional :: npaire, isc
   
   ! Wall-velocity BC, in order to differentiate between different directions
   real(mytype), intent(in), optional :: wall_vel

   ! LOCAL
   integer :: is, code, ierror

#ifdef DEBG 
   real(mytype) avg_param
#endif

#ifdef DEBG
   avg_param = zero
   call avg3d (var1, avg_param)
   if (nrank == 0) write(*,*)'## SUB intt VAR var1 (start) AVG ', avg_param
   avg_param = zero
   call avg3d (dvar1(:,:,:,1), avg_param)
   if (nrank == 0) write(*,*)'## SUB intt VAR dvar1(1) (start) AVG ', avg_param
   avg_param = zero
   call avg3d (dvar1(:,:,:,2), avg_param)
   if (nrank == 0) write(*,*)'## SUB intt VAR dvar1(2) (start) AVG ', avg_param
#endif

   ! (Semi) implicit Y-diffusion 
   ! (all integrations are performed inside inttimp subroutine)
   if (iimplicit.ge.1) then

      if (present(isc)) then
         is = isc
      else
         is = 0
      endif
      
      if (present(npaire).and.present(forcing1)) then
         call inttimp(var1, dvar1, npaire=npaire, isc=is, forcing1=forcing1, wall_vel=wall_vel)
      else if (present(npaire)) then
         call inttimp(var1, dvar1, npaire=npaire, isc=is, wall_vel=wall_vel)
      else
         if (nrank  == 0) write(*,*) "Error in intt call."
         call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
      endif
   
   !--- Fully-explicit time integration ---!
   
   ! Euler
   elseif (itimescheme.eq.1) then
      var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
   
   ! Adam-Bashforth second order (AB2)
   elseif(itimescheme.eq.2) then
      
      ! Do first time step with Euler
      if(itime.eq.1.and.irestart.eq.0) then
         var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
      else
         var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+var1(:,:,:)
      endif
      dvar1(:,:,:,2)=dvar1(:,:,:,1)

   ! Adams-Bashforth third order (AB3)
   elseif(itimescheme.eq.3) then
     
      ! Do first time step with Euler
      if(itime.eq.1.and.irestart.eq.0) then
         var1(:,:,:)=dt*dvar1(:,:,:,1)+var1(:,:,:)
      elseif(itime.eq.2.and.irestart.eq.0) then
         ! Do second time step with AB2
         var1(:,:,:)=onepfive*dt*dvar1(:,:,:,1)-half*dt*dvar1(:,:,:,2)+var1(:,:,:)
         dvar1(:,:,:,3)=dvar1(:,:,:,2)
      else
         ! Finally using AB3
         var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+cdt(itr)*dvar1(:,:,:,3)+var1(:,:,:)
         dvar1(:,:,:,3)=dvar1(:,:,:,2)
      endif
      dvar1(:,:,:,2)=dvar1(:,:,:,1)
   
   ! Adams-Bashforth fourth order (AB4)
   elseif(itimescheme.eq.4) then

      if (nrank==0) then
         write(*,*) "AB4 not implemented!"
         stop
      endif

   ! Runge-Kutta (low storage) RK3
   elseif(itimescheme.eq.5) then
      if(itr.eq.1) then
         var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
      else
         var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+var1(:,:,:)
      endif
      dvar1(:,:,:,2)=dvar1(:,:,:,1)
   
   ! Runge-Kutta (low storage) RK4
   elseif(itimescheme.eq.6) then

      if (nrank==0) then
         write(*,*) "RK4 not implemented!"
         STOP
      endif

   else

      if (nrank==0) then
         write(*,*) "Unrecognised itimescheme: ", itimescheme
         STOP
      endif

   endif

#ifdef DEBG
   avg_param = zero
   call avg3d (var1, avg_param)
   if (nrank == 0) write(*,*)'## SUB intt VAR var1 AVG ', avg_param
   avg_param = zero
   call avg3d (dvar1(:,:,:,1), avg_param)
   if (nrank == 0) write(*,*)'## SUB intt VAR dvar1 AVG ', avg_param
   if (nrank   ==  0) write(*,*)'# intt done'
#endif

   return

 end subroutine intt

end module time_integrators
