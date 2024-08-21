!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!-----------------------------------------------------------------------------!                                                                               !
!  SUBROUTINE: parameter                                                       
! DESCRIPTION: Reads the input.i3d file and sets the parameters of the         
!              simulation.                                                     
!      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>                                                                                      
!-----------------------------------------------------------------------------!
subroutine parameter(input_i3d)

  use MPI
  
  use iso_fortran_env

  use param
  use variables
  use complex_geometry
  use decomp_2d
  use ibm_param

  use visu,  only : output2D
  use tools, only : program_header

  implicit none

  character(len=80), intent(in) :: input_i3d
  real(mytype) :: theta,cfl,cf2
  integer :: longueur,impi,j,is,total

  NAMELIST /BasicParam/ itype, p_row, p_col, nx, ny, nz,          &
                        istret, beta, xlx, yly, zlz,              &
                        re, dt, ifirst, ilast, irestart,          &
                        numscalar, iin, init_noise,               &
                        nclx1, nclxn, ncly1, nclyn, nclz1, nclzn, &
                        iibm, ilmn, ilesmod, iscalar,             &
                        ipost, ifilter, C_filter,                 &
                        gravx, gravy, gravz
       
  NAMELIST /NumOptions/ ifirstder, isecondder, ipinter, itimescheme, iimplicit, &
                        nu0nu, cnu
                        
  NAMELIST /InOutParam/ icheckpoint, ioutput, ioutput_cf, ioutput_plane, &
                        ilist, isnap, ivisu, nvisu, output2D, start_output
  
  NAMELIST /AdditionalControls/ iswitch_wo
  NAMELIST /WallOscillations/ ifeedback_control, a_wo, t_wo, in_phase
  
  NAMELIST /ChannelParam/ cpg, idir_stream, wrotation, spinup_time
  NAMELIST /TemporalTBLParam/ uwall, twd, uln, lln, phiwall 
  
  NAMELIST /ScalarParam/ sc, ri, uset, cp, &
                         nclxS1, nclxSn, nclyS1, nclySn, nclzS1, nclzSn, &
                         scalar_lbound, scalar_ubound, sc_even, sc_skew, &
                         alpha_sc, beta_sc, g_sc, Tref
       
  NAMELIST /LESModel/ jles, smagcst, smagwalldamp, nsmag, walecst, maxdsmagcst, iwall
  
  NAMELIST /IbmStuff/ cex, cey, cez, ra, nobjmax, nraf, npif, izap, ianal, &
                      imove, thickness, chord, omega, ubcx, ubcy, ubcz, rads, c_air
  
  NAMELIST /LMN/ dens1, dens2, prandtl, ilmn_bound, ivarcoeff, ilmn_solve_temp, &
                 massfrac, mol_weight, imultispecies, primary_species, &
                 Fr, ibirman_eos
  
  ! Not used at the moment
  NAMELIST /Statistics/ nstat, initstat
  NAMELIST /ExtraNumControl/ icfllim, cfl_limit
  
#ifdef DEBG
  if (nrank == 0) write(*,*) '# parameter start'
#endif

  ! Program header with credits
  call program_header()

  call parameter_defaults()

  ! Read parameters
  open(10, file=input_i3d)

  ! These are the 'essential' parameters
  read(10, nml=BasicParam); rewind(10)
  read(10, nml=NumOptions); rewind(10)
  read(10, nml=InOutParam); rewind(10)
  
  ! Additional controls
  read(10, nml=AdditionalControls); rewind(10)
  
  ! Controls for wall oscillations
  if(iswitch_wo .eq. 1) then
      read(10, nml=WallOscillations); rewind(10); 
  end if
  
  ! Read parameters for Channel case
  if (itype.eq.itype_channel) then
     read(10, nml=ChannelParam); rewind(10);   
  end if
  
  ! Read parameters for temporal TBL case
  if (itype.eq.itype_ttbl) then
     read(10, nml=TemporalTBLParam); rewind(10);   
  end if
  
  ! Immersed boundary method
  if (iibm.ne.0) then
      read(10, nml=ibmstuff); rewind(10)
  endif
    
  ! Set Scalar BCs same as fluid (may be overridden) [DEFAULT]
  nclxS1 = nclx1; nclxSn = nclxn
  nclyS1 = ncly1; nclySn = nclyn
  nclzS1 = nclz1; nclzSn = nclzn
  
  if (numscalar.ne.0) then
     iscalar = 1

     ! Allocate scalar arrays and set sensible defaults
     allocate(massfrac(numscalar))
     allocate(mol_weight(numscalar))
     massfrac(:) = .FALSE.
     mol_weight(:) = one
     allocate(sc(numscalar), ri(numscalar), uset(numscalar), cp(numscalar))
     ri(:) = zero
     uset(:) = zero
     cp(:) = zero
     if (iimplicit.gt.0) then
        allocate(xcst_sc(numscalar))
        xcst_sc(:) = zero
        allocate(alpha_sc(numscalar,2), beta_sc(numscalar,2), g_sc(numscalar,2))
        ! Default scalar BC : dirichlet BC, zero value
        alpha_sc = one
        beta_sc = zero
        g_sc = zero
     endif

     ! In case of symmetry, scalars are even by default
     allocate(sc_even(numscalar))
     sc_even(:) = .true.

     ! Skew-symmetric convection of scalars, off by default
     allocate(sc_skew(numscalar))
     sc_skew(:) = .false.

     allocate(scalar_lbound(numscalar), scalar_ubound(numscalar))
     scalar_lbound(:) = -huge(one)
     scalar_ubound(:) = huge(one)
  endif

  if (ilmn) then
     if (istret.ne.0) then
        if (nrank.eq.0) then
           print *, "WARNING: LMN solver does not currently support stretching!"
           stop
        endif
     endif
 
  read(10, nml=LMN); rewind(10)
     do is = 1, numscalar
        if (massfrac(is)) then
           imultispecies = .TRUE.
        endif
     enddo

     if (imultispecies) then
        if (primary_species.lt.1) then
           if (nrank==0) then
              write(*,*)  "Error: you must set a primary species for multispecies flow"
              write(*,*)  "       solver will enforce Y_p = 1 - sum_s Y_s, s != p."
              stop
           endif
        else if (.not.massfrac(primary_species)) then
           if (nrank==0) then
              write(*,*)  "Error: primary species must be a massfraction!"
           endif
        endif
     endif
  endif
  
  if (numscalar.ne.0) then
     read(10, nml=ScalarParam); rewind(10)
  endif
    
  if(ilesmod.ne.0) then
     read(10, nml=LESModel); rewind(10)
  endif
    
  !read(10, nml=Statistics); rewind(10)
    
  ! Read extra numerics control (Adjustable time-step)
  !read(10, nml=ExtraNumControl); rewind(10);
       
  close(10)

  ! allocate(sc(numscalar),cp(numscalar),ri(numscalar),group(numscalar))

  if (nclx1.eq.0.and.nclxn.eq.0) then
     nclx=.true.
     nxm=nx
  else
     nclx=.false.
     nxm=nx-1
  endif
  if (ncly1.eq.0.and.nclyn.eq.0) then
     ncly=.true.
     nym=ny
  else
     ncly=.false.
     nym=ny-1
  endif
  if (nclz1.eq.0.and.nclzn.eq.0) then
     nclz=.true.
     nzm=nz
  else
     nclz=.false.
     nzm=nz-1
  endif

  ! Mesh spacings
  dx=xlx/real(nxm,mytype)
  dy=yly/real(nym,mytype)
  dz=zlz/real(nzm,mytype)

  dx2 = dx * dx
  dy2 = dy * dy
  dz2 = dz * dz

  xnu = one/re
  
  ! Calculation of correct viscosity and Re numbers for a channel flow case
  if (itype .eq. itype_channel) then
  
      ! Constant pressure gradient, re = Re_tau --> used to compute Re_centerline
      if (cpg) then
  
          ! Calculate Reynolds centerline of a laminar Poiseuille flow (re = Re_tau in this case)
          re_cent = (re/0.116_mytype)**(1.0_mytype/0.88_mytype)
    
          ! Viscosity based on Re_cent to keep same scaling as CFR
          xnu = one/re_cent 
    
          !
          fcpg = two/yly * (re/re_cent)**2
    
          ! Calculate the related bulk Reynolds number (Pope, "Turbulent Flows")
          re_bulk = (re/0.09_mytype)**(1.0_mytype/0.88_mytype)
      
      else
      
          ! Calculate Re_tau in the case of CFR (re = Re_0 in this case)
          re_tau = 0.116_mytype*(re**0.88_mytype)
          
          ! Calculate the related bulk Reynolds number (Pope, "Turbulent Flows")
          re_bulk = (re_tau/0.09_mytype)**(1.0_mytype/0.88_mytype)
          
      end if
  
  end if

  if (ilmn) then
     if (ivarcoeff) then
        npress = 2 !! Need current pressure and previous iterate
     else
        npress = 1
     endif
  endif

  ! 2D snapshot is not compatible with coarse visualization
  if (output2D.ne.0) nvisu = 1
#ifdef ADIOS2
  if (nvisu .ne. 1) then
     if (nrank .eq. 0) then
        print *, "ADIOS2 output is not compatible with coarse visualisation"
        print *, "disabling coarse visualisation"
        print *, "To compress the IO, see ADIOS2 options"
     endif
     nvisu = 1
  endif
#if defined(DOUBLE_PREC) && defined(SAVE_SINGLE)
  print *, "ADIOS2 does not support mixing the simulation and output precision"
  call MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
#endif
#endif

  if (iimplicit.ne.0) then
     
     !if ((itimescheme==5).or.(itimescheme==6)) then
     !   if (nrank==0) write(*,*) 'Error: implicit Y diffusion not yet compatible with RK time schemes'
     !   stop
     !endif
     
     if (isecondder==5) then
        if (nrank==0) write(*,*)  "Warning : support for implicit Y diffusion and "
        if (nrank==0) write(*,*)  "isecondder=5 is experimental. "
      endif
     if (iimplicit==1) then
        xcst = dt * xnu
     else if (iimplicit==2) then
        xcst = dt * xnu * half
     else
        if (nrank==0) write(*,*)  'Error: wrong value for iimplicit ', iimplicit
        stop
     endif
     if (iscalar.eq.1) xcst_sc = xcst / sc
  endif
  
  ! Creating /data folder
  if (nrank==0) call execute_command_line('mkdir -p data')
  
  ! Creating /restart_info folder inside /data
  if (nrank==0) call execute_command_line('mkdir -p data/restart_info')
  
  ! Creating /monitoring folder inside /data
  if (nrank==0) call execute_command_line('mkdir -p data/monitoring')
  
  ! Creating /planes folder inside /data
  if (nrank==0) call execute_command_line('mkdir -p data/planes')
  
  ! Creating /checkpoints folder inside /data
  if (nrank==0) call execute_command_line('mkdir -p data/checkpoints')
  
#ifdef DEBG
  if (nrank == 0) write(*,*) '# parameter input.i3d done'
#endif
  if (nrank==0) then
     print *,'==========================================================='
     if (itype.eq.itype_channel) then
        print *,'Simulating channel'
     elseif (itype.eq.itype_ttbl) then
        print *,'Temporal TBL'
     else
        print *,'Unknown itype: ', itype
        stop
     endif
     print *,'==========================================================='
     if (itype.eq.itype_channel) then
       if (.not.cpg) then
         write(*,*) 'Channel forcing with constant flow rate (CFR)'
         write(*,"(' Re_0 (centerline)             : ',F17.3)") re
         write(*,"(' Re_B (bulk, estimated)        : ',F17.3)") re_bulk
         write(*,"(' Re_tau (estimated)            : ',F17.3)") re_tau
       else 
         write(*,*) 'Channel forcing with constant pressure gradient (CPG)'
         write(*,"(' Re_tau                        : ',F17.3)") re
         write(*,"(' Re_0 (centerline, estimated)  : ',F17.3)") re_cent
         write(*,"(' Re_B (bulk, estimated)        : ',F17.3)") re_bulk
         write(*,"(' fcpg                          : ',F17.8)") fcpg
       end if
     else
       write(*,"(' Reynolds number Re              : ',F17.3)") re
     endif
     
     write(*,"(' xnu                           : ',F17.8)") xnu
     
     ! Displaying if we are using wall oscillations and if open or closed loop strategy
     if (iswitch_wo .eq. 1) then
         write(*,*) '==========================================================='
         write(*,*) 'Spanwise wall oscillations enabled'
     
         if (ifeedback_control .eq. 0) then
             write(*,*) 'Open-loop control'
         else if (ifeedback_control .eq. 1) then
             write(*,*) 'Closed-loop control (feedback control)'
         end if    
     end if
     
     write(*,*) '==========================================================='
     write(*,"(' p_row, p_col           : ',I9, I8)") p_row, p_col
     write(*,*) '==========================================================='
     if(icfllim .eq. 0) then
     write(*,"(' Time step dt           : ',F17.8)") dt
     else if(icfllim .eq. 1) then
     write(*,"(' CFL max                : ',F17.8)") cfl_limit     
     end if
     
     ! Time schemes
     if (itimescheme.eq.1) then
       write(*,"(' Temporal scheme        : ',A20)") "Forward Euler"
     elseif (itimescheme.eq.2) then
       write(*,"(' Temporal scheme        : ',A20)") "Adams-Bashforth 2"
     elseif (itimescheme.eq.3) then
       write(*,"(' Temporal scheme        : ',A20)") "Adams-Bashforth 3"
     elseif (itimescheme.eq.4) then
       write(*,"(' Temporal scheme        : ',A20)") "Adams-Bashforth 4"
       print *,'Error: Adams-bashforth 4 not implemented!'
       stop
     elseif (itimescheme.eq.5) then
       write(*,"(' Temporal scheme        : ',A20)") "Runge-Kutta 3"
     elseif (itimescheme.eq.6) then
       write(*,"(' Temporal scheme        : ',A20)") "Runge-Kutta 4"
       print *,'Error: Runge-kutta 4 not implemented!'
       stop
     else
       print *,'Error: itimescheme must be specified as 1-6'
       stop
     endif
     
     ! Semi-implicit y-diffusion     
     if (iimplicit.ne.0) then
       if (iimplicit.eq.1) then
         write(*,"('            ',A40)") "With backward Euler for Y diffusion"
       else if (iimplicit.eq.2) then
         write(*,"('            ',A40)") "With CN for Y diffusion"
       endif
     endif
     
     ! Displaying the specific model adopted
     write(*,*) '==========================================================='   
     if (ilesmod==0) then
          write(*,"(' Turbulence closure     : ',A17)") "DNS"
     else
       if (jles==1) then
          write(*,"(' Turbulence closure     : ',A17)") "Phys Smag"
       else if (jles==2) then
          write(*,"(' Turbulence closure     : ',A17)") "Phys WALE"      
       else if (jles==3) then
          write(*,"(' Turbulence closure     : ',A17)") "Phys dyn. Smag"
       else if (jles==4) then
          write(*,"(' Turbulence closure     : ',A17)") "iSVV"
       endif
     endif
     
     write(*,*) '==========================================================='
     write(*,"(' ifirst                 : ',I17)") ifirst
     write(*,"(' ilast                  : ',I17)") ilast
     write(*,*) '==========================================================='
     write(*,"(' Lx                     : ',F17.8)") xlx
     write(*,"(' Ly                     : ',F17.8)") yly
     write(*,"(' Lz                     : ',F17.8)") zlz
     write(*,"(' nx                     : ',I17)") nx
     write(*,"(' ny                     : ',I17)") ny
     write(*,"(' nz                     : ',I17)") nz
     write(*,*) '==========================================================='
     write(*,"(' istret                 : ',I17)") istret
     write(*,"(' beta                   : ',F17.8)") beta
     write(*,*) '==========================================================='
     write(*,"(' nu0nu                  : ',F17.8)") nu0nu
     write(*,"(' cnu                    : ',F17.8)") cnu
     write(*,*) '==========================================================='
     if (iscalar==0) write(*,"(' Scalar                 : ',A17)") "off"
     if (iscalar==1) write(*,"(' Scalar                 : ',A17)") "on"
     write(*,"(' numscalar              : ',I17)") numscalar
     if (iscalar.eq.1) then
       do is=1, numscalar
          write(*,"(' Schmidt number sc(',I2,')  : ',F17.8)") is, sc(is)
          write(*,"(' Richardson n.  ri(',I2,')  : ',F17.8)") is, ri(is)
          if (scalar_lbound(is).gt.-huge(one)) then
             write(*,"(' Lower bound      (',I2,')  : ',F17.8)") is, scalar_lbound(is)
          else
             ! This is the default option, no information printed in the listing
          endif
          if (scalar_ubound(is).lt.huge(one)) then
             write(*,"(' Upper bound      (',I2,')  : ',F17.8)") is, scalar_ubound(is)
          else
             ! This is the default option, no information printed in the listing
          endif
          if (iscalar.eq.1) then
             if (nclxS1.eq.1 .or. nclxSn.eq.1 .or. &
                 nclyS1.eq.1 .or. nclySn.eq.1 .or. &
                 nclzS1.eq.1 .or. nclzSn.eq.1) then
                if (sc_even(is)) then
                   ! This is the default option, no information printed in the listing
                else
                   write(*,"(' Scalar ',I2,' is odd')") is
                endif
             endif
             if (sc_skew(is)) then
                write(*,"(' Scalar ',I2,' with skew-symmetric convection')") is
             else
                ! This is the default option, no information printed in the listing
             endif
          endif
       end do
     endif
     write(*,*) '==========================================================='
     write(*,"(' spinup_time            : ',I17)") spinup_time
     write(*,"(' wrotation              : ',F17.8)") wrotation
     write(*,*) '==========================================================='
     if (iibm==0) write(*,"(' Immersed boundary      : ',A17)") "off"
     if (iibm.gt.1) then
      write(*,"(' Immersed boundary      : ',A17)") "on"
      write(*,"(' iibm                   : ',I17)") iibm
     end if
     if (iibm==1) write(*,*) 'Simple immersed boundary method'
     if (iibm==2) then
       write(*,*) 'Lagrangian polynomial reconstruction'
       write(*,*) '==========================================================='
       write(*,"(' npif                   : ',I17)") npif
       write(*,"(' izap                   : ',I17)") izap
       write(*,"(' nraf                   : ',I17)") nraf
       write(*,"(' nobjmax                : ',I17)") nobjmax
     end if
     write(*,*) '==========================================================='
     write(*,"(' Boundary condition velocity field: ')")
     write(*,"(' nclx1, nclxn           : ',I15,',',I1 )") nclx1,nclxn
     write(*,"(' ncly1, nclyn           : ',I15,',',I1 )") ncly1,nclyn
     write(*,"(' nclz1, nclzn           : ',I15,',',I1 )") nclz1,nclzn
     write(*,*) '==========================================================='
     if ((iscalar==1).or.(ilmn)) then
       write(*,"(' Boundary condition scalar field: ')")
       write(*,"(' nclxS1, nclxSn         : ',I15,',',I1 )") nclxS1,nclxSn
       write(*,"(' nclyS1, nclySn         : ',I15,',',I1 )") nclyS1,nclySn
       write(*,"(' nclzS1, nclzSn         : ',I15,',',I1 )") nclzS1,nclzSn
       write(*,*) '==========================================================='
     endif

#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
     write(*,*) 'Numerical precision: Double, saving in single'
#else
     print *,'Numerical precision: Double'
#endif
#else
     write(*,*) 'Numerical precision: Single'
#endif
     write(*,*) '==========================================================='
     write(*,"(' Gravity vector     : (gx, gy, gz)=(',F15.8,',',F15.8,',',F15.8,')')") gravx, gravy, gravz
     if (ilmn) then
        write(*,*)  "LMN                : Enabled"
        if (ivarcoeff) then
           write(*,*)  "LMN-Poisson solver : Variable-coefficient"
        else
           write(*,*)  "LMN-Poisson solver : Constant-coefficient"
        endif
        if (ilmn_bound) then
           write(*,*)  "LMN boundedness    : Enforced"
        else
           write(*,*)  "LMN boundedness    : Not enforced"
        endif
        write(*,"(' dens1 and dens2    : ',F6.2,',',F6.2)") dens1, dens2
        write(*,"(' Prandtl number Re  : ',F15.8)") prandtl
     endif
     write(*,*) ' '
     write(*,*) '==========================================================='
  endif
  
  if (iibm.eq.3) then ! This is only for the Cubic Spline Reconstruction
     npif=npif+1
  endif

#ifdef DEBG
  if (nrank == 0) write(*,*) '# parameter done'
#endif

  return
end subroutine parameter

!-----------------------------------------------------------------------------!
!  SUBROUTINE: parameter_defaults                                          
! DESCRIPTION: Sets the default simulation parameters.                     
!      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>        
!-----------------------------------------------------------------------------!
subroutine parameter_defaults()

  use param
  use variables
  use decomp_2d
  use complex_geometry

  use visu, only : output2D
  
  implicit none

  integer :: i
  
  ! BasicParam
  istret      = 0
  beta        = 0
  irestart    = 0
  iin         = 0
  init_noise  = zero
  iibm        = 0                      
  ilmn        = .FALSE.                      
  ilesmod     = 0
  iscalar     = 0
  ipost       = 0
  ifilter     = 0
  C_filter    = 0.49_mytype
  gravx       = zero
  gravy       = zero
  gravz       = zero 
  
  ! NumOptions
  ipinter     = 3
  itimescheme = 3
  iimplicit   = 0                 
  nu0nu       = four
  cnu         = 0.44_mytype
  
  ! InOutParam
  icheckpoint = 40000   ! Frequency for writing backup files
  ioutput = 40000       ! Frequency for saving snapshots
  ioutput_cf = 200      ! Frequency for saving cf and related quantities
  ioutput_plane = 200   ! Frequency for saving planes for visualization
  ilist = 25            ! Frequency for writing to screen (out/log file)
  isnap = 1             ! Save snapshots (0: no, 1: yes)
  ivisu = 1             ! Save case-specific field for visualization (e.g. Q-criterion) (0: no, 1: yes)        
  nvisu = 1             ! Size for visualisation collection (2: every 2 mesh nodes, 4: every 4 mesh nodes)
  output2D = 0          ! Writing snapshots on a plane (0: no, 1: x-dir, 2: y-dir, 3: z-dir)
  start_output = 1      ! Time-step at which we start to save snapshots (valid for both 3d and 2d snapshots)                                             

  ! AdditionalControls
  iswitch_wo  = 0  ! wall oscillations (0: no, 1: yes)
  
  ! WallOscillations 
  ifeedback_control = 0  ! Switcher to enable feedback control from run-time streamwise shear velocity (closed loop) (0: no, 1: yes)
  a_wo = twelve          ! Amplitude of spanwise wall oscillations (in friction units if feedback control enabled) 
  t_wo = onehundred      ! Period of spanwise wall oscillations (in friction units if feedback control enabled) 
  in_phase = zero        ! Initial phase of the wall oscillations, given as fraction of pi [rad]

  ! ChannelParam 
  cpg         = .FALSE.
  idir_stream = 1
  wrotation   = zero

  ! TemporalTBLParam
  uwall       = one          ! Velocity of translating bottom wall (U_wall)   
  twd         = one          ! Trip wire diameter (D)
  uln         = 0.01_mytype  ! Upper limit of the noise; (uwall - um) < uln*uwall; (default value as Kozul et al. (2016))
  lln         = 0.5_mytype   ! Lower limit of the noise; y+ restriction, based on the mean gradient of the IC
  phiwall     = one          ! Scalar value at the wall
   
  ! ScalarParam
  ! ...
  
  ! LESModel 
  smagwalldamp = 0
  nsmag        = 1
  
  ! IbmStuff
  nobjmax = 0
  nraf    = 0
  npif    = 2
  izap    = 1
  
  ! LMN
  dens1           = one
  dens2           = one
  prandtl         = one
  ilmn_bound      = .TRUE.
  ivarcoeff       = .FALSE.
  ilmn_solve_temp = .FALSE.
  imultispecies   = .FALSE.
  primary_species = -1
  Fr              = zero    
  ibirman_eos     = .FALSE. 
  
  ! Statistics               
  initstat = huge(i)                
  
  ! ExtraNumControl 
  icfllim   = 0     ! Switcher to enable CFL limit constraint (0: no, 1: yes)
  cfl_limit = 0.95  ! CFL limit to adjust time-step 
  
  !-- Additional parameters not present in namelists --!                  
  imodulo2  = 1
  filepath  = './data/'
  datapath  = './data/'
  itime0    = 0
  t0        = zero
  pressure0 = one
  irotation = 0
  itest     = 1
  npress    = 1  ! By default only one pressure field is needed

end subroutine parameter_defaults




