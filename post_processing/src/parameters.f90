!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################
!###########################################################################
!
!  SUBROUTINE: parameter
! DESCRIPTION: Reads the input.i3d file and sets the parameters of the
!              simulation.
!      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>
!
!###########################################################################
subroutine parameter(input_i3d)

  use mpi
  
  use iso_fortran_env

  use param
  use variables
  use complex_geometry
  use decomp_2d
  use ibm_param
  use dbg_schemes, only: sin_prec, cos_prec

  use var, only : dphi1

  use visu, only : output2D

  implicit none

  character(len=80), intent(in) :: input_i3d
  real(mytype) :: theta, cfl,cf2
  integer :: longueur ,impi,j, is, total

  NAMELIST /BasicParam/ p_row, p_col, nx, ny, nz, istret, beta, xlx, yly, zlz, &
       itype, iin, re, u1, u2, init_noise, inflow_noise, &
       dt, ifirst, ilast, &
       numscalar, iibm, ilmn, &
       ilesmod, iscalar, &
       nclx1, nclxn, ncly1, nclyn, nclz1, nclzn, &
       ivisu, ipost, &
       gravx, gravy, gravz, &
       cpg, idir_stream, &
       ifilter, C_filter, iturbine
  NAMELIST /NumOptions/ ifirstder, isecondder, itimescheme, iimplicit, &
       nu0nu, cnu, ipinter
  NAMELIST /InOutParam/ irestart, icheckpoint, ioutput, nvisu, ilist, iprocessing, &
       ninflows, ntimesteps, inflowpath, ioutflow, output2D
  NAMELIST /Statistics/ wrotation,spinup_time, nstat, initstat
  NAMELIST /ScalarParam/ sc, ri, uset, cp, &
       nclxS1, nclxSn, nclyS1, nclySn, nclzS1, nclzSn, &
       scalar_lbound, scalar_ubound, sc_even, sc_skew, &
       alpha_sc, beta_sc, g_sc, Tref
  NAMELIST /LESModel/ jles, smagcst, smagwalldamp, nSmag, walecst, maxdsmagcst, iwall
  NAMELIST /WallModel/ smagwalldamp
  NAMELIST /Tripping/ itrip,A_tr,xs_tr_tbl,ys_tr_tbl,ts_tr_tbl,x0_tr_tbl
  NAMELIST /ibmstuff/ cex,cey,cez,ra,nobjmax,nraf, npif, izap, ianal, imove, thickness, chord, omega ,ubcx,ubcy,ubcz,rads, c_air
  NAMELIST /LMN/ dens1, dens2, prandtl, ilmn_bound, ivarcoeff, ilmn_solve_temp, &
       massfrac, mol_weight, imultispecies, primary_species, &
       Fr, ibirman_eos
  NAMELIST /ABL/ z_zero, iwallmodel, k_roughness, ustar, dBL, &
       imassconserve, ibuoyancy, iPressureGradient, iCoriolis, CoriolisFreq, &
       istrat, idamping, iheight, TempRate, TempFlux, itherm, gravv, UG, T_wall, T_top, ishiftedper, iconcprec, pdl 
  NAMELIST /CASE/ tgv_twod
  NAMELIST /ALMParam/ iturboutput,NTurbines,TurbinesPath,NActuatorlines,ActuatorlinesPath,eps_factor,rho_air
  NAMELIST /ADMParam/ Ndiscs,ADMcoords,C_T,aind,iturboutput,rho_air
  
#ifdef DEBG
  if (nrank == 0) write(*,*) '# parameter start'
#endif

  if (nrank==0) then
     write(*,*) '==========================================================='
     write(*,*) '===================== PostIncompact3d ====================='
     write(*,*) '== Copyright (c) 2018 Eric Lamballais and Sylvain Laizet =='
     write(*,*) '== Modified by Felipe Schuch and Ricardo Frantz ==========='
     write(*,*) '== Modified by Paul Bartholomew, Georgios Deskos and ======'
     write(*,*) '== Sylvain Laizet -- 2018- ================================'
     write(*,*) '==========================================================='
     write(*,*) '== Modified by Filippo Moroni -- 2024 ====================='
     write(*,*) '== From original files by Roberto Corsini ================='
     write(*,*) '==========================================================='
     
#if defined(VERSION)
     write(*,*)'Git version        : ', VERSION
#else
     write(*,*)'Git version        : unknown'
#endif
  endif

  call parameter_defaults()

  !! Read parameters
  open(10, file=input_i3d)

  !! These are the 'essential' parameters
  read(10, nml=BasicParam); rewind(10)
  read(10, nml=NumOptions); rewind(10)
  read(10, nml=InOutParam); rewind(10)
  read(10, nml=Statistics); rewind(10)
  
  if (iibm.ne.0) then
     read(10, nml=ibmstuff); rewind(10)
  endif
   
  !! Set Scalar BCs same as fluid (may be overridden) [DEFAULT]
  nclxS1 = nclx1; nclxSn = nclxn
  nclyS1 = ncly1; nclySn = nclyn
  nclzS1 = nclz1; nclzSn = nclzn
  
  if (numscalar.ne.0) then
     iscalar = 1

     !! Allocate scalar arrays and set sensible defaults
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
  
  ! These are the 'optional'/model parameters
  if(ilesmod==0) then
     nu0nu=four
     cnu=0.44_mytype
  endif
  if(ilesmod.ne.0) then
     read(10, nml=LESModel); rewind(10)
  endif
  if (itype.eq.itype_tbl) then
     read(10, nml=Tripping); rewind(10)
  endif
  if (itype.eq.itype_abl) then
     read(10, nml=ABL); rewind(10)
  endif
  
  if (iturbine.eq.1) then
     read(10, nml=ALMParam); rewind(10)
  else if (iturbine.eq.2) then
     read(10, nml=ADMParam); rewind(10)
  endif
  
  !read(10, nml=TurbulenceWallModel); rewind(10)
  
  read(10, nml=CASE); rewind(10)                 ! Read case-specific variables
  
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

  dx=xlx/real(nxm,mytype)
  dy=yly/real(nym,mytype)
  dz=zlz/real(nzm,mytype)

  dx2 = dx * dx
  dy2 = dy * dy
  dz2 = dz * dz

  xnu=one/re
  
  ! Constant pressure gradient, re = Re_tau -> use to compute Re_centerline
  if (cpg) then
    re_cent = (re/0.116_mytype)**(1.0_mytype/0.88_mytype)
    xnu = one/re_cent ! viscosity based on Re_cent to keep same scaling as CFR
    !
    fcpg = two/yly * (re/re_cent)**2
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
     if ((itimescheme==5).or.(itimescheme==6)) then
        if (nrank==0) write(*,*) 'Error: implicit Y diffusion not yet compatible with RK time schemes'
        stop
     endif
     if (isecondder==5) then
        if (nrank==0) write(*,*)  "Warning : support for implicit Y diffusion and isecondder=5 is experimental"
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

  if (itype==itype_tbl.and.A_tr .gt. zero.and.nrank==0)  write(*,*)  "TBL tripping is active"

  anglex = sin_prec(pi*angle/onehundredeighty)
  angley = cos_prec(pi*angle/onehundredeighty)
  !###########################################################################
  ! Log-output
  !###########################################################################
#ifdef DEBG
  if (nrank == 0) write(*,*) '# parameter input.i3d done'
#endif
  if (nrank==0) then
     print *,'==========================================================='
     if (itype.eq.itype_user) then
        print *,'User-defined simulation'
     elseif (itype.eq.itype_lockexch) then
        print *,'Simulating lock-exchange'
     elseif (itype.eq.itype_tgv) then
        print *,'Simulating TGV'
     elseif (itype.eq.itype_channel) then
        print *,'Simulating channel'
     elseif (itype.eq.itype_hill) then
        print *,'Simulating periodic hill'
     elseif (itype.eq.itype_cyl) then
        print *,'Simulating cylinder'
     elseif (itype.eq.itype_dbg) then
        print *,'Debug schemes'
     elseif (itype.eq.itype_mixlayer) then
        print *,'Mixing layer'
     elseif (itype.eq.itype_jet) then
        print *,'Jet'
     elseif (itype.eq.itype_tbl) then
        print *,'Turbulent boundary layer'
     elseif (itype.eq.itype_abl) then
        print *,'Atmospheric boundary layer'
     elseif (itype.eq.itype_uniform) then
        print *,'Uniform flow'
     elseif (itype.eq.itype_sandbox) then
           print *,'Sandbox'
     else
        print *,'Unknown itype: ', itype
        stop
     endif
     print *,'==========================================================='
     if (itype.eq.itype_channel) then
       if (.not.cpg) then
         write(*,*) 'Channel forcing with constant flow rate (CFR)'
         write(*,"(' Re_cl                  : ',F17.3)") re
       else 
         write(*,*) 'Channel forcing with constant pressure gradient (CPG)'
         write(*,"(' Re_tau                 : ',F17.3)") re
         write(*,"(' Re_cl (estimated)      : ',F17.3)") re_cent
         write(*,"(' fcpg                   : ',F17.8)") fcpg
       end if
     else
       write(*,"(' Reynolds number Re     : ',F17.3)") re
     endif
     write(*,"(' xnu                    : ',F17.8)") xnu
     write(*,*) '==========================================================='
     write(*,"(' p_row, p_col           : ',I9, I8)") p_row, p_col
     write(*,*) '==========================================================='
     write(*,"(' Time step dt           : ',F17.8)") dt
     
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
     write(*,"(' High and low speed : u1=',F6.2,' and u2=',F6.2)") u1,u2
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
        write(*,"(' dens1 and dens2    : ',F6.2' ',F6.2)") dens1, dens2
        write(*,"(' Prandtl number Re  : ',F15.8)") prandtl
     endif
     if (angle.ne.0.) write(*,"(' Solid rotation     : ',F6.2)") angle
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

!###########################################################################
!
!  SUBROUTINE: parameter_defaults
! DESCRIPTION: Sets the default simulation parameters.
!      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>
!
!###########################################################################
subroutine parameter_defaults()

  use param
  use variables
  use decomp_2d
  use complex_geometry

  use visu, only : output2D
  
  implicit none

  integer :: i

  ro = 99999999._mytype
  angle = zero
  u1 = 2
  u2 = 1
  init_noise = zero
  inflow_noise = zero
  iin = 0
  itimescheme = 4
  iimplicit = 0
  istret = 0
  ipinter=3
  beta = 0
  ilesmod = 0
  iscalar = 0
  cont_phi = 0
  filepath = './data/'
  irestart = 0
  itime0 = 0
  t0 = zero
  datapath = './data/'

  !! LES stuff
  SmagWallDamp=0
  nSmag=1

  !! IBM stuff
  nraf = 0
  nobjmax = 0

  itrip = 0
  wrotation = zero
  irotation = 0
  itest=1

  !! Gravity field
  gravx = zero
  gravy = zero
  gravz = zero

  !! LMN stuff
  ilmn = .FALSE.
  ilmn_bound = .TRUE.
  pressure0 = one
  prandtl = one
  dens1 = one
  dens2 = one
  ivarcoeff = .FALSE.
  npress = 1 !! By default people only need one pressure field
  ilmn_solve_temp = .FALSE.
  imultispecies = .FALSE.
  Fr = zero
  ibirman_eos = .FALSE.

  primary_species = -1

  !! Channel
  cpg = .false.
  idir_stream = 1

  !! Filter
  ifilter=0
  C_filter=0.49_mytype

  !! ABL
  z_zero=zpone
  k_roughness=zpfour
  ustar=0.45_mytype
  dBL=250._mytype
  iPressureGradient=1
  iwallmodel=1
  imassconserve=0
  ibuoyancy=1
  iheight=0
  itherm=1
  idamping=0
  gravv=9.81_mytype
  TempRate=-zptwofive/3600_mytype
  TempFlux=0.24_mytype
  UG=[zero,zero,zero]
  ishiftedper=0
  iconcprec=0
  pdl=zero
  
  !! Turbine modelling
  iturbine=0
  rho_air=one

  !! IO
  ivisu = 1
  ipost = 0
  iprocessing = huge(i)
  initstat = huge(i)
  ninflows=1
  ntimesteps=1
  inflowpath='./'
  ioutflow=0
  output2D = 0

  ipost = 0
  iibm=0
  npif = 2
  izap = 1

  imodulo2 = 1

  !! CASE specific variables
  tgv_twod = .FALSE.

  !! TRIPPING
  A_tr=zero
  xs_tr_tbl=1.402033_mytype
  ys_tr_tbl=0.350508_mytype
  ts_tr_tbl=1.402033_mytype
  x0_tr_tbl=3.505082_mytype
  
end subroutine parameter_defaults