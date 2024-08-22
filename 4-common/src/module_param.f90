!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!-----------------------------------------------------------------------------!
! This module stores (mainly) work variables.
!-----------------------------------------------------------------------------!
module variables
    
  use decomp_2d, only : mytype

  ! Boundary conditions : ncl = 2 --> Dirichlet
  ! Boundary conditions : ncl = 1 --> Free-slip
  ! Boundary conditions : ncl = 0 --> Periodic
  ! l: power of 2,3,4,5 and 6
  ! if ncl = 1 or 2, --> n  = 2l+ 1
  !                  --> nm = n - 1
  !                  --> m  = n + 1
  ! If ncl = 0,      --> n  = 2*l
  !                  --> nm = n
  !                  --> m  = n + 2
  
  ! nstat = size arrays for statistic collection
  ! 2-->every 2 mesh nodes
  ! 4-->every 4 mesh nodes
  ! nvisu = size arrays for visualization collection

  ! Possible n points: 3 5 7 9 11 13 17 19 21 25 31 33 37 41 49 51 55 61 65 73 81 91 97 101 109 121 129 145 151 161 163 181 193 201 217 241 251 257 271 289 301 321 325 361 385 401 433 451 481 487 501 513 541 577 601 641 649 721 751 769 801 811 865 901 961 973 1001 1025 1081 1153 1201 1251 1281 1297 1351 1441 1459 1501 1537 1601 1621 1729 1801 1921 1945 2001 2049 2161 2251 2305 2401 2431 2501 2561 2593 2701 2881 2917 3001 3073 3201 3241 3457 3601 3751 3841 3889 4001 4051 4097 4321 4375 4501 4609 4801 4861 5001 5121 5185 5401 5761 5833 6001 6145 6251 6401 6481 6751 6913 7201 7291 7501 7681 7777 8001 8101 8193 8641 8749 9001 9217 9601 9721 enough

  integer :: nx,ny,nz,numscalar,p_row,p_col,nxm,nym,nzm
  integer :: nstat=1,nprobe=1,nlength=1  

  real(mytype),allocatable,dimension(:) :: sc,uset,cp,ri,group
  
#ifndef DOUBLE_PREC
  integer,parameter :: prec = 4
#else
#ifdef SAVE_SINGLE
  integer,parameter :: prec = 4
#else
  integer,parameter :: prec = 8
#endif
#endif

  ! module filter
  real(mytype),dimension(200) :: idata
  real(mytype),allocatable,dimension(:) :: fiffx, fifcx, fifbx, fisfx, fiscx, fisbx,fifsx,fifwx,fissx,fiswx
  real(mytype),allocatable,dimension(:) :: fiffxp,fifsxp,fifwxp,fisfxp,fissxp,fiswxp
  real(mytype),allocatable,dimension(:) :: fiffy, fifcy, fifby, fisfy, fiscy, fisby,fifsy,fifwy,fissy,fiswy
  real(mytype),allocatable,dimension(:) :: fiffyp,fifsyp,fifwyp,fisfyp,fissyp,fiswyp
  real(mytype),allocatable,dimension(:) :: fiffz, fifcz, fifbz, fisfz, fiscz, fisbz,fifsz,fifwz,fissz,fiswz
  real(mytype),allocatable,dimension(:) :: fiffzp,fifszp,fifwzp,fisfzp,fisszp,fiswzp

  real(mytype),allocatable,dimension(:,:) :: fisx,fivx
  real(mytype),allocatable,dimension(:,:) :: fisy,fivy
  real(mytype),allocatable,dimension(:,:) :: fisz,fivz

  ! module derivative
  real(mytype),allocatable,dimension(:) :: ffx,sfx,fsx,fwx,ssx,swx
  real(mytype),allocatable,dimension(:) :: ffxp,sfxp,fsxp,fwxp,ssxp,swxp
  real(mytype),allocatable,dimension(:) :: ffy,sfy,fsy,fwy,ssy,swy
  real(mytype),allocatable,dimension(:) :: ffyp,sfyp,fsyp,fwyp,ssyp,swyp
  real(mytype),allocatable,dimension(:) :: ffz,sfz,fsz,fwz,ssz,swz
  real(mytype),allocatable,dimension(:) :: ffzp,sfzp,fszp,fwzp,sszp,swzp

  real(mytype),allocatable,dimension(:) :: ffxS,sfxS,fsxS,fwxS,ssxS,swxS
  real(mytype),allocatable,dimension(:) :: ffxpS,sfxpS,fsxpS,fwxpS,ssxpS,swxpS
  real(mytype),allocatable,dimension(:) :: ffyS,sfyS,fsyS,fwyS,ssyS,swyS
  real(mytype),allocatable,dimension(:) :: ffypS,sfypS,fsypS,fwypS,ssypS,swypS
  real(mytype),allocatable,dimension(:) :: ffzS,sfzS,fszS,fwzS,sszS,swzS
  real(mytype),allocatable,dimension(:) :: ffzpS,sfzpS,fszpS,fwzpS,sszpS,swzpS

  real(mytype), save, allocatable, dimension(:,:) :: sx,vx
  real(mytype), save, allocatable, dimension(:,:) :: sy,vy
  real(mytype), save, allocatable, dimension(:,:) :: sz,vz

  ! module implicit
  real(mytype), dimension(:), pointer :: gg, hh, ss, rr, vv, ww, zz, lo1, lo2, lo3, up1, up2, up3
  ! velocity, ncly1 = 2, nclyn = 2
  real(mytype), allocatable, target, dimension(:) :: aam,bbm,ccm,ddm,eem,ggm,hhm,wwm,zzm
  real(mytype), allocatable, target, dimension(:) :: rrm,qqm,vvm,ssm
  real(mytype), allocatable, target, dimension(:) :: sssm, zzzm, ttm, uum  !!Nona
  ! velocity, ncly1 = 1, nclyn = 1, npaire = 0
  real(mytype), allocatable, target, dimension(:) :: aam10,bbm10,ccm10,ddm10,eem10,ggm10,hhm10,wwm10,zzm10
  real(mytype), allocatable, target, dimension(:) :: rrm10,qqm10,vvm10,ssm10
  ! velocity, ncly1 = 1, nclyn = 1, npaire = 1
  real(mytype), allocatable, target, dimension(:) :: aam11,bbm11,ccm11,ddm11,eem11,ggm11,hhm11,wwm11,zzm11
  real(mytype), allocatable, target, dimension(:) :: rrm11,qqm11,vvm11,ssm11
  ! velocity, ncly1 = 0, nclyn = 0
  real(mytype), allocatable, target, dimension(:) :: aam0,bbm0,ccm0,ddm0,eem0,ggm0,hhm0,wwm0,zzm0
  real(mytype), allocatable, target, dimension(:) :: rrm0,qqm0,vvm0,ssm0,l1m,l2m,l3m,u1m,u2m,u3m
  ! velocity, ncly1 = 1, nclyn = 2, npaire = 0
  real(mytype), allocatable, target, dimension(:) :: aam120,bbm120,ccm120,ddm120,eem120,ggm120,hhm120,wwm120,zzm120
  real(mytype), allocatable, target, dimension(:) :: rrm120,qqm120,vvm120,ssm120
  ! velocity, ncly1 = 1, nclyn = 2, npaire = 1
  real(mytype), allocatable, target, dimension(:) :: aam121,bbm121,ccm121,ddm121,eem121,ggm121,hhm121,wwm121,zzm121
  real(mytype), allocatable, target, dimension(:) :: rrm121,qqm121,vvm121,ssm121
  ! velocity, ncly1 = 2, nclyn = 1, npaire = 0
  real(mytype), allocatable, target, dimension(:) :: aam210,bbm210,ccm210,ddm210,eem210,ggm210,hhm210,wwm210,zzm210
  real(mytype), allocatable, target, dimension(:) :: rrm210,qqm210,vvm210,ssm210
  ! velocity, ncly1 = 2, nclyn = 1, npaire = 1
  real(mytype), allocatable, target, dimension(:) :: aam211,bbm211,ccm211,ddm211,eem211,ggm211,hhm211,wwm211,zzm211
  real(mytype), allocatable, target, dimension(:) :: rrm211,qqm211,vvm211,ssm211
  ! scalar, ncly1 = 2, nclyn = 2
  real(mytype), allocatable, target, dimension(:,:) :: aamt,bbmt,ccmt,ddmt,eemt,ggmt,hhmt,wwmt,zzmt
  real(mytype), allocatable, target, dimension(:,:) :: rrmt,qqmt,vvmt,ssmt
  real(mytype), allocatable, target, dimension(:,:) :: uumt,ttmt,sssmt,zzzmt !! Nona
  ! scalar, ncly1 = 1, nclyn = 1, npaire = 0
  real(mytype), allocatable, target, dimension(:,:) :: aam10t,bbm10t,ccm10t,ddm10t,eem10t,ggm10t,hhm10t,wwm10t,zzm10t
  real(mytype), allocatable, target, dimension(:,:) :: rrm10t,qqm10t,vvm10t,ssm10t
  ! scalar, ncly1 = 1, nclyn = 1, npaire = 1
  real(mytype), allocatable, target, dimension(:,:) :: aam11t,bbm11t,ccm11t,ddm11t,eem11t,ggm11t,hhm11t,wwm11t,zzm11t
  real(mytype), allocatable, target, dimension(:,:) :: rrm11t,qqm11t,vvm11t,ssm11t
  ! scalar, ncly1 = 0, nclyn = 0
  real(mytype), allocatable, target, dimension(:,:) :: aam0t,bbm0t,ccm0t,ddm0t,eem0t,ggm0t,hhm0t,wwm0t,zzm0t
  real(mytype), allocatable, target, dimension(:,:) :: rrm0t,qqm0t,vvm0t,ssm0t,l1mt,l2mt,l3mt,u1mt,u2mt,u3mt
  ! scalar, ncly1 = 1, nclyn = 2, npaire = 0
  real(mytype), allocatable, target, dimension(:,:) :: aam120t,bbm120t,ccm120t,ddm120t,eem120t,ggm120t,hhm120t,wwm120t,zzm120t
  real(mytype), allocatable, target, dimension(:,:) :: rrm120t,qqm120t,vvm120t,ssm120t
  ! scalar, ncly1 = 1, nclyn = 2, npaire = 1
  real(mytype), allocatable, target, dimension(:,:) :: aam121t,bbm121t,ccm121t,ddm121t,eem121t,ggm121t,hhm121t,wwm121t,zzm121t
  real(mytype), allocatable, target, dimension(:,:) :: rrm121t,qqm121t,vvm121t,ssm121t
  ! scalar, ncly1 = 2, nclyn = 1, npaire = 0
  real(mytype), allocatable, target, dimension(:,:) :: aam210t,bbm210t,ccm210t,ddm210t,eem210t,ggm210t,hhm210t,wwm210t,zzm210t
  real(mytype), allocatable, target, dimension(:,:) :: rrm210t,qqm210t,vvm210t,ssm210t
  ! scalar, ncly1 = 2, nclyn = 1, npaire = 1
  real(mytype), allocatable, target, dimension(:,:) :: aam211t,bbm211t,ccm211t,ddm211t,eem211t,ggm211t,hhm211t,wwm211t,zzm211t
  real(mytype), allocatable, target, dimension(:,:) :: rrm211t,qqm211t,vvm211t,ssm211t

  ! Interface for derivatives subroutines
  ABSTRACT INTERFACE
     SUBROUTINE DERIVATIVE_X(t,u,r,s,ff,fs,fw,nx,ny,nz,npaire,lind)
       use decomp_2d, only : mytype
       integer :: nx,ny,nz,npaire
       real(mytype), dimension(nx,ny,nz) :: t,u,r
       real(mytype), dimension(ny,nz):: s
       real(mytype), dimension(nx):: ff,fs,fw
       real(mytype) :: lind
     END SUBROUTINE DERIVATIVE_X
     SUBROUTINE DERIVATIVE_Y(t,u,r,s,ff,fs,fw,pp,nx,ny,nz,npaire,lind)
       use decomp_2d, only : mytype
       integer :: nx,ny,nz,npaire
       real(mytype), dimension(nx,ny,nz) :: t,u,r
       real(mytype), dimension(nx,nz):: s
       real(mytype), dimension(ny):: ff,fs,fw,pp
       real(mytype) :: lind
     END SUBROUTINE DERIVATIVE_Y
     SUBROUTINE DERIVATIVE_YY(t,u,r,s,ff,fs,fw,nx,ny,nz,npaire,lind)
       use decomp_2d, only : mytype
       integer :: nx,ny,nz,npaire
       real(mytype), dimension(nx,ny,nz) :: t,u,r
       real(mytype), dimension(nx,nz):: s
       real(mytype), dimension(ny):: ff,fs,fw
       real(mytype) :: lind
     END SUBROUTINE DERIVATIVE_YY
     SUBROUTINE DERIVATIVE_Z(t,u,r,s,ff,fs,fw,nx,ny,nz,npaire,lind)
       use decomp_2d, only : mytype
       integer :: nx,ny,nz,npaire
       real(mytype), dimension(nx,ny,nz) :: t,u,r
       real(mytype), dimension(nx,ny):: s
       real(mytype), dimension(nz):: ff,fs,fw
       real(mytype) :: lind
     END SUBROUTINE DERIVATIVE_Z
     
     ! Additional subroutines for 1-dimensional derivatives in y-direction (R. Corsini)
     SUBROUTINE DERIVATIVE_Y_1D(t,u,r,s,ff,fs,fw,pp,ny,npaire)
       use decomp_2d, only : mytype
       integer :: ny,npaire
       real(mytype), dimension(ny) :: t,u,r
       real(mytype) :: s
       real(mytype), dimension(ny):: ff,fs,fw,pp
     END SUBROUTINE DERIVATIVE_Y_1D
     SUBROUTINE DERIVATIVE_YY_1D(t,u,r,s,ff,fs,fw,ny,npaire)
       use decomp_2d, only : mytype
       integer :: ny,npaire
       real(mytype), dimension(ny) :: t,u,r
       real(mytype) :: s
       real(mytype), dimension(ny):: ff,fs,fw
     END SUBROUTINE DERIVATIVE_YY_1D
     
  END INTERFACE

  PROCEDURE (DERIVATIVE_X) derx_00,derx_11,derx_12,derx_21,derx_22,&
       derxx_00,derxx_11,derxx_12,derxx_21,derxx_22
  PROCEDURE (DERIVATIVE_X), POINTER :: derx,derxx,derxS,derxxS
  PROCEDURE (DERIVATIVE_Y) dery_00,dery_11,dery_12,dery_21,dery_22
  PROCEDURE (DERIVATIVE_Y), POINTER :: dery,deryS
  PROCEDURE (DERIVATIVE_YY) &
       deryy_00,deryy_11,deryy_12,deryy_21,deryy_22
  PROCEDURE (DERIVATIVE_YY), POINTER :: deryy,deryyS
  PROCEDURE (DERIVATIVE_Z) derz_00,derz_11,derz_12,derz_21,derz_22,&
       derzz_00,derzz_11,derzz_12,derzz_21,derzz_22
  PROCEDURE (DERIVATIVE_Z), POINTER :: derz,derzz,derzS,derzzS
  
  ! Additional subroutines for 1-dimensional derivatives in y-direction (R. Corsini)
  PROCEDURE (DERIVATIVE_Y_1D) &
       dery1D_00,dery1D_11,dery1D_12,dery1D_21,dery1D_22
  PROCEDURE (DERIVATIVE_Y_1D), POINTER :: dery1D,deryS1D
  PROCEDURE (DERIVATIVE_YY_1D) &
       deryy1D_00,deryy1D_11,deryy1D_12,deryy1D_21,deryy1D_22
  PROCEDURE (DERIVATIVE_YY_1D), POINTER :: deryy1D,deryyS1D

  ! O6SVV (Order 6th Spectral Vanishing Viscosity)
  real(mytype),allocatable,dimension(:) :: newsm,newtm,newsmt,newtmt
  real(mytype),allocatable,dimension(:) :: newrm,newrmt

  ABSTRACT INTERFACE
     SUBROUTINE FILTER_X(t,u,r,s,ff,fs,fw,nx,ny,nz,npaire,lind)
       use decomp_2d, only : mytype
       integer :: nx,ny,nz,npaire
       real(mytype), dimension(nx,ny,nz) :: t,u,r
       real(mytype), dimension(ny,nz):: s
       real(mytype), dimension(nx):: ff,fs,fw
       real(mytype) :: lind
     END SUBROUTINE FILTER_X
     SUBROUTINE FILTER_Y(t,u,r,s,ff,fs,fw,nx,ny,nz,npaire,lind)
       use decomp_2d, only : mytype
       integer :: nx,ny,nz,npaire
       real(mytype), dimension(nx,ny,nz) :: t,u,r
       real(mytype), dimension(nx,nz):: s
       real(mytype), dimension(ny):: ff,fs,fw
       real(mytype) :: lind
     END SUBROUTINE FILTER_Y
     SUBROUTINE FILTER_Z(t,u,r,s,ff,fs,fw,nx,ny,nz,npaire,lind)
       use decomp_2d, only : mytype
       integer :: nx,ny,nz,npaire
       real(mytype), dimension(nx,ny,nz) :: t,u,r
       real(mytype), dimension(nx,ny):: s
       real(mytype), dimension(nz):: ff,fs,fw
       real(mytype) :: lind
     END SUBROUTINE FILTER_Z
  END INTERFACE

  PROCEDURE (FILTER_X) filx_00,filx_11, filx_12, filx_21, filx_22
  PROCEDURE (FILTER_X), POINTER :: filx,filxS
  PROCEDURE (FILTER_Y) fily_00,fily_11, fily_12, fily_21, fily_22
  PROCEDURE (FILTER_Y), POINTER :: fily,filyS
  PROCEDURE (FILTER_Z) filz_00,filz_11, filz_12, filz_21, filz_22
  PROCEDURE (FILTER_Z), POINTER :: filz,filzS

  ! module pressure
  real(mytype), save, allocatable, dimension(:,:) :: dpdyx1,dpdyxn,dpdzx1,dpdzxn
  real(mytype), save, allocatable, dimension(:,:) :: dpdxy1,dpdxyn,dpdzy1,dpdzyn
  real(mytype), save, allocatable, dimension(:,:) :: dpdxz1,dpdxzn,dpdyz1,dpdyzn

  ! module inflow (BCs)
  real(mytype), save, allocatable, dimension(:,:) :: bxx1,bxy1,bxz1,bxxn,bxyn,bxzn,bxo,byo,bzo
  real(mytype), save, allocatable, dimension(:,:) :: byx1,byy1,byz1,byxn,byyn,byzn
  real(mytype), save, allocatable, dimension(:,:) :: bzx1,bzy1,bzz1,bzxn,bzyn,bzzn

  ! module derpres
  real(mytype),allocatable,dimension(:) :: cfx6,ccx6,cbx6,cfxp6,ciwxp6,csxp6,&
       cwxp6,csx6,cwx6,cifx6,cicx6,cisx6
  real(mytype),allocatable,dimension(:) :: cibx6,cifxp6,cisxp6,ciwx6
  real(mytype),allocatable,dimension(:) :: cfi6,cci6,cbi6,cfip6,csip6,cwip6,csi6,&
       cwi6,cifi6,cici6,cibi6,cifip6
  real(mytype),allocatable,dimension(:) :: cisip6,ciwip6,cisi6,ciwi6
  real(mytype),allocatable,dimension(:) :: cfy6,ccy6,cby6,cfyp6,csyp6,cwyp6,csy6
  real(mytype),allocatable,dimension(:) :: cwy6,cify6,cicy6,ciby6,cifyp6,cisyp6,&
       ciwyp6,cisy6,ciwy6
  real(mytype),allocatable,dimension(:) :: cfi6y,cci6y,cbi6y,cfip6y,csip6y,cwip6y,&
       csi6y,cwi6y,cifi6y,cici6y
  real(mytype),allocatable,dimension(:) :: cibi6y,cifip6y,cisip6y,ciwip6y,cisi6y,ciwi6y
  real(mytype),allocatable,dimension(:) :: cfz6,ccz6,cbz6,cfzp6,cszp6,cwzp6,csz6
  real(mytype),allocatable,dimension(:) :: cwz6,cifz6,cicz6,cibz6,cifzp6,ciszp6,&
       ciwzp6,cisz6,ciwz6
  real(mytype),allocatable,dimension(:) :: cfi6z,cci6z,cbi6z,cfip6z,csip6z,cwip6z,&
       csi6z,cwi6z,cifi6z,cici6z
  real(mytype),allocatable,dimension(:) :: cibi6z,cifip6z,cisip6z,ciwip6z,cisi6z,ciwi6z

  ! module waves (spectral Poisson solver)
  complex(mytype),allocatable,dimension(:) :: zkz,zk2,ezs
  complex(mytype),allocatable,dimension(:) :: yky,yk2,eys
  complex(mytype),allocatable,dimension(:) :: xkx,xk2,exs

  ! module mesh
  real(mytype),allocatable,dimension(:) :: ppy,pp2y,pp4y
  real(mytype),allocatable,dimension(:) :: ppyi,pp2yi,pp4yi
  real(mytype),allocatable,dimension(:) :: xp,xpi,yp,ypi,dyp,zp,zpi,del
  real(mytype),allocatable,dimension(:) :: yeta,yetai
  real(mytype) :: alpha,beta

end module variables
!-----------------------------------------------------------------------------!
! This module stores (mainly) variables that do not change runtime
! (they are parameters, even if they are not in Fortran meaning).
!-----------------------------------------------------------------------------!
module param

  use decomp_2d, only : mytype

  ! BCs (velocity and scalar fields)
  integer :: nclx1,nclxn,ncly1,nclyn,nclz1,nclzn
  integer :: nclxS1,nclxSn,nclyS1,nclySn,nclzS1,nclzSn

  ! Logical variable for boundary condition that is true in periodic case
  ! and false otherwise
  logical :: nclx,ncly,nclz

  ! Flow cases
  integer :: itype
  integer, parameter :: itype_channel = 3,  &
                        itype_ttbl    = 13

  ! Time step, first and last
  integer :: itime,ifirst,ilast
  
  ! How many (sub)timestpeps do we need to store?
  integer :: ntime 
  integer :: itime0
  integer :: iscalar,iadvance_time
  
  real(mytype) :: xlx,yly,zlz,dx,dy,dz,dx2,dy2,dz2,t,t0
  real(mytype) :: dt,re,xnu,init_noise
  
  ! Basic Param
  integer :: irestart

  ! Mesh and domain
  integer :: istret
  
  ! Type of initialization (it depends on the specific flowcase)
  integer :: iin
  
  ! Iterations of sub-time steps for RK schemes
  integer :: itr
  
  ! Low-Mach Number
  real(mytype) :: dens1, dens2
  
  ! Compact filter
  real(mytype) :: C_filter
  integer      :: ifilter
  
  ! Logical, true when synchronization is needed
  logical, save :: sync_vel_needed = .true.
  logical, save :: sync_scal_needed = .true.
  
  ! Numerics control
  integer      :: ifirstder,isecondder,ipinter,itimescheme,iimplicit
  real(mytype) :: nu0nu, cnu

  ! In/Out Param
  integer :: icheckpoint    ! Frequency for writing backup files
  integer :: ioutput        ! Frequency for saving snapshots 
  integer :: ioutput_cf     ! Frequency for saving cf and related quantities
  integer :: ioutput_plane  ! Frequency for saving planes for visualization
  integer :: ilist          ! Frequency for writing to screen (out/log file)
  integer :: isnap          ! Save snapshots (0: no, 1: yes)
  integer :: ivisu          ! Save case-specific field for visualization (e.g. Q-criterion) (0: no, 1: yes)
  integer :: nvisu          ! Size for visualisation collection (2: every 2 mesh nodes, 4: every 4 mesh nodes)
  integer :: start_output   ! Time-step at which we start to save snapshots (valid for both 3d and 2d snapshots)
  
  !--- Channel flow ---!
  logical            :: cpg
  integer            :: idir_stream, spinup_time
  real(mytype)       :: wrotation
  real(mytype)       :: re_cent, fcpg, re_bulk, re_tau
  
  ! Bulk velocity
  real(mytype), save :: ubulk
  !--------------------!

  ! Numerical Fourier (CFL diffusion parameter)
  real(mytype) :: cfl_diff_x,cfl_diff_y,cfl_diff_z,cfl_diff_sum

  ! Coefficient for the implicit y-diffusion terms (it involves viscosity and time-step)
  real(mytype)                            :: xcst     ! For velocity field
  real(mytype), allocatable, dimension(:) :: xcst_sc  ! For scalar fields
  
  real(mytype), allocatable, dimension(:,:) :: alpha_sc, beta_sc, g_sc

  !! Scalars
  logical, allocatable, dimension(:) :: sc_even, sc_skew
  real(mytype), allocatable, dimension(:) :: scalar_lbound, scalar_ubound
  real(mytype) :: Tref
  
  ! IBM
  integer :: npif,izap,ianal,iibm

  ! LES
  integer :: itest           ! frequency to print to screen LES constants 
  integer :: iles            ! to enable LES modelling (0: no, 1: yes)
  integer :: ilesmod, iwall  ! modelling options for LES
  integer :: jles            ! LES model
  integer :: smagwalldamp
  real(mytype) :: smagcst,nSmag,walecst,FSGS,pr_t,maxdsmagcst

  !! Gravity field (vector components)
  real(mytype) :: gravx, gravy, gravz

  !! LMN
  logical :: ilmn, ilmn_bound, ilmn_solve_temp
  real(mytype) :: pressure0, prandtl, Fr
  integer :: nrhotime, npress
  logical :: ivarcoeff

  logical :: imultispecies
  logical, allocatable, dimension(:) :: massfrac
  real(mytype), allocatable, dimension(:) :: mol_weight
  integer :: primary_species

  logical :: ibirman_eos
  
  ! Coefficients for time-integration schemes
  real(mytype), dimension(5) :: adt,bdt,cdt,ddt,gdt
    
  ! Shear quantities at the wall (used for Channel and TTBL)
  real(mytype), save :: sh_vel        ! Total shear velocity
  real(mytype), save :: sh_velx       ! Shear velocity along x
  real(mytype), save :: sh_velz       ! Shear velocity along z
  real(mytype), save :: fric_coeff    ! Total skin friction coefficient
  real(mytype), save :: fric_coeffx   ! Skin friction coefficient along x
  real(mytype), save :: fric_coeffz   ! Skin friction coefficient along z
  real(mytype), save :: t_viscous     ! Viscous time unit (based on shear velocity and viscous length)
  
  real(mytype), save :: deltaxplus    ! Delta x^+
  real(mytype), save :: deltayplusw   ! Delta y^+ at the wall
  real(mytype), save :: deltazplus    ! Delta z^+
  real(mytype), save :: deltayplusd   ! Delta y^+ at the TTBL mean interface or at channel half-height
    
  real(mytype), save :: xlxplus       ! Lx^+
  real(mytype), save :: ylyplus       ! Ly^+
  real(mytype), save :: zlzplus       ! Lz^+
  real(mytype), save :: delta_nu      ! Viscous unit (length)
    
  !--- Additional controls namelist and related quantities ---!
  
  ! Spanwise wall oscillations
  integer      :: iswitch_wo          ! Switcher to enable the reading of wall-oscillation parameters and the application of wall oscillations
  real(mytype) :: span_vel            ! Spanwise velocity at the wall due to imposed wall oscillations, variable calculated
  
  !--- TTBL namelist and related quantities ---!
  
  ! Temporal TBL input parameters
  real(mytype) :: uwall               ! Velocity of translating bottom wall (U_wall) 
  real(mytype) :: twd                 ! Trip wire diameter (D)
  real(mytype) :: uln                 ! Upper limit of the noise; (uwall - um) < uln*uwall; (default value as Kozul et al.)
  real(mytype) :: lln                 ! Lower limit of the noise; y+ restriction, based on the mean gradient of the IC 
  real(mytype) :: phiwall             ! Scalar value at the wall
  
  ! Temporal TBL output(s)
  real(mytype) :: powerin             ! Power input to the TTBL 
  
  ! Boundary layer thickness parameters (TTBLs only)
  real(mytype), save :: delta_99      ! BL thickness (1% of velocity of the wall)   
  real(mytype), save :: re_tau_tbl    ! Friction Reynolds number based on mean streamwise gradient
  integer,      save :: counter       ! Index to save the position of the yp array at which there is the BL mean interface 
  real(mytype), save :: mean_phigwtot ! Mean scalar gradient at the wall total (all processors)
  real(mytype), save :: A_fact        ! Analogy factor to monitor runtime the degree of dissimilarity
  
  !--- Extra numerics control namelist and related quantities ---!
  
  ! Extra controls for numerics 
  integer      :: icfllim             ! Index or switcher for enabling CFL limit constraint (0: no, 1: yes)
  real(mytype) :: cfl_limit           ! CFL limit to adjust the time-step
  
  ! Variable to save and show the maximum CFL in the restart file
  real(mytype) :: cflmax
  
  !--- Wall oscillations namelist and related quantities ---!
  
  ! Parameters for wall oscillations (used for channel flows and TTBLs)
  real(mytype) :: a_wo                ! Amplitude of spanwise wall oscillations (in friction units if feedback control enabled)  
  real(mytype) :: t_wo                ! Period of spanwise wall oscillations (in friction units if feedback control enabled) 
  integer      :: ifeedback_control   ! Switcher to enable feedback control from run-time streamwise shear velocity (closed loop)
  real(mytype) :: in_phase            ! Initial phase of the wall oscillations, given as fraction of pi [rad]
    
  ! Numbers
  real(mytype),parameter :: zpzeroone=0.01_mytype
  real(mytype),parameter :: zpone=0.1_mytype
  real(mytype),parameter :: zptwo=0.2_mytype
  real(mytype),parameter :: zptwoone=0.21_mytype
  real(mytype),parameter :: zptwofive=0.25_mytype
  real(mytype),parameter :: zpthree=0.3_mytype
  real(mytype),parameter :: zpfour=0.4_mytype
  real(mytype),parameter :: zpfive=0.5_mytype
  real(mytype),parameter :: zpsix=0.6_mytype
  real(mytype),parameter :: zpseven=0.7_mytype
  real(mytype),parameter :: zpeight=0.8_mytype
  real(mytype),parameter :: zpnine=0.9_mytype
  real(mytype),parameter :: zpninefive=0.95_mytype

  real(mytype),parameter :: zero=0.0_mytype
  real(mytype),parameter :: half=0.5_mytype
  real(mytype),parameter :: twothird=2.0_mytype/3.0_mytype
  real(mytype),parameter :: one=1.0_mytype
  real(mytype),parameter :: onepfive=1.5_mytype
  real(mytype),parameter :: two=2.0_mytype
  real(mytype),parameter :: twopfive=2.5_mytype
  real(mytype),parameter :: three=3.0_mytype
  real(mytype),parameter :: threepfive=3.5_mytype
  real(mytype),parameter :: four=4.0_mytype
  real(mytype),parameter :: five=5.0_mytype
  real(mytype),parameter :: six=6.0_mytype
  real(mytype),parameter :: seven=7.0_mytype
  real(mytype),parameter :: eight=8.0_mytype
  real(mytype),parameter :: nine=9.0_mytype

  real(mytype),parameter :: ten=10.0_mytype
  real(mytype),parameter :: eleven=11.0_mytype
  real(mytype),parameter :: twelve=12.0_mytype
  real(mytype),parameter :: thirteen=13.0_mytype
  real(mytype),parameter :: fourteen=14.0_mytype
  real(mytype),parameter :: fifteen=15.0_mytype
  real(mytype),parameter :: sixteen=16.0_mytype
  real(mytype),parameter :: seventeen=17.0_mytype
  real(mytype),parameter :: eighteen=18.0_mytype

  real(mytype),parameter :: twenty=20.0_mytype
  real(mytype),parameter :: twentyone=21.0_mytype
  real(mytype),parameter :: twentythree=23.0_mytype
  real(mytype),parameter :: twentyfour=24.0_mytype
  real(mytype),parameter :: twentyfive=25.0_mytype
  real(mytype),parameter :: twentyseven=27.0_mytype
  real(mytype),parameter :: twentyeight=28.0_mytype
  
  real(mytype),parameter :: thirty=30.0_mytype
  real(mytype),parameter :: thirtytwo=32.0_mytype
  real(mytype),parameter :: thirtyfour=34.0_mytype
  real(mytype),parameter :: thirtysix=36.0_mytype
  real(mytype),parameter :: thirtyseven=37.0_mytype
  
  real(mytype),parameter :: forty=40.0_mytype
  real(mytype),parameter :: fortyfour=44.0_mytype
  real(mytype),parameter :: fortyfive=45.0_mytype
  real(mytype),parameter :: fortyeight=48.0_mytype
  
  real(mytype),parameter :: fifty=50.0_mytype
  real(mytype),parameter :: fiftyfour=54.0_mytype
  real(mytype),parameter :: fiftyfive=55.0_mytype
  real(mytype),parameter :: fiftynine=59.0_mytype
  
  real(mytype),parameter :: sixty=60.0_mytype
  real(mytype),parameter :: sixtytwo=62.0_mytype
  real(mytype),parameter :: sixtythree=63.0_mytype
  
  real(mytype),parameter :: seventy=70.0_mytype
  real(mytype),parameter :: seventyfive=75.0_mytype
  
  real(mytype),parameter :: onehundred=100.0_mytype
  real(mytype),parameter :: onehundredtwentysix=126.0_mytype
  real(mytype),parameter :: onehundredtwentyeight=128.0_mytype
  real(mytype),parameter :: onehundredeighty=180.0_mytype
  
  real(mytype),parameter :: twohundredsix=206.0_mytype
  real(mytype),parameter :: twohundredeight=208.0_mytype
  real(mytype),parameter :: twohundredfiftysix=256.0_mytype
  real(mytype),parameter :: twohundredseventytwo=272.0_mytype
  
  real(mytype),parameter :: onethousand=1000.0_mytype
  real(mytype),parameter :: twothousand=2000.0_mytype
  real(mytype),parameter :: threethousandsixhundred=3600.0_mytype
  
  complex(mytype),parameter :: cx_one_one=cmplx(one, one, kind=mytype)

  real(mytype),parameter :: pi=acos(-one)
  real(mytype),parameter :: twopi=two*acos(-one)

end module param
!-----------------------------------------------------------------------------!
module complex_geometry

  use decomp_2d,only : mytype
  use variables,only : nx,ny,nz,nxm,nym,nzm

  integer     ,allocatable,dimension(:,:)   :: nobjx,nobjy,nobjz
  integer     ,allocatable,dimension(:,:,:) :: nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif
  real(mytype),allocatable,dimension(:,:,:) :: xi,xf,yi,yf,zi,zf
  integer :: nxraf,nyraf,nzraf,nraf,nobjmax
end module complex_geometry
!-----------------------------------------------------------------------------!
module derivX

  use decomp_2d, only : mytype

  real(mytype) :: alcaix6,acix6,bcix6
  real(mytype) :: ailcaix6,aicix6,bicix6,cicix6,dicix6
  real(mytype) :: alfa1x,af1x,bf1x,cf1x,df1x,alfa2x,af2x,alfanx,afnx,bfnx
  real(mytype) :: cfnx,dfnx,alfamx,afmx,alfaix,afix,bfix,alsa1x,as1x,bs1x
  real(mytype) :: cs1x,ds1x,alsa2x,as2x,alsanx,asnx,bsnx,csnx,dsnx,alsamx
  real(mytype) :: asmx,alsa3x,as3x,bs3x,alsatx,astx,bstx

  ! O6SVV (Order 6th Spectral Vanishing Viscosity)
  real(mytype) :: alsa4x,as4x,bs4x,cs4x
  real(mytype) :: alsattx,asttx,bsttx,csttx
  real(mytype) :: alsaix,asix,bsix,csix,dsix
  real(mytype) :: alsaixt,asixt,bsixt,csixt,dsixt

end module derivX
!-----------------------------------------------------------------------------!
module derivY

  use decomp_2d, only : mytype

  real(mytype) :: alcaiy6,aciy6,bciy6
  real(mytype) :: ailcaiy6,aiciy6,biciy6,ciciy6,diciy6
  real(mytype) :: alfa1y,af1y,bf1y,cf1y,df1y,alfa2y,af2y,alfany,afny,bfny
  real(mytype) :: cfny,dfny,alfamy,afmy,alfajy,afjy,bfjy,alsa1y,as1y,bs1y
  real(mytype) :: cs1y,ds1y,alsa2y,as2y,alsany,asny,bsny,csny,dsny,alsamy
  real(mytype) :: asmy,alsa3y,as3y,bs3y,alsaty,asty,bsty

  ! O6SVV (Order 6th Spectral Vanishing Viscosity)
  real(mytype) :: alsa4y,as4y,bs4y,cs4y
  real(mytype) :: alsatty,astty,bstty,cstty
  real(mytype) :: alsajy,asjy,bsjy,csjy,dsjy
  real(mytype) :: alsajyt,asjyt,bsjyt,csjyt,dsjyt

end module derivY
!-----------------------------------------------------------------------------!
module derivZ

  use decomp_2d, only : mytype

  real(mytype) :: alcaiz6,aciz6,bciz6
  real(mytype) :: ailcaiz6,aiciz6,biciz6,ciciz6,diciz6
  real(mytype) :: alfa1z,af1z,bf1z,cf1z,df1z,alfa2z,af2z,alfanz,afnz,bfnz
  real(mytype) :: cfnz,dfnz,alfamz,afmz,alfakz,afkz,bfkz,alsa1z,as1z,bs1z
  real(mytype) :: cs1z,ds1z,alsa2z,as2z,alsanz,asnz,bsnz,csnz,dsnz,alsamz
  real(mytype) :: asmz,alsa3z,as3z,bs3z,alsatz,astz,bstz

  ! O6SVV (Order 6th Spectral Vanishing Viscosity)
  real(mytype) :: alsa4z,as4z,bs4z,cs4z
  real(mytype) :: alsattz,asttz,bsttz,csttz
  real(mytype) :: alsakz,askz,bskz,cskz,dskz
  real(mytype) :: alsakzt,askzt,bskzt,cskzt,dskzt

end module derivZ
!-----------------------------------------------------------------------------!
! Describes the parameters for the discrete filters in x-pencils
module parfiX
  use decomp_2d, only : mytype
  real(mytype) :: fial1x, fia1x, fib1x, fic1x, fid1x, fie1x, fif1x  ! Coefficients for filter at boundary point 1
  real(mytype) :: fial2x, fia2x, fib2x, fic2x, fid2x, fie2x, fif2x  ! Coefficients for filter at boundary point 2
  real(mytype) :: fial3x, fia3x, fib3x, fic3x, fid3x, fie3x, fif3x  ! Coefficients for filter at boundary point 3
  real(mytype) :: fialix, fiaix, fibix, ficix, fidix                ! Coefficient for filter at interior points
  real(mytype) :: fialnx, fianx, fibnx, ficnx, fidnx, fienx, fifnx  ! Coefficient for filter at boundary point n
  real(mytype) :: fialmx, fiamx, fibmx, ficmx, fidmx, fiemx, fifmx  ! Coefficient for filter at boundary point m=n-1
  real(mytype) :: fialpx, fiapx, fibpx, ficpx, fidpx, fiepx, fifpx  ! Coefficient for filter at boundary point p=n-2
end module parfiX
!-----------------------------------------------------------------------------!
! Describes the parameters for the discrete filters in y-pencils
module parfiY
  use decomp_2d, only : mytype
  real(mytype) :: fial1y, fia1y, fib1y, fic1y, fid1y, fie1y, fif1y ! Coefficients for filter at boundary point 1
  real(mytype) :: fial2y, fia2y, fib2y, fic2y, fid2y, fie2y, fif2y ! Coefficients for filter at boundary point 2
  real(mytype) :: fial3y, fia3y, fib3y, fic3y, fid3y, fie3y, fif3y ! Coefficients for filter at boundary point 3
  real(mytype) :: fialjy, fiajy, fibjy, ficjy, fidjy               ! Coefficient for filter at interior points
  real(mytype) :: fialny, fiany, fibny, ficny, fidny, fieny, fifny ! Coefficient for filter at boundary point n
  real(mytype) :: fialmy, fiamy, fibmy, ficmy, fidmy, fiemy, fifmy ! Coefficient for filter at boundary point m=n-1
  real(mytype) :: fialpy, fiapy, fibpy, ficpy, fidpy, fiepy, fifpy ! Coefficient for filter at boundary point p=n-2
end module parfiY
!-----------------------------------------------------------------------------!
! Describes the parameters for the discrete filters in z-pencils
module parfiZ
  use decomp_2d, only : mytype
  real(mytype) :: fial1z, fia1z, fib1z, fic1z, fid1z, fie1z, fif1z ! Coefficients for filter at boundary point 1
  real(mytype) :: fial2z, fia2z, fib2z, fic2z, fid2z, fie2z, fif2z ! Coefficients for filter at boundary point 2
  real(mytype) :: fial3z, fia3z, fib3z, fic3z, fid3z, fie3z, fif3z ! Coefficients for filter at boundary point 3
  real(mytype) :: fialkz, fiakz, fibkz, fickz, fidkz               ! Coefficient for filter at interior points
  real(mytype) :: fialnz, fianz, fibnz, ficnz, fidnz, fienz, fifnz ! Coefficient for filter at boundary point n
  real(mytype) :: fialmz, fiamz, fibmz, ficmz, fidmz, fiemz, fifmz ! Coefficient for filter at boundary point m=n-1
  real(mytype) :: fialpz, fiapz, fibpz, ficpz, fidpz, fiepz, fifpz ! Coefficient for filter at boundary point p=n-2
end module parfiZ
!-----------------------------------------------------------------------------!
module simulation_stats
  real(8) :: tstart,time1,trank,tranksum,ttotal,tremaining,telapsed
end module simulation_stats
!-----------------------------------------------------------------------------!
module ibm_param
  use decomp_2d, only : mytype
  real(mytype) :: cex,cey,cez,ra,ubcx,ubcy,ubcz,rads, c_air
  real(mytype) :: chord,thickness,omega
  integer :: inana ! Analytical BC as Input
  integer :: imove
end module ibm_param
!-----------------------------------------------------------------------------!
