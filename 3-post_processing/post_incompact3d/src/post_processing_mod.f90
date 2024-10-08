
!This file is not part of standard Xcompact3d releases (xcompact3d.com).

!-----------------------------------------------------------------------------!
! DESCRIPTION: This module is used to store post-processing variables used in 
!              'post_incompact3d' and basic subroutines to allocate memory 
!              and to reset to zero arrays. 
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
!              Roberto Corsini <roberto.corsini@unimore.it> 
!-----------------------------------------------------------------------------!

module post_processing

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  use variables
  use param
  use var

  implicit none
  
  ! Variables declaration
  !integer :: FS
  !character(len=100) :: fileformat
  !character(len=1), parameter :: NL=char(10) !new line character

  ! Logicals for if conditions during post-processing
  logical, save :: post_mean,post_vort,post_diss,post_corz,post_tke_eq
  logical, save :: read_vel,read_pre,read_phi 

  !--- Arrays for statistic collection ---!

  ! Explanation of nomenclature:
  ! - mean:   simple arithmetic mean of data at the same (x,y,z) location across different
  !           flow realizations and / or time units;
  ! - meanH1: sum over homogeneous directions, on a single processor ('1') ('H': Homogeneous);
  ! - meanHT: sum over different processors (through MPI) (T: Total).
  
  ! Point-value  
  real(mytype), save, allocatable, dimension(:,:,:) :: u1mean,v1mean,w1mean
  real(mytype), save, allocatable, dimension(:,:,:) :: u2mean,v2mean,w2mean
  real(mytype), save, allocatable, dimension(:,:,:) :: u3mean,v3mean,w3mean
  real(mytype), save, allocatable, dimension(:,:,:) :: u4mean,v4mean,w4mean
  real(mytype), save, allocatable, dimension(:,:,:) :: uvmean,uwmean,vwmean
  real(mytype), save, allocatable, dimension(:,:,:) :: pre1mean,pre2mean
  real(mytype), save, allocatable, dimension(:,:,:) :: vpremean
  
  ! Sum at each processor
  real(mytype), save, allocatable, dimension(:) :: u1meanH1,v1meanH1,w1meanH1
  real(mytype), save, allocatable, dimension(:) :: u2meanH1,v2meanH1,w2meanH1
  real(mytype), save, allocatable, dimension(:) :: u3meanH1,v3meanH1,w3meanH1
  real(mytype), save, allocatable, dimension(:) :: u4meanH1,v4meanH1,w4meanH1
  real(mytype), save, allocatable, dimension(:) :: uvmeanH1,uwmeanH1,vwmeanH1
  real(mytype), save, allocatable, dimension(:) :: pre1meanH1,pre2meanH1
  real(mytype), save, allocatable, dimension(:) :: vpremeanH1
  
  ! Total sum
  real(mytype), save, allocatable, dimension(:) :: u1meanHT,v1meanHT,w1meanHT
  real(mytype), save, allocatable, dimension(:) :: u2meanHT,v2meanHT,w2meanHT
  real(mytype), save, allocatable, dimension(:) :: u3meanHT,v3meanHT,w3meanHT
  real(mytype), save, allocatable, dimension(:) :: u4meanHT,v4meanHT,w4meanHT
  real(mytype), save, allocatable, dimension(:) :: uvmeanHT,uwmeanHT,vwmeanHT
  real(mytype), save, allocatable, dimension(:) :: pre1meanHT,pre2meanHT
  real(mytype), save, allocatable, dimension(:) :: vpremeanHT
  
  !--- Scalar field (only 1 supported at the moment) ---!
  
  ! Point-value     
  real(mytype), save, allocatable, dimension(:,:,:) :: phi1mean,phi2mean
  real(mytype), save, allocatable, dimension(:,:,:) :: uphimean,vphimean,wphimean

  ! Sum at each processor   
  real(mytype), save, allocatable, dimension(:) :: phi1meanH1,phi2meanH1
  real(mytype), save, allocatable, dimension(:) :: uphimeanH1,vphimeanH1,wphimeanH1

  ! Total sum  
  real(mytype), save, allocatable, dimension(:) :: phi1meanHT,phi2meanHT
  real(mytype), save, allocatable, dimension(:) :: uphimeanHT,vphimeanHT,wphimeanHT
            
  !--- Arrays for vorticity ---!
  real(mytype), save, allocatable, dimension(:,:,:) :: vortxmean,vortymean,vortzmean
  real(mytype), save, allocatable, dimension(:)     :: vortxmeanH1,vortymeanH1,vortzmeanH1
  real(mytype), save, allocatable, dimension(:)     :: vortxmeanHT,vortymeanHT,vortzmeanHT
    
  ! Arrays for mean streamwise gradient (dU/dy) (x: streamwise)
  real(mytype), save, allocatable, dimension(:,:,:) :: mean_gradientx
  real(mytype), save, allocatable, dimension(:)     :: mean_gradientxH1
  real(mytype), save, allocatable, dimension(:)     :: mean_gradientxHT
  
  ! Arrays for mean spanwise gradient (dW/dy) (z: spanwise)
  real(mytype), save, allocatable, dimension(:,:,:) :: mean_gradientz
  real(mytype), save, allocatable, dimension(:)     :: mean_gradientzH1
  real(mytype), save, allocatable, dimension(:)     :: mean_gradientzHT
  
  ! Arrays for scalar gradient (dPhi/dy)
  real(mytype), save, allocatable, dimension(:,:,:) :: mean_gradphi
  real(mytype), save, allocatable, dimension(:)     :: mean_gradphiH1
  real(mytype), save, allocatable, dimension(:)     :: mean_gradphiHT
  
  ! Arrays for total dissipation rate eps
  real(mytype), save, allocatable, dimension(:,:,:) :: epsmean
  real(mytype), save, allocatable, dimension(:)     :: epsmeanH1
  real(mytype), save, allocatable, dimension(:)     :: epsmeanHT
  
  ! Arrays for correlation functions (velocity and scalar (p: phi)) 
  real(mytype), save, allocatable, dimension(:,:)   :: RuuzH1,RvvzH1,RwwzH1,RuvzH1,RsszH1
  real(mytype), save, allocatable, dimension(:,:)   :: RuuzHT,RvvzHT,RwwzHT,RuvzHT,RsszHT
  
  !--- Arrays for fluctuating terms for TKE equation ---!
  real(mytype), save, allocatable, dimension(:,:,:) :: kvprime_mean,  pseudo_eps_tke_mean
  real(mytype), save, allocatable, dimension(:)     :: kvprime_meanH1,pseudo_eps_tke_meanH1
  real(mytype), save, allocatable, dimension(:)     :: kvprime_meanHT,pseudo_eps_tke_meanHT
  
  !--- Arrays for TKE equation ---!
  real(mytype), save, allocatable, dimension(:)     :: tke_diffHT      ! diffusive transport term: -nu * d^2 (<k>) / dy^2
  real(mytype), save, allocatable, dimension(:)     :: tke_prodHT      ! production term: - <u'v'> dU/dy
  
  ! Work arrays for TKE
  real(mytype), save, allocatable, dimension(:)     :: temp_dery       ! temporary variable to store the derivative in y of a generic quantity 
  real(mytype), save, allocatable, dimension(:)     :: di1d
  real(mytype) :: sy1d
     
contains

  !-----------------------------------------------------------------------------!
  ! DESCRIPTION: Subroutine to allocate memory for working arrays used only
  !              during post-processing. 
  !   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
  !-----------------------------------------------------------------------------!
  subroutine init_post_variables
  
    use var
    
    implicit none
    
    ! Print that we are initializing post-processing work variables
    if (nrank == 0) write(*,*) '==========================================================='
    if (nrank == 0) write(*,*) 'Initializing post-processing work variables ...'
       
    ! Allocate y-pencil pressure array (not allocated in the solver)
    call alloc_y(pre2)
    pre2=zero      

  end subroutine init_post_variables
  
  !-----------------------------------------------------------------------------!
  ! DESCRIPTION: Subroutine to allocate memory for statistics arrays.
  !   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
  !              Roberto Corsini <roberto.corsini@unimore.it> 
  !-----------------------------------------------------------------------------!
  subroutine init_statistics()
    
    implicit none
  
    ! Print that we are initializing post-processing statistics variables
    if (nrank == 0) write(*,*) '==========================================================='
    if (nrank == 0) write(*,*) 'Initializing post-processing statistics variables ...'
      
    if (post_mean) then
       
       ! Point-value
       allocate(u1mean  (ysize(1),ysize(2),ysize(3))); u1mean   = zero    
       allocate(v1mean  (ysize(1),ysize(2),ysize(3))); v1mean   = zero   
       allocate(w1mean  (ysize(1),ysize(2),ysize(3))); w1mean   = zero  
       allocate(u2mean  (ysize(1),ysize(2),ysize(3))); u2mean   = zero   
       allocate(v2mean  (ysize(1),ysize(2),ysize(3))); v2mean   = zero    
       allocate(w2mean  (ysize(1),ysize(2),ysize(3))); w2mean   = zero   
       allocate(u3mean  (ysize(1),ysize(2),ysize(3))); u3mean   = zero    
       allocate(v3mean  (ysize(1),ysize(2),ysize(3))); v3mean   = zero   
       allocate(w3mean  (ysize(1),ysize(2),ysize(3))); w3mean   = zero    
       allocate(u4mean  (ysize(1),ysize(2),ysize(3))); u4mean   = zero 
       allocate(v4mean  (ysize(1),ysize(2),ysize(3))); v4mean   = zero   
       allocate(w4mean  (ysize(1),ysize(2),ysize(3))); w4mean   = zero      
       allocate(uvmean  (ysize(1),ysize(2),ysize(3))); uvmean   = zero      
       allocate(uwmean  (ysize(1),ysize(2),ysize(3))); uwmean   = zero       
       allocate(vwmean  (ysize(1),ysize(2),ysize(3))); vwmean   = zero    
       allocate(pre1mean(ysize(1),ysize(2),ysize(3))); pre1mean = zero   
       allocate(pre2mean(ysize(1),ysize(2),ysize(3))); pre2mean = zero
       allocate(vpremean(ysize(1),ysize(2),ysize(3))); vpremean = zero 
       
       ! Sum at each processor
       allocate(u1meanH1  (ysize(2))); u1meanH1   = zero  
       allocate(v1meanH1  (ysize(2))); v1meanH1   = zero 
       allocate(w1meanH1  (ysize(2))); w1meanH1   = zero
       allocate(u2meanH1  (ysize(2))); u2meanH1   = zero
       allocate(v2meanH1  (ysize(2))); v2meanH1   = zero
       allocate(w2meanH1  (ysize(2))); w2meanH1   = zero
       allocate(u3meanH1  (ysize(2))); u3meanH1   = zero  
       allocate(v3meanH1  (ysize(2))); v3meanH1   = zero 
       allocate(w3meanH1  (ysize(2))); w3meanH1   = zero
       allocate(u4meanH1  (ysize(2))); u4meanH1   = zero  
       allocate(v4meanH1  (ysize(2))); v4meanH1   = zero
       allocate(w4meanH1  (ysize(2))); w4meanH1   = zero 
       allocate(uvmeanH1  (ysize(2))); uvmeanH1   = zero 
       allocate(uwmeanH1  (ysize(2))); uwmeanH1   = zero 
       allocate(vwmeanH1  (ysize(2))); vwmeanH1   = zero
       allocate(pre1meanH1(ysize(2))); pre1meanH1 = zero 
       allocate(pre2meanH1(ysize(2))); pre2meanH1 = zero
       allocate(vpremeanH1(ysize(2))); vpremeanH1 = zero
   
       ! Total sum
       allocate(u1meanHT  (ysize(2))); u1meanHT   = zero 
       allocate(v1meanHT  (ysize(2))); v1meanHT   = zero
       allocate(w1meanHT  (ysize(2))); w1meanHT   = zero
       allocate(u2meanHT  (ysize(2))); u2meanHT   = zero  
       allocate(v2meanHT  (ysize(2))); v2meanHT   = zero 
       allocate(w2meanHT  (ysize(2))); w2meanHT   = zero
       allocate(u3meanHT  (ysize(2))); u3meanHT   = zero  
       allocate(v3meanHT  (ysize(2))); v3meanHT   = zero
       allocate(w3meanHT  (ysize(2))); w3meanHT   = zero
       allocate(u4meanHT  (ysize(2))); u4meanHT   = zero  
       allocate(v4meanHT  (ysize(2))); v4meanHT   = zero 
       allocate(w4meanHT  (ysize(2))); w4meanHT   = zero
       allocate(uvmeanHT  (ysize(2))); uvmeanHT   = zero
       allocate(uwmeanHT  (ysize(2))); uwmeanHT   = zero
       allocate(vwmeanHT  (ysize(2))); vwmeanHT   = zero 
       allocate(pre1meanHT(ysize(2))); pre1meanHT = zero
       allocate(pre2meanHT(ysize(2))); pre2meanHT = zero
       allocate(vpremeanHT(ysize(2))); vpremeanHT = zero

       !--- With the current version, we can deal only with 1 scalar field (modify to a rank 4 array if needed) ---! 
       
       ! Scalar field      
       
       ! Point-value
       allocate(phi1mean(ysize(1),ysize(2),ysize(3))); phi1mean = zero 
       allocate(phi2mean(ysize(1),ysize(2),ysize(3))); phi2mean = zero     
       allocate(uphimean(ysize(1),ysize(2),ysize(3))); uphimean = zero   
       allocate(vphimean(ysize(1),ysize(2),ysize(3))); vphimean = zero   
       allocate(wphimean(ysize(1),ysize(2),ysize(3))); wphimean = zero    

       ! Sum at each processor
       allocate(phi1meanH1(ysize(2))); phi1meanH1 = zero
       allocate(phi2meanH1(ysize(2))); phi2meanH1 = zero
       allocate(uphimeanH1(ysize(2))); uphimeanH1 = zero
       allocate(vphimeanH1(ysize(2))); vphimeanH1 = zero
       allocate(wphimeanH1(ysize(2))); wphimeanH1 = zero
       
       ! Total sum
       allocate(phi1meanHT(ysize(2))); phi1meanHT = zero
       allocate(phi2meanHT(ysize(2))); phi2meanHT = zero
       allocate(uphimeanHT(ysize(2))); uphimeanHT = zero
       allocate(vphimeanHT(ysize(2))); vphimeanHT = zero
       allocate(wphimeanHT(ysize(2))); wphimeanHT = zero 
            
    end if
    
    if (post_vort) then
       
        ! Vorticity
        allocate(vortxmean(ysize(1),ysize(2),ysize(3))); vortxmean = zero
        allocate(vortymean(ysize(1),ysize(2),ysize(3))); vortymean = zero 
        allocate(vortzmean(ysize(1),ysize(2),ysize(3))); vortzmean = zero           
       
        allocate(vortxmeanH1(ysize(2))); vortxmeanH1 = zero
        allocate(vortymeanH1(ysize(2))); vortymeanH1 = zero
        allocate(vortzmeanH1(ysize(2))); vortzmeanH1 = zero
             
        allocate(vortxmeanHT(ysize(2))); vortxmeanHT = zero
        allocate(vortymeanHT(ysize(2))); vortymeanHT = zero
        allocate(vortzmeanHT(ysize(2))); vortzmeanHT = zero
      
        ! Mean gradients (velocity)
        allocate(mean_gradientx(ysize(1),ysize(2),ysize(3))); mean_gradientx = zero
        allocate(mean_gradientz(ysize(1),ysize(2),ysize(3))); mean_gradientz = zero
       
        allocate(mean_gradientxH1(ysize(2))); mean_gradientxH1 = zero
        allocate(mean_gradientzH1(ysize(2))); mean_gradientzH1 = zero
       
        allocate(mean_gradientxHT(ysize(2))); mean_gradientxHT = zero
        allocate(mean_gradientzHT(ysize(2))); mean_gradientzHT = zero
       
        ! Scalar mean gradient
        allocate(mean_gradphi    (ysize(1),ysize(2),ysize(3))); mean_gradphi   = zero
        allocate(mean_gradphiH1  (ysize(2)));                   mean_gradphiH1 = zero
        allocate(mean_gradphiHT  (ysize(2)));                   mean_gradphiHT = zero
                                                                                                                                                                 
    endif
    
    if (post_diss) then
    
        ! Total dissipation rate
        allocate(epsmean  (ysize(1),ysize(2),ysize(3))); epsmean   = zero
        allocate(epsmeanH1(ysize(2)));                   epsmeanH1 = zero
        allocate(epsmeanHT(ysize(2)));                   epsmeanHT = zero
    
    end if
    
    if (post_corz) then
    
        ! Correlation functions in z-direction (velocity and scalar)
        allocate(RuuzH1(zsize(2),zsize(3))); RuuzH1 = zero
        allocate(RvvzH1(zsize(2),zsize(3))); RvvzH1 = zero
        allocate(RwwzH1(zsize(2),zsize(3))); RwwzH1 = zero
        allocate(RuvzH1(zsize(2),zsize(3))); RuvzH1 = zero      
        allocate(RsszH1(zsize(2),zsize(3))); RsszH1 = zero 
               
        allocate(RuuzHT(zsize(2),zsize(3))); RuuzHT = zero
        allocate(RvvzHT(zsize(2),zsize(3))); RvvzHT = zero
        allocate(RwwzHT(zsize(2),zsize(3))); RwwzHT = zero
        allocate(RuvzHT(zsize(2),zsize(3))); RuvzHT = zero
        allocate(RsszHT(zsize(2),zsize(3))); RsszHT = zero
    
    end if
    
    if (post_tke_eq) then
       
        ! Fluctuating terms for TKE equation (turbulent transport and pseudo dissipation)
        allocate(kvprime_mean       (ysize(1),ysize(2),ysize(3))); kvprime_mean        = zero
        allocate(pseudo_eps_tke_mean(ysize(1),ysize(2),ysize(3))); pseudo_eps_tke_mean = zero
        
        allocate(kvprime_meanH1       (ysize(2))); kvprime_meanH1        = zero
        allocate(pseudo_eps_tke_meanH1(ysize(2))); pseudo_eps_tke_meanH1 = zero
        
        allocate(kvprime_meanHT       (ysize(2))); kvprime_meanHT        = zero
        allocate(pseudo_eps_tke_meanHT(ysize(2))); pseudo_eps_tke_meanHT = zero
        
        ! All other terms
        allocate(tke_diffHT   (ysize(2))); tke_diffHT = zero     
        allocate(tke_prodHT   (ysize(2))); tke_prodHT = zero
        allocate(temp_dery    (ysize(2))); temp_dery  = zero
        allocate(di1d         (ysize(2))); di1d       = zero
                
    end if
    
    ! If we need to calculate fluctuations, allocate memory to read mean statistics    
    if (post_corz .or. post_tke_eq) then
    
        ! Velocity statistics       
        allocate(u1meanHT  (ysize(2))); u1meanHT   = zero   
        allocate(v1meanHT  (ysize(2))); v1meanHT   = zero  
        allocate(w1meanHT  (ysize(2))); w1meanHT   = zero
        allocate(u2meanHT  (ysize(2))); u2meanHT   = zero   
        allocate(v2meanHT  (ysize(2))); v2meanHT   = zero  
        allocate(w2meanHT  (ysize(2))); w2meanHT   = zero
        
        ! Reynolds stresses
        allocate(uvmeanHT  (ysize(2))); uvmeanHT   = zero   
        allocate(uwmeanHT  (ysize(2))); uwmeanHT   = zero  
        allocate(vwmeanHT  (ysize(2))); vwmeanHT   = zero
        
        ! Pressure statistics
        allocate(pre1meanHT(ysize(2))); pre1meanHT = zero
        allocate(vpremeanHT(ysize(2))); vpremeanHT = zero
        
        ! Mean scalar field              
        allocate(phi1meanHT(ysize(2))); phi1meanHT = zero
           
    end if

  end subroutine init_statistics
  
  !-----------------------------------------------------------------------------!
  ! DESCRIPTION: Subroutine to reset to zero the arrays of averages on 
  !              each single mesh point. 
  !   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
  !-----------------------------------------------------------------------------!
  subroutine reset_averages()
    
  implicit none
  
  if (post_mean) then
      u1mean   = zero; v1mean   = zero; w1mean   = zero
      u2mean   = zero; v2mean   = zero; w2mean   = zero
      u3mean   = zero; v3mean   = zero; w3mean   = zero
      u4mean   = zero; v4mean   = zero; w4mean   = zero
      
      uvmean   = zero; uwmean   = zero; vwmean   = zero
      pre1mean = zero; pre2mean = zero; vpremean = zero
  
      ! Scalar field
      phi1mean = zero; phi2mean = zero
      uphimean = zero; vphimean = zero; wphimean = zero
      
  end if
  
  if (post_vort) then
      vortxmean      = zero; vortymean      = zero; vortzmean    = zero
      mean_gradientx = zero; mean_gradientz = zero; mean_gradphi = zero
  end if
  
  if (post_diss) then
      epsmean = zero
  end if
  
  if(post_tke_eq) then
      kvprime_mean        = zero
      pseudo_eps_tke_mean = zero
  end if
     
  end subroutine reset_averages

  !-----------------------------------------------------------------------------!
  ! DESCRIPTION: Subroutine to reset to zero the arrays of averages on subdomains
  !              and on total domain. 
  !   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it>
  !-----------------------------------------------------------------------------!  
  subroutine reset_subdomains_and_domain()
    
  implicit none
  
  if (post_mean) then
  
      ! Subdomains   
      u1meanH1   = zero; v1meanH1   = zero; w1meanH1 = zero
      u2meanH1   = zero; v2meanH1   = zero; w2meanH1 = zero
      u3meanH1   = zero; v3meanH1   = zero; w3meanH1 = zero
      u4meanH1   = zero; v4meanH1   = zero; w4meanH1 = zero
      
      uvmeanH1   = zero; uwmeanH1   = zero; vwmeanH1 = zero
      pre1meanH1 = zero; pre2meanH1 = zero
      vpremeanH1 = zero
      
      ! Total domain
      u1meanHT   = zero; v1meanHT   = zero; w1meanHT = zero
      u2meanHT   = zero; v2meanHT   = zero; w2meanHT = zero
      u3meanHT   = zero; v3meanHT   = zero; w3meanHT = zero
      u4meanHT   = zero; v4meanHT   = zero; w4meanHT = zero
      
      uvmeanHT   = zero; uwmeanHT   = zero; vwmeanHT = zero
      pre1meanHT = zero; pre2meanHT = zero
      vpremeanHT = zero
            
      ! Scalar field
            
      ! Subdomains   
      phi1meanH1 = zero; phi2meanH1 = zero
      uphimeanH1 = zero; vphimeanH1 = zero; wphimeanH1 = zero
      
      ! Total domain
      phi1meanHT = zero; phi2meanHT = zero
      uphimeanHT = zero; vphimeanHT = zero; wphimeanHT = zero
           
  end if
  
  if (post_vort) then
  
      ! Subdomains
      vortxmeanH1      = zero; vortymeanH1      = zero; vortzmeanH1    = zero
      mean_gradientxH1 = zero; mean_gradientzH1 = zero; mean_gradphiH1 = zero
  
      ! Total domain
      vortxmeanHT      = zero; vortymeanHT      = zero; vortzmeanHT    = zero
      mean_gradientxHT = zero; mean_gradientzHT = zero; mean_gradphiHT = zero
  
  end if
  
  if (post_diss) then
  
      ! Subdomains
      epsmeanH1 = zero
  
      ! Total domain
      epsmeanHT = zero
  
  end if
  
  if (post_corz) then
  
      ! Subdomains
      RuuzH1 = zero; RvvzH1 = zero; RwwzH1 = zero; RuvzH1 = zero; RsszH1 = zero
  
      ! Total domain
      RuuzHT = zero; RvvzHT = zero; RwwzHT = zero; RuvzHT = zero; RsszHT = zero
  
  end if
  
  if(post_tke_eq) then
  
      ! Subdomains
      kvprime_meanH1 = zero; pseudo_eps_tke_meanH1 = zero
  
      ! Total domain
      kvprime_meanHT = zero; pseudo_eps_tke_meanHT = zero
      
      ! All the other terms are only in total domain
      tke_diffHT = zero     
      tke_prodHT = zero
      temp_dery  = zero
      di1d       = zero
  
  end if
  
  end subroutine reset_subdomains_and_domain
     
end module post_processing
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
! Additional subroutines used only by 'post_incompact3d'
! (they are equivalent to the ones of the solver).
!-----------------------------------------------------------------------------!

  !----------------------------------------------!
  ! Header of the program printed to the screen. !
  !----------------------------------------------!
  subroutine program_header()
  
  use decomp_2d_mpi, only : nrank
    
  implicit none
  
  if (nrank==0) then
     write(*,*) '!---------------------------------------------------------!'
     write(*,*) '!                 ~  PostIncompact3D  ~                   !'
     write(*,*) '!  Copyright (c) 2018 Eric Lamballais and Sylvain Laizet  !'
     write(*,*) '!  Modified by Felipe Schuch and Ricardo Frantz           !'
     write(*,*) '!  Modified by Paul Bartholomew, Georgios Deskos and      !'
     write(*,*) '!  Sylvain Laizet, 2018                                   !'
     write(*,*) '!                                                         !'
     write(*,*) '!  Modified by Filippo Moroni & Roberto Corsini, 2024     !'
     write(*,*) '!---------------------------------------------------------!'
     
#if defined(VERSION)
     write(*,*)'Git version        : ', VERSION
#else
     write(*,*)'Git version        : unknown'
#endif
  endif
  
  end subroutine program_header
  !-----------------------------------------------------------------------------!
  subroutine reading_input_file()

  use decomp_2d_mpi, only : nrank
  use decomp_2d_io
  use variables
  use param
  use var
  use MPI
  
  implicit none
  
  ! Variables to read the input.i3d file
  integer :: nargin, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase
  
  ! Reading of the input file as Xcompact3d does
  nargin=command_argument_count()
  if (nargin <1) then
     InputFN='input.i3d'
     if (nrank==0) write(*,*) 'PostIncompact3d is run with the default file -->', trim(InputFN)
  elseif (nargin >= 1) then
     call get_command_argument(1,InputFN,FNLength,status)
     back=.true.
     FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
     DecInd=index(FNBase,'.',back)
     if (DecInd >1) then
        FNBase=FNBase(1:(DecInd-1))
     end if
     if (nrank==0) write(*,*) 'PostIncompact3d is run with the provided file -->', trim(InputFN)
  endif
  
  ! Reading the input file for geometry and numerics
  call parameter(InputFN)
  
  end subroutine reading_input_file
  !-----------------------------------------------------------------------------!



