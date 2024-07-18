!********************************************************************
module post_processing

  use decomp_2d
  use variables
  use param
  use var

  implicit none
  
  ! Variables declaration
  integer :: FS
  character(len=100) :: fileformat
  character(len=1), parameter :: NL=char(10) !new line character

  logical, save :: post_mean,post_vort,post_diss,post_corz,post_tke_eq

  !--- Arrays for statistic collection ---!
  
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
  
  !--- Scalar field ---!
  
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
  
  ! Arrays for mean total parallel gradient (dU_parallel/dy) (p: parallel)
  real(mytype), save, allocatable, dimension(:,:,:) :: mean_gradientp
  real(mytype), save, allocatable, dimension(:)     :: mean_gradientpH1
  real(mytype), save, allocatable, dimension(:)     :: mean_gradientpHT
  
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
  real(mytype), save, allocatable, dimension(:,:)   :: RuuzH1,RvvzH1,RwwzH1,RuvzH1,RppzH1
  real(mytype), save, allocatable, dimension(:,:)   :: RuuzHT,RvvzHT,RwwzHT,RuvzHT,RppzHT
  
  !--- Arrays for fluctuating terms for TKE equation ---!
  real(mytype), save, allocatable, dimension(:,:,:) :: kvprime_mean,  pseudo_eps_tke_mean
  real(mytype), save, allocatable, dimension(:)     :: kvprime_meanH1,pseudo_eps_tke_meanH1
  real(mytype), save, allocatable, dimension(:)     :: kvprime_meanHT,pseudo_eps_tke_meanHT
  
  !--- Arrays for TKE equation ---!
  real(mytype), save, allocatable, dimension(:)     :: tke_convHT      ! convective term: d(<k> U) / dy
  real(mytype), save, allocatable, dimension(:)     :: tke_diffHT      ! diffusive transport term: -nu * d^2 (<k>) / dy^2
  real(mytype), save, allocatable, dimension(:)     :: tke_prodHT      ! production term: - <u'v'> dU/dy
  real(mytype), save, allocatable, dimension(:)     :: temp_dery       ! temporary variable to store the derivative in y of a generic quantity 
     
contains

  !******************************************************************
  subroutine init_post_variables

    USE var
    
    integer :: i,j,k

    TYPE(DECOMP_INFO), save :: ph  ! decomposition object

    if (nrank==0) print *,'Initializing post-processing variables...'

    if (nclx) then
       nxmsize = xsize(1)
    else
       nxmsize = xsize(1) -1
    endif
    if (ncly) then
       nymsize = ysize(2)
    else
       nymsize = ysize(2) -1
    endif
    if (nclz) then
       nzmsize = zsize(3)
    else
       nzmsize = zsize(3) -1
    endif
    
    call decomp_info_init(nxmsize, nymsize, nzmsize, ph)
        
    !xsize(i), ysize(i), zsize(i), i=1,2,3 - sizes of the sub-domains held by the current process. The first letter refers to the pencil orientation and the three 1D array elements contain the sub-domain sizes in X, Y and Z directions, respectively. In a 2D pencil decomposition, there is always one dimension which completely resides in local memory. So by definition xsize(1)==nx_global, ysize(2)==ny_global and zsize(3)==nz_global.

    !xstart(i), ystart(i), zstart(i), xend(i), yend(i), zend(i), i=1,2,3 - the starting and ending indices for each sub-domain, as in the global coordinate system. Obviously, it can be seen that xsize(i)=xend(i)-xstart(i)+1. It may be convenient for certain applications to use global coordinate (for example when extracting a 2D plane from a 3D domain, it is easier to know which process owns the plane if global index is used).
    
    ! Allocate x-pencils arrays using global indices (temporary array is not necessary as x-pencil)
    call alloc_x(ux1, opt_global=.true.)  !global indices
    ux1 = zero
    call alloc_x(uy1, opt_global=.true.)  !global indices
    uy1 = zero
    call alloc_x(uz1, opt_global=.true.)  !global indices
    uz1 = zero
    call alloc_x(pre1, opt_global=.true.) !global indices
    pre1 = zero    
    allocate(phi1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices
    phi1 = zero
    
    ! Allocate y-pencils arrays
    call alloc_y(ux2)
    ux2=zero
    call alloc_y(uy2)
    uy2=zero
    call alloc_y(uz2)
    uz2=zero
    call alloc_y(pre2)
    pre2=zero      
    allocate(phi2(ysize(1),ysize(2),ysize(3)))
    phi2=zero   
    call alloc_y(ta2)
    ta2 = zero
    
    ! Allocate z-pencils arrays
    call alloc_z(ux3)
    ux3=zero
    call alloc_z(uy3)
    uy3=zero
    call alloc_z(uz3)
    uz3=zero

    ! Allocate memory for vorticity & mean gradient calculations (these are gradients)
    if (post_vort) then
       if (.not.allocated(ta1)) call alloc_x(ta1)
       if (.not.allocated(tb1)) call alloc_x(tb1)
       if (.not.allocated(tc1)) call alloc_x(tc1)
       if (.not.allocated(td1)) call alloc_x(td1)
       if (.not.allocated(te1)) call alloc_x(te1)
       if (.not.allocated(tf1)) call alloc_x(tf1)
       if (.not.allocated(tg1)) call alloc_x(tg1)
       if (.not.allocated(th1)) call alloc_x(th1)
       if (.not.allocated(ti1)) call alloc_x(ti1)
       if (.not.allocated(di1)) call alloc_x(di1)

       if (.not.allocated(ta2)) call alloc_y(ta2)
       if (.not.allocated(tb2)) call alloc_y(tb2)
       if (.not.allocated(tc2)) call alloc_y(tc2)
       if (.not.allocated(td2)) call alloc_y(td2)
       if (.not.allocated(te2)) call alloc_y(te2)
       if (.not.allocated(tf2)) call alloc_y(tf2)
       if (.not.allocated(di2)) call alloc_y(di2)
     
       if (.not.allocated(ta3)) call alloc_z(ta3)
       if (.not.allocated(tb3)) call alloc_z(tb3)
       if (.not.allocated(tc3)) call alloc_z(tc3)
       if (.not.allocated(td3)) call alloc_z(td3)
       if (.not.allocated(te3)) call alloc_z(te3)
       if (.not.allocated(tf3)) call alloc_z(tf3)
       if (.not.allocated(di3)) call alloc_z(di3)
    endif
    
    ! Module derivative
    allocate(ffx(nx),sfx(nx),fsx(nx),fwx(nx),ssx(nx),swx(nx))
    allocate(ffxp(nx),sfxp(nx),fsxp(nx),fwxp(nx),ssxp(nx),swxp(nx))
    allocate(ffy(ny),sfy(ny),fsy(ny),fwy(ny),ssy(ny),swy(ny))
    allocate(ffyp(ny),sfyp(ny),fsyp(ny),fwyp(ny),ssyp(ny),swyp(ny))
    allocate(ffz(nz),sfz(nz),fsz(nz),fwz(nz),ssz(nz),swz(nz))
    allocate(ffzp(nz),sfzp(nz),fszp(nz),fwzp(nz),sszp(nz),swzp(nz))

    allocate(ffxS(nx),sfxS(nx),fsxS(nx),fwxS(nx),ssxS(nx),swxS(nx))
    allocate(ffxpS(nx),sfxpS(nx),fsxpS(nx),fwxpS(nx),ssxpS(nx),swxpS(nx))
    allocate(ffyS(ny),sfyS(ny),fsyS(ny),fwyS(ny),ssyS(ny),swyS(ny))
    allocate(ffypS(ny),sfypS(ny),fsypS(ny),fwypS(ny),ssypS(ny),swypS(ny))
    allocate(ffzS(nz),sfzS(nz),fszS(nz),fwzS(nz),sszS(nz),swzS(nz))
    allocate(ffzpS(nz),sfzpS(nz),fszpS(nz),fwzpS(nz),sszpS(nz),swzpS(nz))

    allocate(sx(xsize(2),xsize(3)),vx(xsize(2),xsize(3)))
    allocate(sy(ysize(1),ysize(3)),vy(ysize(1),ysize(3)))
    allocate(sz(zsize(1),zsize(2)),vz(zsize(1),zsize(2)))

    ! Module derpres
    allocate(cfx6(nxm),ccx6(nxm),cbx6(nxm),cfxp6(nxm),ciwxp6(nxm),csxp6(nxm),&
         cwxp6(nxm),csx6(nxm),cwx6(nxm),cifx6(nxm),cicx6(nxm),cisx6(nxm))
    allocate(cibx6(nxm),cifxp6(nxm),cisxp6(nxm),ciwx6(nxm))
    allocate(cfi6(nx),cci6(nx),cbi6(nx),cfip6(nx),csip6(nx),cwip6(nx),csi6(nx),&
         cwi6(nx),cifi6(nx),cici6(nx),cibi6(nx),cifip6(nx))
    allocate(cisip6(nx),ciwip6(nx),cisi6(nx),ciwi6(nx))
    allocate(cfy6(nym),ccy6(nym),cby6(nym),cfyp6(nym),csyp6(nym),cwyp6(nym),csy6(nym))
    allocate(cwy6(nym),cify6(nym),cicy6(nym),ciby6(nym),cifyp6(nym),cisyp6(nym),&
         ciwyp6(nym),cisy6(nym),ciwy6(nym))
    allocate(cfi6y(ny),cci6y(ny),cbi6y(ny),cfip6y(ny),csip6y(ny),cwip6y(ny),&
         csi6y(ny),cwi6y(ny),cifi6y(ny),cici6y(ny))
    allocate(cibi6y(ny),cifip6y(ny),cisip6y(ny),ciwip6y(ny),cisi6y(ny),ciwi6y(ny))
    allocate(cfz6(nzm),ccz6(nzm),cbz6(nzm),cfzp6(nzm),cszp6(nzm),cwzp6(nzm),csz6(nzm))
    allocate(cwz6(nzm),cifz6(nzm),cicz6(nzm),cibz6(nzm),cifzp6(nzm),ciszp6(nzm),&
         ciwzp6(nzm),cisz6(nzm),ciwz6(nzm))
    allocate(cfi6z(nz),cci6z(nz),cbi6z(nz),cfip6z(nz),csip6z(nz),cwip6z(nz),&
         csi6z(nz),cwi6z(nz),cifi6z(nz),cici6z(nz))
    allocate(cibi6z(nz),cifip6z(nz),cisip6z(nz),ciwip6z(nz),cisi6z(nz),ciwi6z(nz))
       
    ! Module mesh 
    allocate(ppy(ny))
    ppy=zero
    allocate(pp2y(ny))
    pp2y=zero
    allocate(pp4y(ny))
    pp4y=zero

    allocate(ppyi(ny))
    ppyi=zero
    allocate(pp2yi(ny))
    pp2yi=zero
    allocate(pp4yi(ny))
    pp4yi=zero

    allocate(yp(ny))
    yp=zero
    allocate(ypi(ny))
    ypi=zero
    allocate(del(ny))
    del=zero

    allocate(yeta(ny))
    yeta=zero
    allocate(yetai(ny))
    yetai=zero

    ! y-position
    if (istret.eq.0) then
       do j=1,ny
          yp(j)=real(j-1,mytype)*dy
          ypi(j)=(real(j,mytype)-half)*dy
          ppy(j) = one
       enddo
       if (ncly1.eq.1 .or. ncly1.eq.2) then
          ppy(1) = two
       endif
       if (nclyn.eq.1 .or. nclyn.eq.2) then
          ppy(ny) = two
       endif
    else
       call stretching()

       allocate(dyp(ny))
       ! compute dy for stretched mesh - Kay
       do j=2,ny-1
          dyp(j) = half*(yp(j+1)-yp(j-1))
       enddo
       dyp(1)  = yp(2) -yp(1)
       dyp(ny) = yp(ny)-yp(ny-1)
    endif
  
  end subroutine init_post_variables
  !******************************************************************
  ! Subroutine for allocating the memory for statistics arrays
  subroutine init_statistics()
  
  USE MPI
  
  implicit none
      
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
        allocate(mean_gradientp(ysize(1),ysize(2),ysize(3))); mean_gradientp = zero       
        allocate(mean_gradientx(ysize(1),ysize(2),ysize(3))); mean_gradientx = zero
        allocate(mean_gradientz(ysize(1),ysize(2),ysize(3))); mean_gradientz = zero
       
        allocate(mean_gradientpH1(ysize(2))); mean_gradientpH1 = zero
        allocate(mean_gradientxH1(ysize(2))); mean_gradientxH1 = zero
        allocate(mean_gradientzH1(ysize(2))); mean_gradientzH1 = zero
       
        allocate(mean_gradientpHT(ysize(2))); mean_gradientpHT = zero
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
        allocate(RppzH1(zsize(2),zsize(3))); RppzH1 = zero 
               
        allocate(RuuzHT(zsize(2),zsize(3))); RuuzHT = zero
        allocate(RvvzHT(zsize(2),zsize(3))); RvvzHT = zero
        allocate(RwwzHT(zsize(2),zsize(3))); RwwzHT = zero
        allocate(RuvzHT(zsize(2),zsize(3))); RuvzHT = zero
        allocate(RppzHT(zsize(2),zsize(3))); RppzHT = zero
    
    end if
    
    if (post_tke_eq) then
       
        ! Fluctuating terms for TKE equation
        allocate(kvprime_mean       (ysize(1),ysize(2),ysize(3))); kvprime_mean        = zero
        allocate(pseudo_eps_tke_mean(ysize(1),ysize(2),ysize(3))); pseudo_eps_tke_mean = zero
        
        allocate(kvprime_meanH1       (ysize(2))); kvprime_meanH1        = zero
        allocate(pseudo_eps_tke_meanH1(ysize(2))); pseudo_eps_tke_meanH1 = zero
        
        allocate(kvprime_meanHT       (ysize(2))); kvprime_meanHT        = zero
        allocate(pseudo_eps_tke_meanHT(ysize(2))); pseudo_eps_tke_meanHT = zero
        
        ! All other terms
        allocate(tke_convHT   (ysize(2))); tke_convHT    = zero
        allocate(tke_diffHT   (ysize(2))); tke_diffHT    = zero     
        allocate(tke_prodHT   (ysize(2))); tke_prodHT    = zero
        allocate(temp_dery    (ysize(2))); temp_dery     = zero
                
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
  
  !******************************************************************
  ! Subroutine to reset to zero the arrays of averages on same position 
  
  subroutine reset_averages()
  
  USE MPI
  
  implicit none
  
  if (post_mean) then
      u1mean=zero;v1mean=zero;w1mean=zero
      u2mean=zero;v2mean=zero;w2mean=zero
      u3mean=zero;v3mean=zero;w3mean=zero
      u4mean=zero;v4mean=zero;w4mean=zero
      uvmean=zero;uwmean=zero;vwmean=zero
      pre1mean=zero;pre2mean=zero
      vpremean=zero
  
      ! Scalar field
      phi1mean=zero;phi2mean=zero
      uphimean=zero;vphimean=zero;wphimean=zero
      
  end if
  
  if (post_vort) then
      vortxmean=zero; vortymean=zero; vortzmean=zero
      mean_gradientp=zero; mean_gradientx=zero; mean_gradientz=zero
      mean_gradphi=zero
  end if
  
  if (post_diss) then
      epsmean=zero
  end if
  
  if(post_tke_eq) then
      kvprime_mean=zero;pseudo_eps_tke_mean=zero
  end if
     
  end subroutine reset_averages

  !******************************************************************
  ! Subroutine to reset to zero the arrays of averages on subdomains
  ! and on total domain 
  
  subroutine reset_subdomains_and_domain()
  
  USE MPI
  
  implicit none
  
  if (post_mean) then
  
      ! Subdomains   
      u1meanH1=zero;v1meanH1=zero;w1meanH1=zero
      u2meanH1=zero;v2meanH1=zero;w2meanH1=zero
      u3meanH1=zero;v3meanH1=zero;w3meanH1=zero
      u4meanH1=zero;v4meanH1=zero;w4meanH1=zero
      uvmeanH1=zero;uwmeanH1=zero;vwmeanH1=zero
      pre1meanH1=zero;pre2meanH1=zero
      vpremeanH1=zero
      
      ! Total domain
      u1meanHT=zero;v1meanHT=zero;w1meanHT=zero
      u2meanHT=zero;v2meanHT=zero;w2meanHT=zero
      u3meanHT=zero;v3meanHT=zero;w3meanHT=zero
      u4meanHT=zero;v4meanHT=zero;w4meanHT=zero
      uvmeanHT=zero;uwmeanHT=zero;vwmeanHT=zero
      pre1meanHT=zero;pre2meanHT=zero
      vpremeanHT=zero
            
      ! Scalar field
            
      ! Subdomains   
      phi1meanH1=zero;phi2meanH1=zero
      uphimeanH1=zero;vphimeanH1=zero;wphimeanH1=zero
      
      ! Total domain
      phi1meanHT=zero;phi2meanHT=zero
      uphimeanHT=zero;vphimeanHT=zero;wphimeanHT=zero
           
  end if
  
  if (post_vort) then
  
      ! Subdomains
      vortxmeanH1=zero; vortymeanH1=zero; vortzmeanH1=zero
      mean_gradientpH1=zero; mean_gradientxH1=zero; mean_gradientzH1=zero
      mean_gradphiH1=zero
  
      ! Total domain
      vortxmeanHT=zero; vortymeanHT=zero; vortzmeanHT=zero
      mean_gradientpHT=zero; mean_gradientxHT=zero; mean_gradientzHT=zero
      mean_gradphiHT=zero
  
  end if
  
  if (post_diss) then
  
      ! Subdomains
      epsmeanH1=zero
  
      ! Total domain
      epsmeanHT=zero
  
  end if
  
  if (post_corz) then
  
      ! Subdomains
      RuuzH1=zero;RvvzH1=zero;RwwzH1=zero;RuvzH1=zero;RppzH1=zero
  
      ! Total domain
      RuuzHT=zero;RvvzHT=zero;RwwzHT=zero;RuvzHT=zero;RppzHT=zero
  
  end if
  
  if(post_tke_eq) then
  
      ! Subdomains
      kvprime_meanH1=zero;pseudo_eps_tke_meanH1=zero
  
      ! Total domain
      kvprime_meanHT=zero;pseudo_eps_tke_meanHT=zero
      
      ! All the other terms are only in total domain
      tke_convHT    = zero
      tke_diffHT    = zero     
      tke_prodHT    = zero
      temp_dery     = zero
  
  end if
  
  end subroutine reset_subdomains_and_domain
     
end module post_processing
!********************************************************************




