!********************************************************************
module post_processing

  USE decomp_2d
  USE variables
  USE param

  implicit none
  
  ! Variables declaration
  integer :: FS
  character(len=100) :: fileformat
  character(len=1), parameter :: NL=char(10) !new line character

  logical, save :: post_mean,post_vort

  ! Arrays for statistic collection
  real(mytype), save, allocatable, dimension(:,:,:) :: u1mean,v1mean,w1mean
  real(mytype), save, allocatable, dimension(:,:,:) :: u2mean,v2mean,w2mean
  real(mytype), save, allocatable, dimension(:,:,:) :: u3mean,v3mean,w3mean
  real(mytype), save, allocatable, dimension(:,:,:) :: u4mean,v4mean,w4mean
  real(mytype), save, allocatable, dimension(:,:,:) :: uvmean,uwmean,vwmean
  real(mytype), save, allocatable, dimension(:,:,:) :: pre1mean,pre2mean
  real(mytype), save, allocatable, dimension(:,:,:,:) :: phi1mean,phi2mean
  real(mytype), save, allocatable, dimension(:,:,:,:) :: uphimean,vphimean,wphimean

  real(mytype), save, allocatable, dimension(:) :: u1meanH1,v1meanH1,w1meanH1
  real(mytype), save, allocatable, dimension(:) :: u2meanH1,v2meanH1,w2meanH1
  real(mytype), save, allocatable, dimension(:) :: u3meanH1,v3meanH1,w3meanH1
  real(mytype), save, allocatable, dimension(:) :: u4meanH1,v4meanH1,w4meanH1
  real(mytype), save, allocatable, dimension(:) :: uvmeanH1,uwmeanH1,vwmeanH1
  real(mytype), save, allocatable, dimension(:) :: pre1meanH1,pre2meanH1
  real(mytype), save, allocatable, dimension(:,:) :: phi1meanH1,phi2meanH1
  real(mytype), save, allocatable, dimension(:,:) :: uphimeanH1,vphimeanH1,wphimeanH1

  real(mytype), save, allocatable, dimension(:) :: u1meanHT,v1meanHT,w1meanHT
  real(mytype), save, allocatable, dimension(:) :: u2meanHT,v2meanHT,w2meanHT
  real(mytype), save, allocatable, dimension(:) :: u3meanHT,v3meanHT,w3meanHT
  real(mytype), save, allocatable, dimension(:) :: u4meanHT,v4meanHT,w4meanHT
  real(mytype), save, allocatable, dimension(:) :: uvmeanHT,uwmeanHT,vwmeanHT
  real(mytype), save, allocatable, dimension(:) :: pre1meanHT,pre2meanHT
  real(mytype), save, allocatable, dimension(:,:) :: phi1meanHT,phi2meanHT
  real(mytype), save, allocatable, dimension(:,:) :: uphimeanHT,vphimeanHT,wphimeanHT

contains

  !******************************************************************
  subroutine init_post_variables

    USE var

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


    !Y PENCILS
    call alloc_y(ux1, opt_global=.true.)  !global indices
    ux1 = zero
    call alloc_y(uy1, opt_global=.true.)  !global indices
    uy1 = zero
    call alloc_y(uz1, opt_global=.true.)  !global indices
    uz1 = zero

    call alloc_y(pre1, opt_global=.true.) !global indices
    pre1 = zero

    allocate(phi1(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),1:numscalar)) !global indices
    phi1 = zero

    call alloc_y(ta1)
    ta1 = zero
       
    ! Check if memory is allocated properly
    if (post_mean) then
       if (.not.allocated(pre1)) then
          call alloc_y(pre1, opt_global=.true.) !global indices 
          pre1=zero
       endif
       if (.not.allocated(ta1)) call alloc_y(ta1)
       if (iscalar==1) then
           if (.not.allocated(phi1)) then
              allocate(phi1(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),1:numscalar)) !global indices
              phi1=zero
           endif
       endif
    endif

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

  end subroutine init_post_variables
  !******************************************************************
  ! Subroutine for allocating the memory for arrays
  subroutine init_statistics

    USE MPI
    
    if (post_mean) then
       allocate(u1mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))  ! global indices   
       allocate(v1mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))   
       allocate(w1mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))   
       allocate(u2mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))    
       allocate(v2mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))    
       allocate(w2mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))    
       allocate(u3mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))    
       allocate(v3mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))   
       allocate(w3mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))    
       allocate(u4mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))  
       allocate(v4mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))    
       allocate(w4mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))      
       allocate(uvmean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))      
       allocate(uwmean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))       
       allocate(vwmean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))     
       allocate(pre1mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))    
       allocate(pre2mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))  
       
       allocate(u1meanH1(ysize(2)));   allocate(v1meanH1(ysize(2)));  allocate(w1meanH1(ysize(2)))
       allocate(u2meanH1(ysize(2)));   allocate(v2meanH1(ysize(2)));  allocate(w2meanH1(ysize(2)))
       allocate(u3meanH1(ysize(2)));   allocate(v3meanH1(ysize(2)));  allocate(w3meanH1(ysize(2)))
       allocate(u4meanH1(ysize(2)));   allocate(v4meanH1(ysize(2)));  allocate(w4meanH1(ysize(2)))
       allocate(uvmeanH1(ysize(2)));   allocate(uwmeanH1(ysize(2)));  allocate(vwmeanH1(ysize(2)))
       allocate(pre1meanH1(ysize(2))); allocate(pre2meanH1(ysize(2)))
   
       allocate(u1meanHT(ysize(2)));   allocate(v1meanHT(ysize(2)));  allocate(w1meanHT(ysize(2)))
       allocate(u2meanHT(ysize(2)));   allocate(v2meanHT(ysize(2)));  allocate(w2meanHT(ysize(2)))
       allocate(u3meanHT(ysize(2)));   allocate(v3meanHT(ysize(2)));  allocate(w3meanHT(ysize(2)))
       allocate(u4meanHT(ysize(2)));   allocate(v4meanHT(ysize(2)));  allocate(w4meanHT(ysize(2)))
       allocate(uvmeanHT(ysize(2)));   allocate(uwmeanHT(ysize(2)));  allocate(vwmeanHT(ysize(2)))
       allocate(pre1meanHT(ysize(2))); allocate(pre2meanHT(ysize(2)))

       u1mean=zero;v1mean=zero;w1mean=zero
       u2mean=zero;v2mean=zero;w2mean=zero
       u3mean=zero;v3mean=zero;w3mean=zero
       u4mean=zero;v4mean=zero;w4mean=zero
       uvmean=zero;uwmean=zero;vwmean=zero
       pre1mean=zero;pre2mean=zero
   
       u1meanH1=zero;v1meanH1=zero;w1meanH1=zero
       u2meanH1=zero;v2meanH1=zero;w2meanH1=zero
       u3meanH1=zero;v3meanH1=zero;w3meanH1=zero
       u4meanH1=zero;v4meanH1=zero;w4meanH1=zero
       uvmeanH1=zero;uwmeanH1=zero;vwmeanH1=zero
       pre1meanH1=zero;pre2meanH1=zero
   
       u1meanHT=zero;v1meanHT=zero;w1meanHT=zero
       u2meanHT=zero;v2meanHT=zero;w2meanHT=zero
       u3meanHT=zero;v3meanHT=zero;w3meanHT=zero
       u4meanHT=zero;v4meanHT=zero;w4meanHT=zero
       uvmeanHT=zero;uwmeanHT=zero;vwmeanHT=zero
       pre1meanHT=zero;pre2meanHT=zero

      if (iscalar==1) then
          allocate(phi1mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),1:numscalar))  ! global indices   
          allocate(phi2mean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),1:numscalar))    
          allocate(uphimean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),1:numscalar))    
          allocate(vphimean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),1:numscalar))   
          allocate(wphimean(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3),1:numscalar))    

          allocate(phi1meanH1(ysize(2),1:numscalar)); allocate(phi2meanH1(ysize(2),1:numscalar))
          allocate(uphimeanH1(ysize(2),1:numscalar)); allocate(vphimeanH1(ysize(2),1:numscalar)); allocate(wphimeanH1(ysize(2),1:numscalar))
          allocate(phi1meanHT(ysize(2),1:numscalar)); allocate(phi2meanHT(ysize(2),1:numscalar))
          allocate(uphimeanHT(ysize(2),1:numscalar)); allocate(vphimeanHT(ysize(2),1:numscalar)); allocate(wphimeanHT(ysize(2),1:numscalar))

          phi1mean=zero;phi2mean=zero
          uphimean=zero;vphimean=zero;wphimean=zero
  
          phi1meanH1=zero;phi2meanH1=zero
          uphimeanH1=zero;vphimeanH1=zero;wphimeanH1=zero
  
          phi1meanHT=zero;phi2meanHT=zero
          uphimeanHT=zero;vphimeanHT=zero;wphimeanHT=zero
       endif
    endif

  end subroutine init_statistics
  
end module post_processing
!********************************************************************




