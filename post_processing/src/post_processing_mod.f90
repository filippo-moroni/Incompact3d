!********************************************************************
module post_processing

  USE decomp_2d
  USE variables
  USE param

  implicit none
  !
  integer :: FS
  character(len=100) :: fileformat
  character(len=1), parameter :: NL=char(10) !new line character
  !
  integer, save :: xplane1,xjump,nix,rsize,rsizee

  logical, save :: post_mean,post_vort

  !arrays for statistic collection
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

    integer :: j

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
  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !X PENCILS
    call alloc_x(ux1, opt_global=.true.) !global indices
    call alloc_x(uy1, opt_global=.true.) !global indices
    call alloc_x(uz1, opt_global=.true.) !global indices
    ux1=zero;uy1=zero;uz1=zero

    if (post_mean) then
       if (.not.allocated(pre1)) then
          call alloc_x(pre1, opt_global=.true.) !global indices 
          pre1=zero
       endif
       if (.not.allocated(ta1)) call alloc_x(ta1)
       if (iscalar==1) then
           if (.not.allocated(phi1)) then
              allocate(phi1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:numscalar)) !global indices
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
  subroutine init_statistics

    USE MPI

    if (post_mean) then
       allocate(u1mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(v1mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(w1mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(u2mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(v2mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(w2mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(u3mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(v3mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(w3mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(u4mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(v4mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(w4mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(uvmean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(uwmean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(vwmean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(pre1mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   
       allocate(pre2mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) !global indices   

       allocate(u1meanH1(xsize(1)));allocate(v1meanH1(xsize(1)));allocate(w1meanH1(xsize(1)))
       allocate(u2meanH1(xsize(1)));allocate(v2meanH1(xsize(1)));allocate(w2meanH1(xsize(1)))
       allocate(u3meanH1(xsize(1)));allocate(v3meanH1(xsize(1)));allocate(w3meanH1(xsize(1)))
       allocate(u4meanH1(xsize(1)));allocate(v4meanH1(xsize(1)));allocate(w4meanH1(xsize(1)))
       allocate(uvmeanH1(xsize(1)));allocate(uwmeanH1(xsize(1)));allocate(vwmeanH1(xsize(1)))
       allocate(pre1meanH1(xsize(1)));allocate(pre2meanH1(xsize(1)))
   
       allocate(u1meanHT(xsize(1)));allocate(v1meanHT(xsize(1)));allocate(w1meanHT(xsize(1)))
       allocate(u2meanHT(xsize(1)));allocate(v2meanHT(xsize(1)));allocate(w2meanHT(xsize(1)))
       allocate(u3meanHT(xsize(1)));allocate(v3meanHT(xsize(1)));allocate(w3meanHT(xsize(1)))
       allocate(u4meanHT(xsize(1)));allocate(v4meanHT(xsize(1)));allocate(w4meanHT(xsize(1)))
       allocate(uvmeanHT(xsize(1)));allocate(uwmeanHT(xsize(1)));allocate(vwmeanHT(xsize(1)))
       allocate(pre1meanHT(xsize(1)));allocate(pre2meanHT(xsize(1)))

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
          allocate(phi1mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:numscalar)) !global indices   
          allocate(phi2mean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:numscalar)) !global indices   
          allocate(uphimean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:numscalar)) !global indices   
          allocate(vphimean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:numscalar)) !global indices   
          allocate(wphimean(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3),1:numscalar)) !global indices   

          allocate(phi1meanH1(xsize(1),1:numscalar));allocate(phi2meanH1(xsize(1),1:numscalar))
          allocate(uphimeanH1(xsize(1),1:numscalar));allocate(vphimeanH1(xsize(1),1:numscalar));allocate(wphimeanH1(xsize(1),1:numscalar))
          allocate(phi1meanHT(xsize(1),1:numscalar));allocate(phi2meanHT(xsize(1),1:numscalar))
          allocate(uphimeanHT(xsize(1),1:numscalar));allocate(vphimeanHT(xsize(1),1:numscalar));allocate(wphimeanHT(xsize(1),1:numscalar))

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

