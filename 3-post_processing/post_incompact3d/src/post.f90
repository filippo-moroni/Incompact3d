!----------------------------------------------------------!
!         This program is used for post-processing         !
!                 of Incompact3d Snapshots.                !
!      Adapted from original Incompact3d file (v2.0)       !
!                    of R. Corsini                         !
!----------------------------------------------------------!

PROGRAM post

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE var
  USE MPI
  USE post_processing
  USE tools

  implicit none

  integer :: i,j,k,is
  integer :: ii = 1,ie = 1                         ! internal and external loops 
  integer :: file1,filen,icrfile                   ! indexes for opening snapshots (first, last & increment)
  integer :: nt                                    ! total number of time units
  integer :: nr                                    ! total number of flow realizations
  integer :: ifile                                 ! index to open different snapshots in time
  
  real(mytype) :: tstart=0.0,tend=0.0,ttotal=0.0   ! variables to count time spent to post-process data
   
  integer :: iunit		      		   ! unit for the file to open (assigned by the compiler)
    
  integer,dimension(3) :: sel                      ! index for the number of post-processing subroutines employed (selector index)
  logical :: read_vel,read_pre,read_phi  
 
  character(99):: filename,dirname,snap_index 
  character(1) :: a
  
  ! Integer for MPI
  integer :: code
    
  ! Initialize MPI
  call MPI_INIT(code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,code) 
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,code)    ! same as R. Corsini & Xcompact3d
  
  ! Reading of input.i3d file
  call reading_input_file()
  
  ! Imposing decomposition in slabs
  !p_row=1; p_col=nproc
  
  ! Setting up the 2d decomposition
  !call decomp_2d_init(nx,ny,nz,1,nproc)                      ! modified by R. Corsini
  
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_io_init()
  
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)      !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)      !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.)   !start from 1 == true
  
  call init_post_variables()
  call schemes()
  call decomp_info_init(nxm,nym,nzm,phG)
  
 
  ! Start of the post-processing  
  post_mean=.false.; post_vort=.false.; post_diss=.false.  
  read_vel=.false.;  read_pre=.false.; read_phi=.false.
                  
  ! Reading of the input file of post-processing (post.prm)
  open(10,file='post.prm',status='unknown',form='formatted')
  read (10,'(A1)') a
  read (10,'(A1)') a
  read (10,'(A1)') a
  read (10,*) file1
  read (10,*) filen
  read (10,*) icrfile
  read (10,*) nr 
  read (10,'(A1)') a
  read (10,'(A1)') a
  read (10,'(A1)') a
  do i=1,size(sel)
     read (10,*) sel(i)
  enddo
  close(10)
  
  if (sel(1)==1) post_mean=.true.
  if (sel(2)==1) post_vort=.true.
  if (sel(3)==1) post_diss=.true.

  if (nrank==0) then
     if ((.not.post_mean).and.(.not.post_vort)) &
        call decomp_2d_abort(10,'Invalid post-processing switchers specified, NO WORK TO BE DONE HERE')
  endif
  
  ! Logicals for reading Snapshots
  if (post_mean) then
     read_vel=.true.
     read_pre=.true.
     if (iscalar==1) read_phi=.true. 
  endif

  if (post_vort) read_vel=.true.
  if (post_diss) read_vel=.true.
  
  ! Total number of Snapshots in time
  nt = (filen-file1)/icrfile+1
  
  ! Initialize statistics
  call init_statistics(nt)

  ! Time passed since the program has started
  call cpu_time(tstart)
    
!-----------------------------!
! Post-processing starts here !
!-----------------------------!

!------------Start of the time average cycle---------------!

 do ie=1,nt
      
     ! Number of snapshot  
     ifile = (ie-1)*icrfile+file1
          
     ! Show progress on post-processing    
     if (nrank==0) then
        print *,'----------------------------------------------------'
        write(*,"('We are averaging the realizations of the snapshot =',I3,'/',I3)") ifile,filen
     endif
 
!---------Start of the ensemble average cycle--------------!

  do ii=1,nr 
     
     ! Writing the directories where snapshots are saved (different realizations)
     if(nr .eq. 1) then
     
        ! nr = 1 (only /data folder is present)
        write(dirname,"('./data')")  
     else
        ! nr > 1 (/data_r1, /data_r2, etc. folders) 
        write(dirname,"('./data_r',I1.1)") ii
     end if
     
     ! Reading of velocity, pressure and scalar fields if required     
     if (read_vel) then  
          
        ! Reading the x-pencils (snapshots have metadata of the default saving along x)    
        write(filename,"('ux-',I3.3,'.bin')") ifile           
        call decomp_2d_read_one(1,ux1,dirname,filename,a)
  
        write(filename,"('uy-',I3.3,'.bin')") ifile           
        call decomp_2d_read_one(1,uy1,dirname,filename,a) 

        write(filename,"('uz-',I3.3,'.bin')") ifile           
        call decomp_2d_read_one(1,uz1,dirname,filename,a)    
        
        ! Check if divergent values are present             
        call test_speed_min_max(ux1,uy1,uz1)
     endif
     
     if (read_pre) then    
        write(filename,"('pp-',I3.3,'.bin')") ifile            
        call decomp_2d_read_one(1,pre1,dirname,filename,a)        
     endif
     
     if (read_phi) then      
        write(filename,"('phi01-',I3.3,'.bin')") ifile            
        call decomp_2d_read_one(1,phi1(:,:,:),dirname,filename,a)           
     endif
        
     ! Transpose data to y-pencils 
     call transpose_x_to_y(ux1,ux2)
     call transpose_x_to_y(uy1,uy2)
     call transpose_x_to_y(uz1,uz2)
     call transpose_x_to_y(pre1,pre2)
     call transpose_x_to_y(phi1,phi2)

     ! Statistics computation through external subroutines
     if (post_mean) call stat_mean(ux2,uy2,uz2,pre2,phi2,nr,nt, &
                                   u1mean,v1mean,w1mean,u2mean,v2mean,w2mean, &
                                   u3mean,v3mean,w3mean,u4mean,v4mean,w4mean, &
                                   uvmean,uwmean,vwmean,pre1mean,pre2mean,phi1mean, &
                                   phi2mean,uphimean,vphimean,wphimean)
                                                                          
     if (post_vort) call stat_vorticity(ux1,uy1,uz1,nr,nt,vortxmean,vortymean,vortzmean,mean_gradient)
     
     if (post_diss) call stat_dissipation(ux1,uy1,uz1,nr,nt,epsmean)

  enddo ! closing of the do-loop on the different flow realizations

#ifdef TTBL_MODE

#else  
  ! Closing of the do-loop for the different time units (or snapshots) (ie index)
  enddo 
#endif
  
!----Mean over homogeneous directions (H = Homogeneous)----!

  ! Summation over x and z directions
  if (post_mean) then
     do k=ystart(3),yend(3)
        do i=ystart(1),yend(1)
           do j=ystart(2),yend(2)          
              u1meanH1(j)=u1meanH1(j)+u1mean(i,j,k)/real(nx*nz,mytype)
              v1meanH1(j)=v1meanH1(j)+v1mean(i,j,k)/real(nx*nz,mytype)
              w1meanH1(j)=w1meanH1(j)+w1mean(i,j,k)/real(nx*nz,mytype)
              u2meanH1(j)=u2meanH1(j)+u2mean(i,j,k)/real(nx*nz,mytype)
              v2meanH1(j)=v2meanH1(j)+v2mean(i,j,k)/real(nx*nz,mytype)
              w2meanH1(j)=w2meanH1(j)+w2mean(i,j,k)/real(nx*nz,mytype)
              u3meanH1(j)=u3meanH1(j)+u3mean(i,j,k)/real(nx*nz,mytype)
              v3meanH1(j)=v3meanH1(j)+v3mean(i,j,k)/real(nx*nz,mytype)
              w3meanH1(j)=w3meanH1(j)+w3mean(i,j,k)/real(nx*nz,mytype)
              u4meanH1(j)=u4meanH1(j)+u4mean(i,j,k)/real(nx*nz,mytype)
              v4meanH1(j)=v4meanH1(j)+v4mean(i,j,k)/real(nx*nz,mytype)
              w4meanH1(j)=w4meanH1(j)+w4mean(i,j,k)/real(nx*nz,mytype)
              uvmeanH1(j)=uvmeanH1(j)+uvmean(i,j,k)/real(nx*nz,mytype)
              uwmeanH1(j)=uwmeanH1(j)+uwmean(i,j,k)/real(nx*nz,mytype)
              vwmeanH1(j)=vwmeanH1(j)+vwmean(i,j,k)/real(nx*nz,mytype)
              pre1meanH1(j)=pre1meanH1(j)+pre1mean(i,j,k)/real(nx*nz,mytype)
              pre2meanH1(j)=pre2meanH1(j)+pre2mean(i,j,k)/real(nx*nz,mytype)                                                   
              phi1meanH1(j)=phi1meanH1(j)+phi1mean(i,j,k)/real(nx*nz,mytype)
              phi2meanH1(j)=phi2meanH1(j)+phi2mean(i,j,k)/real(nx*nz,mytype)
              uphimeanH1(j)=uphimeanH1(j)+uphimean(i,j,k)/real(nx*nz,mytype)
              vphimeanH1(j)=vphimeanH1(j)+vphimean(i,j,k)/real(nx*nz,mytype)
              wphimeanH1(j)=wphimeanH1(j)+wphimean(i,j,k)/real(nx*nz,mytype)                                  
           enddo          
        enddo
     enddo
  endif
  
  if (post_vort) then
     do k=ystart(3),yend(3)
        do i=ystart(1),yend(1)
           do j=ystart(2),yend(2)          
              vortxmeanH1(j)=vortxmeanH1(j)+vortxmean(i,j,k)/real(nx*nz,mytype)
              vortymeanH1(j)=vortymeanH1(j)+vortymean(i,j,k)/real(nx*nz,mytype)
              vortzmeanH1(j)=vortzmeanH1(j)+vortzmean(i,j,k)/real(nx*nz,mytype) 
              mean_gradientH1(j)=mean_gradientH1(j)+mean_gradient(i,j,k)/real(nx*nz,mytype)                  
           enddo
        enddo
     enddo
  endif
  
  if (post_diss) then
     do k=ystart(3),yend(3)
        do i=ystart(1),yend(1)
           do j=ystart(2),yend(2)          
              epsmeanH1(j)=epsmeanH1(j)+epsmean(i,j,k)/real(nx*nz,mytype)                  
           enddo
        enddo
     enddo
  endif

#ifdef TTBL_MODE
   ! Reset to zero the arrays used to collect the averages locally
   call reset_averages()  
#endif
 
!---------Mean over all MPI processes (T = Total)----------!

  if (post_mean) then
     call MPI_REDUCE(u1meanH1,u1meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(v1meanH1,v1meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(w1meanH1,w1meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(u2meanH1,u2meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(v2meanH1,v2meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(w2meanH1,w2meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(u3meanH1,u3meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(v3meanH1,v3meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(w3meanH1,w3meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(u4meanH1,u4meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(v4meanH1,v4meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(w4meanH1,w4meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(uvmeanH1,uvmeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(uvmeanH1,uvmeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(vwmeanH1,vwmeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(pre1meanH1,pre1meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(pre2meanH1,pre2meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)   
     call MPI_REDUCE(phi1meanH1,phi1meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(phi2meanH1,phi2meanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(uphimeanH1,uphimeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(vphimeanH1,vphimeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(wphimeanH1,wphimeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)    
  endif
  
  if (post_vort) then
     call MPI_REDUCE(vortxmeanH1,vortxmeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(vortymeanH1,vortymeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(vortzmeanH1,vortzmeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)  
     call MPI_REDUCE(mean_gradientH1,mean_gradientHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code) 
  endif
  
  if (post_diss) then
     call MPI_REDUCE(epsmeanH1,epsmeanHT,ysize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
  end if

!------------- MPI process nrank = 0 at work --------------!

  ! High-order moments (variance, skewness, kurtosis)
  if(nrank.eq.0) then  
     
     if (post_mean) then
        do j=ystart(2),yend(2)
           u2meanHT(j)=u2meanHT(j)-u1meanHT(j)**2
           v2meanHT(j)=v2meanHT(j)-v1meanHT(j)**2
           w2meanHT(j)=w2meanHT(j)-w1meanHT(j)**2
           u3meanHT(j)=u3meanHT(j)-u1meanHT(j)**3-3*u1meanHT(j)*u2meanHT(j)
           v3meanHT(j)=v3meanHT(j)-v1meanHT(j)**3-3*v1meanHT(j)*v2meanHT(j)
           w3meanHT(j)=w3meanHT(j)-w1meanHT(j)**3-3*w1meanHT(j)*w2meanHT(j)
           u4meanHT(j)=u4meanHT(j)-u1meanHT(j)**4-6*(u1meanHT(j)**2)*u2meanHT(j)-4*u1meanHT(j)*u3meanHT(j)
           v4meanHT(j)=v4meanHT(j)-v1meanHT(j)**4-6*(v1meanHT(j)**2)*v2meanHT(j)-4*v1meanHT(j)*v3meanHT(j)
           w4meanHT(j)=w4meanHT(j)-w1meanHT(j)**4-6*(w1meanHT(j)**2)*w2meanHT(j)-4*w1meanHT(j)*w3meanHT(j)
           uvmeanHT(j)=uvmeanHT(j)-u1meanHT(j)*v1meanHT(j)
           uwmeanHT(j)=uwmeanHT(j)-u1meanHT(j)*w1meanHT(j)
           vwmeanHT(j)=vwmeanHT(j)-v1meanHT(j)*w1meanHT(j)
           pre2meanHT(j)=pre2meanHT(j)-pre1meanHT(j)**2 
           phi2meanHT(j)=phi2meanHT(j)-phi1meanHT(j)**2
           uphimeanHT(j)=uphimeanHT(j)-u1meanHT(j)*phi1meanH1(j)
           vphimeanHT(j)=vphimeanHT(j)-v1meanHT(j)*phi1meanH1(j)
           wphimeanHT(j)=wphimeanHT(j)-w1meanHT(j)*phi1meanH1(j)         
        enddo
     endif

!------------------Write formatted data--------------------!
     
     ! New directory for the statistics
     write(dirname,"('data_post/')") 
        
     call system('mkdir -p '//trim(dirname))

     ! Formatted data
     print *,'----------------------------------------------------'
     write(*,"(1x,'Writing output data in formatted .txt files :')")
     
     ! Mean statistics writing
     if (post_mean) then

#ifdef TTBL_MODE  
        ! Writing the snapshot index as character
        write(snap_index,'(I3.3)') ifile 
        snap_index = adjustl(snap_index) 
        
        ! Write the mean_stats filename for TTBL
        write(filename, '(A,A,A)') 'mean_stats-', trim(snap_index), '.txt'
        filename = adjustl(filename)
#else
        ! Write the mean_stats filename for channel flow
        write(filename, '(A)') 'mean_stats.txt'
        filename = adjustl(filename)
#endif
        
        ! Open the file and write
        open(newunit=iunit,file=trim(dirname)//trim(filename),form='formatted')
        
        do j = ystart(2),yend(2) + 1
        
        if (j .eq. 1) then
        
        write (iunit, *) 'mean[u]'  , ',', 'mean[v]'  , ',', 'mean[w]', ',', &
                         'var[u]'   , ',', 'var[v]'   , ',', 'var[w]' , ',', &
                         'skew[u]'  , ',', 'skew[v]'  , ',', 'skew[w]', ',', &
                         'kurt[u]'  , ',', 'kurt[v]'  , ',', 'kurt[w]', ',', &
                         "<u'v'>"   , ',', "<u'w'>"   , ',', "<v'w'>" , ',', &
                         'mean[p]'  , ',', 'var[p]'   , ',',                 &
                         'mean[phi]', ',', 'var[phi]' , ',',                 &
                         "<u'phi'>" , ',', "<v'phi'>" , ',', "<w'phi'>" 
                               
        else
        
        write(iunit, *)  u1meanHT(j-1),         ',', &
                         v1meanHT(j-1),         ',', &       
                         w1meanHT(j-1),         ',', &
                         u2meanHT(j-1),         ',', &
                         v2meanHT(j-1),         ',', &
                         w2meanHT(j-1),         ',', &
                         u3meanHT(j-1),         ',', &
                         v3meanHT(j-1),         ',', &
                         w3meanHT(j-1),         ',', &
                         u4meanHT(j-1),         ',', &
                         v4meanHT(j-1),         ',', &
                         w4meanHT(j-1),         ',', &
                         uvmeanHT(j-1),         ',', &
                         uwmeanHT(j-1),         ',', &  
                         vwmeanHT(j-1),         ',', &                                              
                         pre1meanHT(j-1),       ',', &
                         pre2meanHT(j-1),       ',', &                     
                         phi1meanHT(j-1),       ',', &
                         phi2meanHT(j-1),       ',', &                        
                         uphimeanHT(j-1),       ',', &
                         vphimeanHT(j-1),       ',', &
                         wphimeanHT(j-1)
        
        end if
        
        end do
                               
        close(iunit)
     endif
     
     ! Vorticity mean statistics and mean gradient writing
     if (post_vort) then

#ifdef TTBL_MODE  
        ! Writing the snapshot index as character
        write(snap_index,'(I3.3)') ifile 
        snap_index = adjustl(snap_index) 
        
        ! Write the vort_stats filename for TTBL
        write(filename, '(A,A,A)') 'vort_stats-', trim(snap_index), '.txt'
        filename = adjustl(filename)
#else
        ! Write the vort_stats filename for channel flow
        write(filename, '(A)') 'vort_stats.txt'
        filename = adjustl(filename)
#endif
               
        ! Open the file and write      
        open(newunit=iunit,file=trim(dirname)//trim(filename),form='formatted')
        
        do j = ystart(2),yend(2) + 1
        
        if (j .eq. 1) then
        
        write (iunit, *) 'mean[omega_x]'  , ',', 'mean[omega_y]'  , ',', 'mean[omega_z]', ',', 'dU/dy'
                               
        else
        
        write(iunit, *)  vortxmeanHT(j-1),    ',', &
                         vortymeanHT(j-1),    ',', &       
                         vortzmeanHT(j-1),    ',', &
                         mean_gradientHT(j-1)
        
        end if
        
        end do
                               
        close(iunit)
     endif
     
     ! Mean dissipation writing
     if (post_diss) then

#ifdef TTBL_MODE  
        ! Writing the snapshot index as character
        write(snap_index,'(I3.3)') ifile
        snap_index = adjustl(snap_index) 
        
        ! Write the diss_stats filename for TTBL
        write(filename, '(A,A,A)') 'diss_stats-', trim(snap_index), '.txt'
        filename = adjustl(filename)
#else
        ! Write the diss_stats filename for channel flow
        write(filename, '(A)') 'diss_stats.txt'
        filename = adjustl(filename)
#endif
               
        ! Open the file and write      
        open(newunit=iunit,file=trim(dirname)//trim(filename),form='formatted')
        
        do j = ystart(2),yend(2) + 1
        
        if (j .eq. 1) then
        
        write (iunit, *) 'mean[eps]'
                               
        else
        
        write(iunit, *)  epsmeanHT(j-1)
        
        end if
        
        end do
                               
        close(iunit)
  endif
     
  endif ! closing of the if-statement for processor 0

#ifdef TTBL_MODE   
   ! Reset to zero the average vectors on subdomains (H1) and on total domain (HT)
   call reset_subdomains_and_domain() 
   
   ! Closing of the do-loop for the different time units (or SnapShots) (ie index)
   enddo  
#endif
     
  !-----------------------------!
  !  Post-processing ends here  !
  !-----------------------------!

  ! Time passed since the post-processing has started 
  call cpu_time(tend)
  
  ! Net time used to perform post-processing
  ttotal=tend-tstart

  ! Print useful informations on the screen
  if (nrank==0) then
     print *,'==========================================================='
     print *,''
     print *,'Post-processing finished successfully!'
     print *,''
#ifdef TTBL_MODE
     print *,'!--- Temporal TBL mode ---!'
#else
     print *,'!--- Channel flow mode ---!'
#endif
     print *,''
     print *,'2DECOMP with p_row*p_col=',p_row,p_col
     print *,''
     print *,'nx*ny*nz=',nx*ny*nz
     print *,'nx,ny,nz=',nx,ny,nz
     print *,'dx,dy,dz=',dx,dy,dz
     print *,''
     print *,'Averaged time per snapshot (s):',real(ttotal/nt,4)
     print *,'Total wallclock (s):',real(ttotal,4)
     print *,'Total wallclock (m):',real(ttotal/60.,4)
     print *,'Total wallclock (h):',real(ttotal/3600.,4)
     print *,'Total wallclock (d):',real(ttotal*1.1574e-5,4)
     print *,''
     
     ! Mean statistics
     if (post_mean) then
     print *,'==========================================================='
     print *,''
     print *,'The following statistics have been saved in'
     print *,'"mean_stats" file(s):'
     print *,''
     print *,'mean[u], mean[v], mean[w]'
     print *,' var[u],  var[v],  var[w]'
     print *,'skew[u], skew[v], skew[w]'
     print *,'kurt[u], kurt[v], kurt[w]'
     print *,''
     print *,'mean[uv], mean[uw], mean[vw]'
     print *,''
     print *,'mean[p],   var[p]'
     print *,'mean[phi], var[phi]'
     print *,''
     print *,"mean[u'phi'], mean[v'phi'], mean[w'phi']"
     print *,''    
     endif
     
     ! Vorticity and mean gradient 
     if (post_vort) then
     print *,'==========================================================='
     print *,''
     print *,'The following statistics have been saved in'
     print *,'"vort_stats" file(s):'
     print *,''
     print *,'mean[omega_x], mean[omega_y], mean[omega_z], dU/dy'
     print *,''    
     endif
     
     ! Total dissipation rate 
     if (post_diss) then
     print *,'==========================================================='
     print *,''
     print *,'The following statistics have been saved in'
     print *,'"diss_stats" file(s):'
     print *,''
     print *,'mean[eps]'
     print *,''    
     endif
     
     print *,'==========================================================='
     
  endif

  call decomp_2d_finalize
  call MPI_FINALIZE(code)

end program post





