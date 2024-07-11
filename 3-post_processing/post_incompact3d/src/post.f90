!----------------------------------------------------------!
!         This program is used for post-processing         !
!                 of Incompact3d Snapshots.                !
!      Adapted from original Incompact3d file (v2.0)       !
!                    of R. Corsini                         !
!----------------------------------------------------------!

program post

  use decomp_2d
  use decomp_2d_io
  use variables
  use param
  use var
  use MPI
  use post_processing
  use tools

  implicit none

  integer :: i,j,k,is
  integer :: ii = 1,ie = 1                           ! internal and external loops 
  integer :: file1,filen,icrfile                     ! indexes for opening snapshots (first, last & increment)
  integer :: nt                                      ! total number of time units
  integer :: nr                                      ! total number of flow realizations
  integer :: ifile                                   ! index to open different snapshots in time
  
  real(mytype) :: tstart=0.0,tend=0.0,ttotal=0.0     ! variables to count time spent to post-process data
  real(mytype) :: den                                ! denominator of the divisions
   
  integer :: iunit		      		     ! unit for the file to open (assigned by the compiler)
    
  integer,dimension(4) :: sel                        ! index for the number of post-processing subroutines employed (selector index)
  logical :: read_vel,read_pre,read_phi  
 
  character(99) :: filename,dirname
  character(99) :: snap_index,snap_n_index,printing  ! characters to print to screen and to read snapshots' indexes
  
  character(1)  :: a
  
  ! Format for snapshots numbers
  character(len=9) :: ifilenameformat = '(I4.4)'
  
  ! Format for correlations
  character(len=50) :: format_string
  
  ! Integer for MPI
  integer :: code
  
  ! Save the initial time of post-processing work
  call cpu_time(tstart)
    
  ! Initialize MPI (same as R. Corsini & Xcompact3d)
  call MPI_INIT(code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,code) 
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,code)     
  
  ! Reading of input.i3d file
  call reading_input_file()
  
  ! Imposing decomposition in slabs
  p_row=nproc; p_col=1
  
  ! Setting up the 2d decomposition
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  
  call decomp_2d_io_init()
  
  ! Initialize post-processing variables
  call init_post_variables()
  
  ! Call subroutines from the decomp2d library 
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)      ! start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)      ! start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.)   ! start from 1 == true
  
  call schemes()
  call decomp_info_init(nxm,nym,nzm,phG)
  
  ! Start of the post-processing  
  post_mean=.false.; post_vort=.false.; post_diss=.false.; post_corz=.false.  
  read_vel=.false.;  read_pre=.false.;  read_phi=.false.
                  
  ! Reading of the input file of post-processing (post.prm)
  open(10,file='post.prm',status='old',form='formatted')
  read (10,'(A1)') a
  read (10,'(A1)') a
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
  if (sel(4)==1) post_corz=.true.

  if (nrank==0) then
     if ((.not.post_mean).and.(.not.post_vort).and.(.not.post_diss).and.(.not.post_corz)) &
        call decomp_2d_abort(10,'Invalid post-processing switchers specified, no work to be done here!')
  endif
  
  ! Logicals for reading Snapshots
  if (post_mean) then
     read_vel=.true.
     read_pre=.true.
  endif
  
  ! Reading of velocity only if necessary
  if (post_vort .or. post_diss .or. post_corz) read_vel=.true.
  
  ! Read of scalar field only if necessary
  if (iscalar==1) read_phi=.true. 
  
  ! Total number of Snapshots in time
  nt = (filen-file1)/icrfile+1
  
  ! Initialize statistics
  call init_statistics()
  
#ifndef TTBL_MODE

  ! Reading of previously calculated mean velocity field for a channel,
  ! if correlation needs to be calculated.

  if (post_corz) then
 
  ! Write directory name
  write(dirname,"('data_post/')")
  
  ! Write the mean_stats filename for channel flow
  write(filename, '(A)') 'mean_stats.txt'
  filename = adjustl(filename)
  
  ! Display that we are reading the mean velocity field file
  if (nrank.eq.0) write(*,"(1x,'Reading mean_stats.txt')")
  
      ! Open the file and read
      open(newunit=iunit,file=trim(dirname)//trim(filename),form='formatted',status='old')
  
      read(iunit, *)
  
      do j = 1, ysize(2)
  
          read(iunit, '(3(F13.9, A1, 1X))') u1meanHT(j), a, &
                                            v1meanHT(j), a, &       
                                            w1meanHT(j), a
      end do
                               
      close(iunit)
  
  end if

#endif
  
!-----------------------------!
! Post-processing starts here !
!-----------------------------!
     
! Writing the final snapshot index as character
if (nrank==0) then
    write(snap_n_index, ifilenameformat) filen 
    snap_n_index = adjustl(snap_n_index)
end if

!------------Start of the time average cycle---------------!

 do ie=1,nt
      
     ! Number of snapshot  
     ifile = (ie-1)*icrfile+file1
     
     ! Writing the snapshot index as character
     write(snap_index, ifilenameformat) ifile 
     snap_index = adjustl(snap_index) 
                       
     ! Show progress on post-processing    
     if (nrank==0) then
        print *,'----------------------------------------------------'

        ! Print the snapshot currently being processed         
        write(printing, '(A,A,A,A)') 'We are averaging snapshot = ', trim(snap_index), ' / ', trim(snap_n_index) 
        printing = adjustl(printing) 
        write(*,*) printing
                 
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
         write(filename, '(A,A,A)') 'ux-', trim(snap_index), '.bin' 
         filename = adjustl(filename)         
         call decomp_2d_read_one(1,ux1,dirname,filename,a)
        
         write(filename, '(A,A,A)') 'uy-', trim(snap_index), '.bin' 
         filename = adjustl(filename)         
         call decomp_2d_read_one(1,uy1,dirname,filename,a)
  
         write(filename, '(A,A,A)') 'uz-', trim(snap_index), '.bin' 
         filename = adjustl(filename)         
         call decomp_2d_read_one(1,uz1,dirname,filename,a)
        
     endif
     
     if (read_pre) then
         write(filename, '(A,A,A)') 'pp-', trim(snap_index), '.bin' 
         filename = adjustl(filename)         
         call decomp_2d_read_one(1,pre1,dirname,filename,a)     
     endif
     
     if (read_phi) then   
         write(filename, '(A,A,A)') 'phi01-', trim(snap_index), '.bin' 
         filename = adjustl(filename)         
         call decomp_2d_read_one(1,phi1,dirname,filename,a)       
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
                                                                          
     if (post_vort) call stat_vorticity(ux1,uy1,uz1,phi1,nr,nt,vortxmean,vortymean,vortzmean,mean_gradientp,mean_gradientx,mean_gradientz,mean_gradphi)
     
     if (post_diss) call stat_dissipation(ux1,uy1,uz1,nr,nt,epsmean)

  enddo ! closing of the do-loop on the different flow realizations

#ifndef TTBL_MODE 
   
   !--- Correlation for channel, mean statistics must be calculated in a previous post-processing run ---!
   if(post_corz) then
   
   ! Fluctuations calculation
   do k=1,ysize(3)
       do i=1,ysize(1)
           do j=1,ysize(2)
               ux2(i,j,k) = ux2(i,j,k)-u1meanHT(j)
               uy2(i,j,k) = uy2(i,j,k)-v1meanHT(j)
               uz2(i,j,k) = uz2(i,j,k)-w1meanHT(j)
           enddo
       enddo
   enddo
   
   ! Correlation functions calculation (each subdomain, z-pencils)
   call stat_correlation_z(ux2,uy2,uz2,nx,nz,nr,nt,RuuzH1,RvvzH1,RwwzH1,RuvzH1)
      
   end if
   !-----------------------------------------------------------------------------------------------------!
   
  ! Closing of the do-loop for the different time units (or snapshots) (ie index)
  enddo 
#endif
  
!--- Mean over homogeneous directions (H = Homogeneous) ---!

  ! Denominator of the divisions
  den = real(nx*nz,mytype)

  ! Summation over x and z directions
  if (post_mean) then
     do k=1,ysize(3)
        do i=1,ysize(1)
           do j=1,ysize(2)          
              u1meanH1(j)=u1meanH1(j)+u1mean(i,j,k)/den
              v1meanH1(j)=v1meanH1(j)+v1mean(i,j,k)/den
              w1meanH1(j)=w1meanH1(j)+w1mean(i,j,k)/den
              u2meanH1(j)=u2meanH1(j)+u2mean(i,j,k)/den
              v2meanH1(j)=v2meanH1(j)+v2mean(i,j,k)/den
              w2meanH1(j)=w2meanH1(j)+w2mean(i,j,k)/den
              u3meanH1(j)=u3meanH1(j)+u3mean(i,j,k)/den
              v3meanH1(j)=v3meanH1(j)+v3mean(i,j,k)/den
              w3meanH1(j)=w3meanH1(j)+w3mean(i,j,k)/den
              u4meanH1(j)=u4meanH1(j)+u4mean(i,j,k)/den
              v4meanH1(j)=v4meanH1(j)+v4mean(i,j,k)/den
              w4meanH1(j)=w4meanH1(j)+w4mean(i,j,k)/den
              uvmeanH1(j)=uvmeanH1(j)+uvmean(i,j,k)/den
              uwmeanH1(j)=uwmeanH1(j)+uwmean(i,j,k)/den
              vwmeanH1(j)=vwmeanH1(j)+vwmean(i,j,k)/den
              pre1meanH1(j)=pre1meanH1(j)+pre1mean(i,j,k)/den
              pre2meanH1(j)=pre2meanH1(j)+pre2mean(i,j,k)/den                                                   
              phi1meanH1(j)=phi1meanH1(j)+phi1mean(i,j,k)/den
              phi2meanH1(j)=phi2meanH1(j)+phi2mean(i,j,k)/den
              uphimeanH1(j)=uphimeanH1(j)+uphimean(i,j,k)/den
              vphimeanH1(j)=vphimeanH1(j)+vphimean(i,j,k)/den
              wphimeanH1(j)=wphimeanH1(j)+wphimean(i,j,k)/den                                  
           enddo          
        enddo
     enddo
  endif
  
  if (post_vort) then
     do k=1,ysize(3)
        do i=1,ysize(1)
           do j=1,ysize(2)          
              vortxmeanH1(j)=vortxmeanH1(j)+vortxmean(i,j,k)/den
              vortymeanH1(j)=vortymeanH1(j)+vortymean(i,j,k)/den
              vortzmeanH1(j)=vortzmeanH1(j)+vortzmean(i,j,k)/den 
              mean_gradientpH1(j)=mean_gradientpH1(j)+mean_gradientp(i,j,k)/den
              mean_gradientxH1(j)=mean_gradientxH1(j)+mean_gradientx(i,j,k)/den 
              mean_gradientzH1(j)=mean_gradientzH1(j)+mean_gradientz(i,j,k)/den 
              mean_gradphiH1(j)=mean_gradphiH1(j)+mean_gradphi(i,j,k)/den                   
           enddo
        enddo
     enddo
  endif
  
  if (post_diss) then
     do k=1,ysize(3)
        do i=1,ysize(1)
           do j=1,ysize(2)          
              epsmeanH1(j)=epsmeanH1(j)+epsmean(i,j,k)/den                  
           enddo
        enddo
     enddo
  endif

#ifdef TTBL_MODE
   ! Reset to zero the arrays used to collect the averages locally
   call reset_averages()  
#endif
 
!-------- Mean over all MPI processes (T = Total) ---------!

  if (post_mean) then
     call MPI_ALLREDUCE(u1meanH1,u1meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(v1meanH1,v1meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(w1meanH1,w1meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(u2meanH1,u2meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(v2meanH1,v2meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(w2meanH1,w2meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(u3meanH1,u3meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(v3meanH1,v3meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(w3meanH1,w3meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(u4meanH1,u4meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(v4meanH1,v4meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(w4meanH1,w4meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(uvmeanH1,uvmeanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(uvmeanH1,uvmeanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(vwmeanH1,vwmeanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(pre1meanH1,pre1meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(pre2meanH1,pre2meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)   
     call MPI_ALLREDUCE(phi1meanH1,phi1meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(phi2meanH1,phi2meanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(uphimeanH1,uphimeanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(vphimeanH1,vphimeanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(wphimeanH1,wphimeanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)    
  endif
  
  if (post_vort) then
     call MPI_ALLREDUCE(vortxmeanH1,vortxmeanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(vortymeanH1,vortymeanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(vortzmeanH1,vortzmeanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)  
     call MPI_ALLREDUCE(mean_gradientpH1,mean_gradientpHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code) 
     call MPI_ALLREDUCE(mean_gradientxH1,mean_gradientxHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(mean_gradientzH1,mean_gradientzHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
     call MPI_ALLREDUCE(mean_gradphiH1,mean_gradphiHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
  endif
  
  if (post_diss) then
     call MPI_ALLREDUCE(epsmeanH1,epsmeanHT,ysize(2),real_type,MPI_SUM,MPI_COMM_WORLD,code)
  end if

  if(post_corz) then

! Correlations calculation for TTBL
#ifdef TTBL_MODE  
   ! Fluctuations calculation
   do k=1,ysize(3)
       do i=1,ysize(1)
           do j=1,ysize(2)
               ux2(i,j,k) = ux2(i,j,k)-u1meanHT(j)
               uy2(i,j,k) = uy2(i,j,k)-v1meanHT(j)
               uz2(i,j,k) = uz2(i,j,k)-w1meanHT(j)
           enddo
       enddo
   enddo
   
   ! Correlation functions calculation (each subdomain, z-pencils)
   call stat_correlation_z(ux2,uy2,uz2,nx,nz,nr,nt,RuuzH1,RvvzH1,RwwzH1,RuvzH1)

#endif
   
   ! Summation over all MPI processes (valid for both TTBL and Channel)   
   call MPI_REDUCE(RuuzH1,RuuzHT,zsize(2)*zsize(3),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
   call MPI_REDUCE(RvvzH1,RvvzHT,zsize(2)*zsize(3),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
   call MPI_REDUCE(RwwzH1,RwwzHT,zsize(2)*zsize(3),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
   call MPI_REDUCE(RuvzH1,RuvzHT,zsize(2)*zsize(3),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
   
  end if

!------------- MPI process nrank = 0 at work --------------!

  ! High-order moments (variance, skewness, kurtosis)
  if(nrank.eq.0) then  
     
     if (post_mean) then
        do j=1,ysize(2)
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
           uphimeanHT(j)=uphimeanHT(j)-u1meanHT(j)*phi1meanHT(j)
           vphimeanHT(j)=vphimeanHT(j)-v1meanHT(j)*phi1meanHT(j)
           wphimeanHT(j)=wphimeanHT(j)-w1meanHT(j)*phi1meanHT(j)         
        enddo
     endif

!------------------Write formatted data--------------------!
     
     ! New directory for the statistics
     write(dirname,"('data_post/')") 
        
     call system('mkdir -p '//trim(dirname))

     ! Formatted data
     print *,'----------------------------------------------------'
     write(*,"(1x,'Writing output data in formatted .txt file(s)')")
     
     ! Mean statistics writing
     if (post_mean) then

#ifdef TTBL_MODE  
        ! Writing the snapshot index as character
        write(snap_index, ifilenameformat) ifile 
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
        
        ! Header
        write(iunit, '(22(A13, A1, 1X))') 'mean[u]'  , ',', 'mean[v]'  , ',', 'mean[w]', ',', &
                                          'var[u]'   , ',', 'var[v]'   , ',', 'var[w]' , ',', &
                                          'skew[u]'  , ',', 'skew[v]'  , ',', 'skew[w]', ',', &
                                          'kurt[u]'  , ',', 'kurt[v]'  , ',', 'kurt[w]', ',', &
                                          "<u'v'>"   , ',', "<u'w'>"   , ',', "<v'w'>" , ',', &
                                          'mean[p]'  , ',', 'var[p]'   , ',',                 &
                                          'mean[phi]', ',', 'var[phi]' , ',',                 &
                                          "<u'phi'>" , ',', "<v'phi'>" , ',', "<w'phi'>" 
               
        do j = 1, ysize(2) 
       
            write(iunit, '(22(F13.9, A1, 1X))') u1meanHT(j),   ',', &
                                                v1meanHT(j),   ',', &       
                                                w1meanHT(j),   ',', &
                                                u2meanHT(j),   ',', &
                                                v2meanHT(j),   ',', &
                                                w2meanHT(j),   ',', &
                                                u3meanHT(j),   ',', &
                                                v3meanHT(j),   ',', &
                                                w3meanHT(j),   ',', &
                                                u4meanHT(j),   ',', &
                                                v4meanHT(j),   ',', &
                                                w4meanHT(j),   ',', &
                                                uvmeanHT(j),   ',', &
                                                uwmeanHT(j),   ',', &  
                                                vwmeanHT(j),   ',', &                                              
                                                pre1meanHT(j), ',', &
                                                pre2meanHT(j), ',', &                     
                                                phi1meanHT(j), ',', &
                                                phi2meanHT(j), ',', &                        
                                                uphimeanHT(j), ',', &
                                                vphimeanHT(j), ',', &
                                                wphimeanHT(j)

        end do
                               
        close(iunit)
     endif
     
     ! Vorticity mean statistics and mean gradients writing
     if (post_vort) then

#ifdef TTBL_MODE  
        ! Writing the snapshot index as character
        write(snap_index, ifilenameformat) ifile 
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
        
        ! Header
        write(iunit, '(7(A13, A1, 1X))') 'mean[omega_x]', ',', 'mean[omega_y]', ',', 'mean[omega_z]', ',', &
                                         'dU_par/dy'    , ',', 'dU/dy'        , ',', 'dW/dy',         ',', &
                                         'dPhi/dy'   
        
        do j = 1, ysize(2) 
                
            write(iunit, '(7(F13.9, A1, 1X))') vortxmeanHT(j),      ',',  &
                                               vortymeanHT(j),      ',',  &       
                                               vortzmeanHT(j),      ',',  &
                                               mean_gradientpHT(j), ',',  &
                                               mean_gradientxHT(j), ',',  &
                                               mean_gradientzHT(j), ',',  &
                                               mean_gradphiHT(j)      
        
        end do
                               
        close(iunit)
     endif
     
     ! Mean dissipation writing
     if (post_diss) then

#ifdef TTBL_MODE  
        ! Writing the snapshot index as character
        write(snap_index, ifilenameformat) ifile 
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
        
        ! Header
        write(iunit, '(1(A13, A1, 1X))') 'mean[eps]'
        
        do j = 1, ysize(2) 
      
            write(iunit, '(1(F13.9, A1, 1X))') epsmeanHT(j)
               
        end do
                               
        close(iunit)
     endif
          
     ! Correlation functions along z writing
     if (post_corz) then
     
!--- Streamwise autocorrelation function, Ruuz ---!     
#ifdef TTBL_MODE
        ! Writing the snapshot index as character
        write(snap_index, ifilenameformat) ifile 
        snap_index = adjustl(snap_index) 
        
        ! Write the Ruuz filename for TTBL
        write(filename, '(A,A,A)') 'Ruuz-', trim(snap_index), '.txt'
        filename = adjustl(filename)
#else
        ! Write the Ruuz filename for channel flow
        write(filename, '(A)') 'Ruuz.txt'
        filename = adjustl(filename)
#endif       

        ! Construct the format string
        write(format_string, '(A, I0, A)') '(', nz, '(F13.9, A1, 1X))'

        ! Open the file and write      
        open(newunit=iunit,file=trim(dirname)//trim(filename),form='formatted')
        
        ! Streamwise fluctuations correlation function        
        do j = 1, ysize(2)
                
            write(iunit, format_string) (RuuzHT(j, k), ' ', k = 1, zsize(3))
               
        end do
                                       
        close(iunit)

!--- Vertical autocorrelation function, Rvvz ---!
#ifdef TTBL_MODE
        ! Writing the snapshot index as character
        write(snap_index, ifilenameformat) ifile 
        snap_index = adjustl(snap_index) 
        
        ! Write the Rvvz filename for TTBL
        write(filename, '(A,A,A)') 'Rvvz-', trim(snap_index), '.txt'
        filename = adjustl(filename)
#else
        ! Write the Rvvz filename for channel flow
        write(filename, '(A)') 'Rvvz.txt'
        filename = adjustl(filename)
#endif       

        ! Construct the format string
        write(format_string, '(A, I0, A)') '(', nz, '(F13.9, A1, 1X))'

        ! Open the file and write      
        open(newunit=iunit,file=trim(dirname)//trim(filename),form='formatted')
        
        ! Vertical fluctuations correlation function        
        do j = 1, ysize(2)
                
            write(iunit, format_string) (RvvzHT(j, k), ' ', k = 1, zsize(3))
               
        end do
                                       
        close(iunit)

!--- Spanwise autocorrelation function, Rwwz ---!
#ifdef TTBL_MODE
        ! Writing the snapshot index as character
        write(snap_index, ifilenameformat) ifile 
        snap_index = adjustl(snap_index) 
        
        ! Write the Rwwz filename for TTBL
        write(filename, '(A,A,A)') 'Rwwz-', trim(snap_index), '.txt'
        filename = adjustl(filename)
#else
        ! Write the Rwwz filename for channel flow
        write(filename, '(A)') 'Rwwz.txt'
        filename = adjustl(filename)
#endif       

        ! Construct the format string
        write(format_string, '(A, I0, A)') '(', nz, '(F13.9, A1, 1X))'

        ! Open the file and write      
        open(newunit=iunit,file=trim(dirname)//trim(filename),form='formatted')
        
        ! Vertical fluctuations correlation function        
        do j = 1, ysize(2)
                
            write(iunit, format_string) (RwwzHT(j, k), ' ', k = 1, zsize(3))
               
        end do
                                       
        close(iunit)

!--- Mixed fluctuations correlation function, Ruvz ---!
#ifdef TTBL_MODE
        ! Writing the snapshot index as character
        write(snap_index, ifilenameformat) ifile 
        snap_index = adjustl(snap_index) 
        
        ! Write the Ruvz filename for TTBL
        write(filename, '(A,A,A)') 'Ruvz-', trim(snap_index), '.txt'
        filename = adjustl(filename)
#else
        ! Write the Ruvz filename for channel flow
        write(filename, '(A)') 'Ruvz.txt'
        filename = adjustl(filename)
#endif       

        ! Construct the format string
        write(format_string, '(A, I0, A)') '(', nz, '(F13.9, A1, 1X))'

        ! Open the file and write      
        open(newunit=iunit,file=trim(dirname)//trim(filename),form='formatted')
        
        ! Vertical fluctuations correlation function        
        do j = 1, ysize(2)
                
            write(iunit, format_string) (RuvzHT(j, k), ' ', k = 1, zsize(3))
               
        end do
                                       
        close(iunit)
         
     endif ! closing of if-statement for writing correlations
                     
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
     print *,'mean[omega_x], mean[omega_y], mean[omega_z]'
     print *,'dU_par/dy,     dU/dy,         dW/dy'
     print *,'dPhi/dy'
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
     
     ! Correlation functions 
     if (post_corz) then
     print *,'==========================================================='
     print *,''
     print *,'The following statistics have been saved in'
     print *,'"Riiz and Ruvz" file(s):'
     print *,''
     print *,'Ruu(z), Rvv(z), Rww(z), Ruv(z)'
     print *,''    
     endif
     
     print *,'==========================================================='
     
  endif

  call decomp_2d_finalize
  call MPI_FINALIZE(code)

end program post





