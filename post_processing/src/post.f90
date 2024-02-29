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
  integer :: ifile
  
  real(mytype) :: tstart=0.0,tend=0.0,ttotal=0.0   ! variables to count time spent to post-process data
   
  integer  :: iunit				   ! unit for the file to open (assigned by the compiler)
    
  integer,dimension(2) :: sel                      ! index for the number of post-processing subroutines employed (selector index)
  logical :: read_vel,read_pre,read_phi  
 
  character(99):: filename,dirname 
  character(1) :: a
  
  ! Integer for MPI
  integer :: code

  ! Variables to read the input.i3d file
  integer :: nargin, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase
    
  ! Initialize MPI
  CALL MPI_INIT(code)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,code) 
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,code)  ! same as R. Corsini & Xcompact3d
  
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
  post_mean=.false.; post_vort=.false.  
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
  
  ! Total number of Snapshots in time
  nt = (filen-file1)/icrfile+1
  
  ! Initialize statistics
  call init_statistics()

  ! Time passed since the program has started
  call cpu_time(tstart)
    
!-----------------------------!
! Post-processing starts here !
!-----------------------------!

!------------Start of the time unit do-loop----------------!

 do ie=1,nt
      
     ! Number of snapshot  
     ifile = (ie-1)*icrfile+file1
     
     ! Time-unit (T)
     t=dt*real(ioutput*ifile,mytype)
     
     ! Show progress on post-processing    
     if (nrank==0) then
        print *,'----------------------------------------------------'
        write(*,"(' We are averaging the realizations of the SnapShot =',I3,'/',I3,', Time unit =',F8.3)") ifile,filen,t
     endif
 
!---------Start of the ensemble average cycle--------------!

  do ii=1,nr 
                                         
     ! Writing the directory where snapshots are saved (different realizations)
     write(dirname,"('./data_r',I1.1)") nr
     
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
     if (post_mean) call STAT_MEAN(ux2,uy2,uz2,pre2,phi2,ta2, &
                                   u1mean,v1mean,w1mean,u2mean,v2mean,w2mean, &
                                   u3mean,v3mean,w3mean,u4mean,v4mean,w4mean, &
                                   uvmean,uwmean,vwmean,pre1mean,pre2mean,phi1mean, &
                                   phi2mean,uphimean,vphimean,wphimean,nr)
                                       
     if (post_vort) call STAT_VORTICITY(ux2,uy2,uz2,ifile,nr)

  enddo ! closing of the do-loop on the different flow realizations

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

!--------------MPI process nrank = 0 at work---------------!

  if(nrank.eq.0) then ! only processor 0 is working
     
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

!------------------Write unformatted data------------------!
     
     ! New directory for the statistics
     write(dirname,"('data_post/')") 
        
     call system('mkdir -p '//trim(dirname))

     ! Formatted data
     print *,'----------------------------------------------------'
     write(*,"(1x,'Writing output data in formatted .txt files :')")
     
     if (post_mean) then
     
        write(filename,"('mean_statistics',F5.1,'.txt')") t
        write(*,"(3x,A)") filename
        
        filename = adjustl(filename)
        
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
  
!------------Calculate parameters of the flow--------------!
  
     ! Reading the coordinates in y of faces' elements
     write(filename,"('yp.dat')") 
               
     filename = adjustl(filename)
        
     open(newunit=iunit,file=trim(filename),form='formatted')
     
     do j = ystart(2),yend(2)
     
        read(iunit, *) yp(j)
    
     end do
     
     close(iunit)
     
     ! delta_99
     do j = ystart(2),yend(2)
     
        delta_99(ie) = yp(j)  
        
        if(u1meanHT(j) > 0.99*u1meanHT(yend(2))) exit
               
     end do
     
     ! displacement thickness
     do j = ystart(2),yend(2) - 1
     
     disp_t(ie) = disp_t(ie) + u1meanHT(j)*(yp(j+1) - yp(j))
                    
     end do
     
     disp_t(ie) = yp(ysize(2)) - disp_t(ie)
     
     ! momentum thickness
     do j = ystart(2),yend(2) - 1
     
     mom_t(ie) = mom_t(ie) + (u1meanHT(j) - u1meanHT(j)**2)*(yp(j+1) - yp(j))
                    
     end do
     
     ! friction or shear velocity
     sh_vel(ie) = sqrt(xnu*u1meanHT(2)/yp(2))
     
     ! Reynolds numbers
     re_tau  (ie) = delta_99(ie)*sh_vel(ie)/xnu  ! friction Re number (or delta99^+)
     re_ds   (ie) = disp_t  (ie)*sh_vel(ie)/xnu  ! Re number based on displacement thickness delta star (ds)
     re_theta(ie) = mom_t   (ie)*sh_vel(ie)/xnu  ! Re number based on momentum thickness theta     
     
  endif ! closing of the if-statement for processor 0
  
 enddo  ! closing of the do-loop for the different time units (or SnapShots) (ie index)
 
 if(nrank.eq.0) then
 
 ! New directory for the statistics
     write(dirname,"('data_post_parameters/')") 
        
     call system('mkdir -p '//trim(dirname))
     
     write(filename,"('parameters_statistics.txt')") 
         
        filename = adjustl(filename)
        
        open(newunit=iunit,file=trim(dirname)//trim(filename),form='formatted')
        
        do ie = 1,nt + 1
        
        ! Number of snapshot  
        ifile = (ie-1)*icrfile+file1
     
        ! Time-unit (T)
        t=dt*real(ioutput*ifile,mytype)
        
        if (ie .eq. 1) then
        
        write (iunit, *) 'delta_99'      , ',','disp_thickness', ',','mom_thickness', ',', &
                         're_tau'        , ',','re_delta*'     , ',','re_theta'            &
                         'shear_velocity', ',','time unit'                         
                               
        else
        
        write(iunit, *)  delta_99(ie-1), ',', &
                         disp_t  (ie-1), ',', &
                         mom_t   (ie-1), ',', &
                         re_tau  (ie-1), ',', &
                         re_ds   (ie-1), ',', &
                         re_theta(ie-1), ',', &
                         sh_vel  (ie-1), ',', &
                         t
        
        end if
        
        end do
 
 end if
   
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
     
     ! Averages 
     if (post_mean) then
     print *,'==========================================================='
     print *,''
     print *,'The following statistics have been saved in "MEAN" files:'
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
     
  endif

  call decomp_2d_finalize
  CALL MPI_FINALIZE(code)

end PROGRAM post




