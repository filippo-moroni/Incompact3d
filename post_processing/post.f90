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

  implicit none

  integer :: code,i,j,k,is
  integer :: ii,ie                                 ! internal and external loops 
  integer :: file1,filen,icrfile,nt
  integer :: nr                                    ! total number of flow realizations
  integer :: ifile,ssfile,iwrn,num
  
  real(mytype) :: tstart,t1,trank,tranksum,ttotal,tremaining,telapsed,trstart, trend
  
  ! Added by R. Corsini
  integer :: ix,rr,ipr,jpr,kpr,ierror,ipos  
  integer(8) :: ttsize 
  integer,dimension(2) :: ipost                    ! index for the number of post-processing subroutines employed
  logical :: read_phi,read_vel,read_ibm,read_pre  
  real(mytype) :: xpos 
  character(30) :: filename,dirname 
  character(1) :: a

  TYPE(DECOMP_INFO) :: phG,ph1,ph2,ph3,ph4

  call ft_parameter(.true.)

  CALL MPI_INIT(code)
  
  ! Modified by R. Corsini
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror) 
  call decomp_2d_init(nx,ny,nz,1,nproc) 
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)     !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)     !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.)  !start from 1 == true
  call parameter()

  ! Start of the post-processing
  p_row=1; p_col=nproc
  
  post_mean=.false.; post_vort=.false.  
  read_vel=.false.; read_pre=.false.; read_phi=.false.; read_ibm=.false.
  
  ! Reading of the input file for post-processing
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
  do i=1,size(ipost)
     read (10,*) ipost(i)
  enddo
  close(10)

  ! Total number of Snapshots in time
  nt = (filen-file1)/icrfile+1
  
  if (ipost(1)==1) post_mean=.true.
  if (ipost(2)==1) post_vort=.true.

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

  ! Initialize variables and statistics
  call init_post_variables
  call schemes
  call init_statistics


  ttsize=0
  if (read_vel) ttsize=ttsize+3
  if (read_pre) ttsize=ttsize+1
  if (read_phi) ttsize=ttsize+nphi
  ttsize=ttsize*nx*ny*nz
  tstart=0.;t1=0.;trank=0.;tranksum=0.;ttotal=0.
  call cpu_time(tstart)

!-------An extra cycle for different time units is thus required, external do loop

 do ie=1,nt
 
     call cpu_time(t1)
     ifile = (ii-1)*icrfile+file1
     t=dt*real(imodulo*ifile,mytype)
     itime=imodulo*ifile
    
     if (nrank==0) then
        print *,'----------------------------------------------------'
        write(*,"(' We are averaging the realizations of the SnapShot =',I3,'/',I3,', Time unit =',F9.5)") ifile,filen,t
     endif
 
!---------Start of the ensemble average cycle--------------!

  do ii=1,nr 
                                    
     ! Read data
     call cpu_time(trstart)
     
     if (read_vel) then
        write(filename,"('./data/ux',I3.3,'_',I1.1)") ifile, nr
        call decomp_2d_read_one(1,ux1,filename)       
        write(filename,"('./data/uy',I3.3,'_',I1.1)") ifile, nr
        call decomp_2d_read_one(1,uy1,filename)       
        write(filename,"('./data/uz',I3.3,'_',I1.1)") ifile, nr
        call decomp_2d_read_one(1,uz1,filename)
        call test_speed_min_max(ux1,uy1,uz1)
     endif
     
     if (read_pre) then
        write(filename,"('./data/pp',I3.3,'_',I1.1)") ifile,nr
        call decomp_2d_read_one(1,pre1,filename)
        if (nscheme==2) then
            pre1 = pre1/dt  !IF nscheme = 2
        else
            if (nrank==0) print *, '!!WARNING!! Pressure field not scaled by dt'
        endif
     endif
     
     if (read_phi) then
        do is=1, nphi
           write(filename,"('./data/phi',I1.1,I4.4,'_',I1.1)") is, ifile, nr
           call decomp_2d_read_one(1,phi1(:,:,:,is),filename)
        enddo
     endif
     
     call cpu_time(trend)

     ! Start statistics computation
     if (post_mean) call STAT_MEAN(ux1,uy1,uz1,pre1,phi1,ta1, &
                                   u1mean,v1mean,w1mean,u2mean,v2mean,w2mean, &
                                   u3mean,v3mean,w3mean,u4mean,v4mean,w4mean, &
                                   uvmean,uwmean,vwmean,pre1mean,pre2mean,phi1mean, &
                                   phi2mean,uphimean,vphimean,wphimean,nr)

     if (post_vort) call STAT_VORTICITY(ux1,uy1,uz1,ifile)


  enddo ! closing of the do-loop on the different flow realizations

!----Mean over homogeneous directions (H = Homogeneous)----!

  ! Summation over x and z directions
  if (post_mean) then
     do k=xstart(3),xend(3)
      do i=xstart(1),xend(1)
        do j=xstart(2),xend(2)          
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
              
              if (iscalar==1) then
                 do is=1,nphi
                    phi1meanH1(j,is)=phi1meanH1(j,is)+phi1mean(i,j,k,is)/real(nx*nz,mytype)
                    phi2meanH1(j,is)=phi2meanH1(j,is)+phi2mean(i,j,k,is)/real(nx*nz,mytype)
                    uphimeanH1(j,is)=uphimeanH1(j,is)+uphimean(i,j,k,is)/real(nx*nz,mytype)
                    vphimeanH1(j,is)=vphimeanH1(j,is)+vphimean(i,j,k,is)/real(nx*nz,mytype)
                    wphimeanH1(j,is)=wphimeanH1(j,is)+wphimean(i,j,k,is)/real(nx*nz,mytype)
                 enddo
              endif 
           enddo
        enddo
     enddo
  endif

!---------Mean over all MPI processes (T = Total)----------!

  if (post_mean) then
     call MPI_REDUCE(u1meanH1,u1meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(v1meanH1,v1meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(w1meanH1,w1meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(u2meanH1,u2meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(v2meanH1,v2meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(w2meanH1,w2meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(u3meanH1,u3meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(v3meanH1,v3meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(w3meanH1,w3meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(u4meanH1,u4meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(v4meanH1,v4meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(w4meanH1,w4meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(uvmeanH1,uvmeanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(uvmeanH1,uvmeanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(vwmeanH1,vwmeanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(pre1meanH1,pre1meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     call MPI_REDUCE(pre2meanH1,pre2meanHT,xsize(2),real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    
     if (iscalar==1) then
        call MPI_REDUCE(phi1meanH1,phi1meanHT,xsize(2)*nphi,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
        call MPI_REDUCE(phi2meanH1,phi2meanHT,xsize(2)*nphi,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
        call MPI_REDUCE(uphimeanH1,uphimeanHT,xsize(2)*nphi,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
        call MPI_REDUCE(vphimeanH1,vphimeanHT,xsize(2)*nphi,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
        call MPI_REDUCE(wphimeanH1,wphimeanHT,xsize(2)*nphi,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
     endif
  endif

!--------------MPI process nrank = 0 at work---------------!

  if(nrank.eq.0) then ! only processor 0 is working
     if (post_mean) then
        do j=xstart(2),xend(2)
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
           
           if (iscalar==1) then
              do is=1,nphi
                 phi2meanHT(j,is)=phi2meanHT(j,is)-phi1meanHT(j,is)**2
                 uphimeanHT(j,is)=uphimeanHT(j,is)-u1meanHT(j)*phi1meanH1(j,is)
                 vphimeanHT(j,is)=vphimeanHT(j,is)-v1meanHT(j)*phi1meanH1(j,is)
                 wphimeanHT(j,is)=wphimeanHT(j,is)-w1meanHT(j)*phi1meanH1(j,is)
              enddo
           endif
        enddo
     endif

!------------------Write unformatted data------------------!
     
     write(dirname,"('data_post_TU',I2.1,'/')") t
     
     !write(dirname,*) 'data_post/'
     
     call system('mkdir -p '//trim(dirname))

     ! Unformatted data
     print *,'----------------------------------------------------'
     write(*,"(1x,'Writing output data to binary files :')")
     
     if (post_mean) then
        write(filename,"('MEAN.bin')")
        write(*,"(3x,A)") filename
        open(unit=210,file=trim(dirname)//trim(filename), &
             form='unformatted',access='stream',status='replace')
        write(210,pos=1) u1meanHT,v1meanHT,w1meanHT, &
                         u2meanHT,v2meanHT,w2meanHT, &
                         u3meanHT,v3meanHT,w3meanHT, &
                         u4meanHT,v4meanHT,w4meanHT, &
                         uvmeanHT,uwmeanHT,vwmeanHT, &
                         pre1meanHT,pre2meanHT
        if (iscalar==1) then
           do is=1,nphi
              inquire(unit=210,pos=ipos)
              write(210,pos=ipos) phi1meanHT(:,is),phi2meanHT(:,is), &
                                  uphimeanHT(:,is),vphimeanHT(:,is),wphimeanHT(:,is)
           enddo
        endif
        close(210)
     endif

  endif ! closing of the if-statement for processor 0
  
 enddo  ! closing of the do-loop for the different time units (or SnapShots)
  
  ! End of post-processing

  call cpu_time(trank)
  ttotal=trank-tstart

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
  endif

  call decomp_2d_finalize
  CALL MPI_FINALIZE(code)
end PROGRAM post



