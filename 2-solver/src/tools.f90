!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module tools

  implicit none

  logical, save :: adios2_restart_initialised = .false.

  character(len=*), parameter :: io_restart = "restart-io"
  character(len=*), parameter :: resfile = "checkpoint"
  
  private

  public :: program_header,        &
            test_speed_min_max,    &
            test_scalar_min_max,   &
            simu_stats,            &
            restart,               &
            apply_spatial_filter,  &
            compute_cfldiff,       &
            compute_cfl,           &
            compute_reynolds_cell, &
            compute_stab_param,    &
            update_time_int_coeff, &
            rescale_pressure,      &
            mean_plane_x,          &
            mean_plane_y,          &
            mean_plane_z,          &
            avg3d

contains

  !----------------------------------------------!
  ! Header of the program printed to the screen. !
  !----------------------------------------------!
  subroutine program_header()
  
  use decomp_2d, only : nrank
  
  implicit none
  
  if (nrank==0) then
     write(*,*) '!---------------------------------------------------------!'
     write(*,*) '!                   ~  Incompact3D  ~                     !'
     write(*,*) '!  Copyright (c) 2018 Eric Lamballais and Sylvain Laizet  !'
     write(*,*) '!  Modified by Felipe Schuch and Ricardo Frantz           !'
     write(*,*) '!  Modified by Paul Bartholomew, Georgios Deskos and      !'
     write(*,*) '!  Sylvain Laizet, 2018                                   !'
     write(*,*) '!                                                         !'
     write(*,*) '!  Modified by Filippo Moroni, 2024                       !'
     write(*,*) '!---------------------------------------------------------!'
     
#if defined(VERSION)
     write(*,*)'Git version        : ', VERSION
#else
     write(*,*)'Git version        : unknown'
#endif
  endif
  
  end subroutine program_header
  !-----------------------------------------------------------------------------!
  subroutine test_speed_min_max(ux,uy,uz)

   use decomp_2d
   use variables
   use param
   use var
   use mpi

   implicit none

   integer :: code,ierror,i,j,k
   real(mytype) :: uxmax,uymax,uzmax,uxmin,uymin,uzmin
   real(mytype) :: uxmax1,uymax1,uzmax1,uxmin1,uymin1,uzmin1
   real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
   real(mytype),dimension(6) :: umaxin, umaxout

   if (iibm > 0) then
      ux(:,:,:) = (one - ep1(:,:,:)) * ux(:,:,:)
      uy(:,:,:) = (one - ep1(:,:,:)) * uy(:,:,:)
      uz(:,:,:) = (one - ep1(:,:,:)) * uz(:,:,:)
   endif

   uxmax=-1609.;uymax=-1609.;uzmax=-1609.;uxmin=1609.;uymin=1609.;uzmin=1609.
   
   ! More efficient version
   uxmax=maxval(ux)
   uymax=maxval(uy)
   uzmax=maxval(uz)
   uxmin=-minval(ux)
   uymin=-minval(uy)
   uzmin=-minval(uz)

   umaxin = (/uxmax, uymax, uzmax, uxmin, uymin, uzmin/)
   call MPI_REDUCE(umaxin,umaxout,6,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)

   uxmax1= umaxout(1)
   uymax1= umaxout(2)
   uzmax1= umaxout(3)
   uxmin1=-umaxout(4)
   uymin1=-umaxout(5)
   uzmin1=-umaxout(6)

   if (nrank == 0) then

      write(*,*) 'U,V,W min=',real(uxmin1,4),real(uymin1,4),real(uzmin1,4)
      write(*,*) 'U,V,W max=',real(uxmax1,4),real(uymax1,4),real(uzmax1,4)
      !print *,'CFL=',real(abs(max(uxmax1,uymax1,uzmax1)*dt)/min(dx,dy,dz),4)

      if((abs(uxmax1)>=onehundred).or.(abs(uymax1)>=onehundred).OR.(abs(uzmax1)>=onehundred)) then
        write(*,*) 'Velocity diverged! SIMULATION IS STOPPED!'
        call MPI_ABORT(MPI_COMM_WORLD,code,ierror)
        stop
      endif

   endif

   return
  end subroutine test_speed_min_max
 !-----------------------------------------------------------------------------!
  subroutine test_scalar_min_max(phi)

    use decomp_2d
    use variables
    use param
    use var
    use mpi

    implicit none

    integer :: code,ierror,i,j,k,is,jglob
    real(mytype) :: phimax,phimin,phimax1,phimin1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(2,numscalar) :: phimaxin,phimaxout

    do is=1, numscalar

      ta1(:,:,:) = phi(:,:,:,is)
      ! ibm
      if (iibm > 0) then
        ta1(:,:,:) = (one - ep1(:,:,:)) * ta1(:,:,:)
      endif

      phimax=-1609._mytype
      phimin=1609._mytype
      phimax = maxval(ta1(:,:,:))
      phimin =-minval(ta1(:,:,:))
      phimaxin(:,is) =  (/phimin, phimax /)
    enddo

    call MPI_REDUCE(phimaxin,phimaxout,numscalar*2,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)

    do is=1,numscalar
      if (nrank == 0) then
        phimin1 = -phimaxout(1,is)
        phimax1 =  phimaxout(2,is)

        write(*,*) 'Phi'//char(48+is)//' min max=', real(phimin1,4), real(phimax1,4)

        if (abs(phimax1) > 100._mytype) then !if phi control turned off
           write(*,*) 'Scalar diverged! SIMULATION IS STOPPED!'
           call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
        endif
      endif
    enddo

    return
  end subroutine test_scalar_min_max
  !-----------------------------------------------------------------------------!
  subroutine simu_stats(iwhen)

    use decomp_2d
    use simulation_stats
    use var
    use MPI

    implicit none

    integer :: iwhen

    ! At the start of the simulation
    if (iwhen == 1) then 
       tstart=zero
       time1=zero
       trank=zero
       tranksum=zero
       ttotal=zero
       call cpu_time(tstart)

    ! At the start of a time step  
    else if (iwhen == 2) then 
       if (nrank == 0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime==ilast)) then
          call cpu_time(time1)
          write(*,*) '==========================================================='
          write(*,"(' Time step =',i7,'/',i7,', Time unit =',F12.4)") itime,ilast,t
       endif

    ! At the end of a time step
    else if ((iwhen == 3).and.(itime > ifirst)) then 
       if (nrank == 0.and.(mod(itime, ilist) == 0 .or. itime == ifirst .or. itime==ilast)) then
          call cpu_time(trank)
          if (nrank==0) write(*,*) 'Time for this time step (s):',real(trank-time1)
          telapsed = (trank-tstart)/threethousandsixhundred
          tremaining  = telapsed*(ilast-itime)/(itime-ifirst)
          write(*,"(' Remaining time:',I8,' h ',I2,' min')") int(tremaining), int((tremaining-int(tremaining))*sixty)
          write(*,"(' Elapsed time:  ',I8,' h ',I2,' min')") int(telapsed), int((telapsed-int(telapsed))*sixty)
       endif

    ! At the end of the simulation
    else if (iwhen == 4) then 
       call cpu_time(trank)
       ttotal=trank-tstart
       if (nrank == 0) then
          write(*,*) '==========================================================='
          write(*,*) '                                                           '
          write(*,*) 'Good job! Xcompact3d finished successfully!                '
          write(*,*) '                                                           '
          write(*,*) '2DECOMP with p_row*p_col=',p_row,p_col
          write(*,*) '                                                           '
          write(*,*) 'nx*ny*nz=',nx*ny*nz
          write(*,*) 'nx,ny,nz=',nx,ny,nz
          write(*,*) 'dx,dy,dz=',dx,dy,dz
          write(*,*) '                                                           '
          write(*,*) 'Averaged time per step (s):',real(ttotal/(ilast-(ifirst-1)),4)
          write(*,*) 'Total wallclock (s):',real(ttotal,4)
          write(*,*) 'Total wallclock (m):',real(ttotal/sixty,4)
          write(*,*) 'Total wallclock (h):',real(ttotal/threethousandsixhundred,4)
          write(*,*) '                                                           '
       endif
    endif

  end subroutine simu_stats
  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: restart
  ! DESCRIPTION: reads or writes restart file
  !      AUTHOR: ?
  !    MODIFIED: Kay Schäfer
  !-----------------------------------------------------------------------------!
  subroutine restart(ux1,uy1,uz1,dux1,duy1,duz1,ep1,pp3,phi1,dphi1,px1,py1,pz1,rho1,drho1,mu1,iresflg)

    use decomp_2d
    use decomp_2d_io
    use variables
    use param
    use MPI
    use navier, only : gradp

    implicit none

    integer :: i,j,k,iresflg,is,code
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1
    real(mytype), dimension(phG%zst(1):phG%zen(1),phG%zst(2):phG%zen(2),phG%zst(3):phG%zen(3)) :: pp3
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3),ntime) :: drho1
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: mu1
    real(mytype) :: tfield

    logical :: fexists
    character(len=90) :: filestart
    character(len=7)  :: fmt1
    character(len=32) :: fmt2,fmt3,fmt4
    character(len=80) :: varname
    
    character(99) :: filename, ts_index
    integer       :: iunit
    logical       :: exists
  
    NAMELIST /Time/ tfield, itime
    NAMELIST /NumParam/ nx, ny, nz, istret, beta, dt, itimescheme

    write(filename,"('restart',I7.7)") itime
    write(filestart,"('restart',I7.7)") ifirst-1

    ! Writing restart (restart flag)
    if (iresflg == 1) then 
       if (mod(itime, icheckpoint) /= 0) then
          return
       endif

       if (nrank==0) then
          write(*,*) '===========================================================<<<<<'
          write(*,*) 'Writing restart point ',filename !itime/icheckpoint
       endif
    end if

    if (.not. adios2_restart_initialised) then
       call init_restart_adios2()
       adios2_restart_initialised = .true.
    end if

    ! Write in 'checkpoint' file
    if (iresflg==1) then 
               
       call decomp_2d_open_io(io_restart, resfile, decomp_2d_write_mode)
       call decomp_2d_start_io(io_restart, resfile)

       ! Write velocity field
       call decomp_2d_write_one(1,ux1,resfile,"ux",0,io_restart,reduce_prec=.false.)
       call decomp_2d_write_one(1,uy1,resfile,"uy",0,io_restart,reduce_prec=.false.)
       call decomp_2d_write_one(1,uz1,resfile,"uz",0,io_restart,reduce_prec=.false.)
       
       ! Write previous time-step if necessary for AB2 or AB3
       if ((itimescheme==2).or.(itimescheme==3)) then
          call decomp_2d_write_one(1,dux1(:,:,:,2),resfile,"dux-2",0,io_restart,reduce_prec=.false.)
          call decomp_2d_write_one(1,duy1(:,:,:,2),resfile,"duy-2",0,io_restart,reduce_prec=.false.)
          call decomp_2d_write_one(1,duz1(:,:,:,2),resfile,"duz-2",0,io_restart,reduce_prec=.false.)
       end if

       ! One more previous time-step for AB3
       if (itimescheme==3) then
          call decomp_2d_write_one(1,dux1(:,:,:,3),resfile,"dux-3",0,io_restart,reduce_prec=.false.)
          call decomp_2d_write_one(1,duy1(:,:,:,3),resfile,"duy-3",0,io_restart,reduce_prec=.false.)
          call decomp_2d_write_one(1,duz1(:,:,:,3),resfile,"duz-3",0,io_restart,reduce_prec=.false.)
       end if
       
       ! Pressure
       call decomp_2d_write_one(3,pp3,resfile,"pp",0,io_restart,phG,reduce_prec=.false.)
       
       ! Scalar
       if (iscalar==1) then
          do is=1, numscalar
             write(varname, *) "phi-", is
             call decomp_2d_write_one(1,phi1(:,:,:,is),resfile,varname,0,io_restart,reduce_prec=.false.)
             
             ! Previous time-step, AB2 or AB3
             if ((itimescheme==2).or.(itimescheme==3)) then 
                write(varname, *) "dphi-", is, "-2"
                call decomp_2d_write_one(1,dphi1(:,:,:,2,is),resfile,varname,0,io_restart,reduce_prec=.false.)
             end if
             
             ! Older time step, AB3
             if (itimescheme==3) then
               write(varname, *) "dphi-", is, "-3"
               call decomp_2d_write_one(1,dphi1(:,:,:,3,is),resfile,varname,0,io_restart,reduce_prec=.false.)
             end if
          end do
       endif
       
       !if (ilmn) then
       !   do is = 1, nrhotime
       !      write(varname, *) "rho-", is
       !      call decomp_2d_write_one(1,rho1(:,:,:,is),resfile,varname,0,io_restart,reduce_prec=.false.)
       !   enddo
       !   do is = 1, ntime
       !      write(varname, *) "drho-", is
       !      call decomp_2d_write_one(1,drho1(:,:,:,is),resfile,varname,0,io_restart,reduce_prec=.false.)
       !   enddo
       !   call decomp_2d_write_one(1,mu1(:,:,:),resfile,"mu",0,io_restart,reduce_prec=.false.)
       !endif

       call decomp_2d_end_io(io_restart, resfile)
       call decomp_2d_close_io(io_restart, resfile)
             
       !--- Store checkpoint in a file copy for safer backup ---!
       if(nrank.eq.0) then
      
           ! Writing the last time step index as character
           write(ts_index,'(I8.8)') itime
           ts_index = adjustl(ts_index) 
        
           ! Write the filename for the checkpoint file 
           write(filename, '(A,A)') 'checkpoint-', trim(ts_index)
           filename = adjustl(filename)
  
           ! Copy and store the checkpoint file with a different name
           call execute_command_line('cp ' // 'checkpoint' // ' ' // filename)
      
           ! Move the created file inside /checkpoints folder
           call execute_command_line('mv ' // filename // ' ' // 'data/checkpoints')
      
       end if
       
       ! Write info file for restart - Kay Schäfer
       if (nrank == 0) then
         write(filename,"('data/restart_info/restart',I7.7,'.info')") itime
         write(fmt2,'("(A,I16)")')
         write(fmt3,'("(A,F16.4)")')
         write(fmt4,'("(A,F16.12)")')
         
         ! Open and write
         open (111,file=filename,action='write',status='replace')
         write(111,'(A)')'!==========================='
         write(111,'(A)')'&Time'
         write(111,'(A)')'!==========================='
         write(111,fmt3) 'tfield=     ',t
         write(111,fmt2) 'itime=      ',itime
         write(111,'(A)')'!==========================='
         write(111,'(A)')'&NumParam'
         write(111,'(A)')'!==========================='
         write(111,fmt2) 'nx=         ',nx
         write(111,fmt2) 'ny=         ',ny
         write(111,fmt2) 'nz=         ',nz
         write(111,fmt3) 'Lx=         ',xlx
         write(111,fmt3) 'Ly=         ',yly
         write(111,fmt3) 'Lz=         ',zlz
         write(111,fmt2) 'istret=     ',istret
         write(111,fmt4) 'beta=       ',beta
         write(111,fmt2) 'numscalar=  ',numscalar
         write(111,fmt2) 'itimescheme=',itimescheme
         write(111,fmt2) 'iimplicit=  ',iimplicit
         write(111,'(A)')'!==========================='
         write(111,'(A)')'&NumStability'
         write(111,'(A)')'!==========================='
         write(111,fmt3) 'dt=         ',dt
         write(111,fmt3) 'CFL,max,sum=',cflmax        
         write(111,'(A)')'!==========================='
         close(111)
       end if
    else
       if (nrank==0) then
         write(*,*)'==========================================================='
         write(*,*)'RESTART from file:', filestart
         write(*,*)'==========================================================='
       end if
       call decomp_2d_open_io(io_restart, resfile, decomp_2d_read_mode)
       call decomp_2d_start_io(io_restart, resfile)

       call decomp_2d_read_one(1,ux1,resfile,"ux",io_restart,reduce_prec=.false.)
       call decomp_2d_read_one(1,uy1,resfile,"uy",io_restart,reduce_prec=.false.)
       call decomp_2d_read_one(1,uz1,resfile,"uz",io_restart,reduce_prec=.false.)
       
       ! Read previous time-step if necessary for AB2 or AB3
       if ((itimescheme==2).or.(itimescheme==3)) then
          call decomp_2d_read_one(1,dux1(:,:,:,2),resfile,"dux-2",io_restart,reduce_prec=.false.)
          call decomp_2d_read_one(1,duy1(:,:,:,2),resfile,"duy-2",io_restart,reduce_prec=.false.)
          call decomp_2d_read_one(1,duz1(:,:,:,2),resfile,"duz-2",io_restart,reduce_prec=.false.)
       end if
       ! Ror AB3 one more previous time-step
       if (itimescheme==3) then 
          call decomp_2d_read_one(1,dux1(:,:,:,3),resfile,"dux-3",io_restart,reduce_prec=.false.)
          call decomp_2d_read_one(1,duy1(:,:,:,3),resfile,"duy-3",io_restart,reduce_prec=.false.)
          call decomp_2d_read_one(1,duz1(:,:,:,3),resfile,"duz-3",io_restart,reduce_prec=.false.)
       end if
       
       call decomp_2d_read_one(3,pp3,resfile,"pp",io_restart,phG,reduce_prec=.false.)
       
       if (iscalar==1) then
         do is=1, numscalar
            write(varname, *) "phi-", is
            call decomp_2d_read_one(1,phi1(:,:,:,is),resfile,varname,io_restart,reduce_prec=.false.)
           ! Previous time-steps, AB2 or AB3
           if ((itimescheme==2).or.(itimescheme==3)) then
             write(varname, *) "dphi-", is, "-2"
             call decomp_2d_read_one(1,dphi1(:,:,:,2,is),resfile,varname,io_restart,reduce_prec=.false.)
           end if
           ! AB3
           if (itimescheme==3) then 
              write(varname, *) "dphi-", is, "-3"
              call decomp_2d_read_one(1,dphi1(:,:,:,3,is),resfile,varname,io_restart,reduce_prec=.false.)
           end if
           
         end do
       endif
       if (ilmn) then
          do is = 1, nrhotime
             write(varname, *) "rho-", is
             call decomp_2d_read_one(1,rho1(:,:,:,is),resfile,varname,io_restart,reduce_prec=.false.)
          enddo
          do is = 1, ntime
             write(varname, *) "drho-", is
             call decomp_2d_read_one(1,drho1(:,:,:,is),resfile,varname,io_restart,reduce_prec=.false.)
          enddo
          call decomp_2d_read_one(1,mu1,resfile,"mu",io_restart,reduce_prec=.false.)
       end if

       call decomp_2d_end_io(io_restart, resfile)
       call decomp_2d_close_io(io_restart, resfile)

       ! Read time of restart file
       write(filename,"('data/restart_info/restart',I7.7,'.info')") ifirst-1
       inquire(file=filename, exist=fexists)
       if (nrank==0) write(*,*) filename
       ! Check if file exists
       if (fexists) then
         open(111, file=filename)
         read(111, nml=Time)
         close(111)
         t0 = tfield
         itime0 = 0
       else
         t0 = zero
         itime0 = ifirst-1
       end if
       
    endif

    ! Reconstruction of the dp/dx, dp/dy and dp/dz from pp3
    if (iresflg==0) then
       if (itimescheme <= 4) itr=1
       if (itimescheme == 5) itr=3
       if (itimescheme == 6) itr=5
       call gradp(px1,py1,pz1,pp3)
       if (nrank == 0) write(*,*) 'reconstruction pressure gradients done!'
    end if

    ! Writing restart
    if (iresflg==1) then 
       if (nrank==0) then
          write(fmt1,"(I7.7)") itime
          write(*,*) 'Restart point restart',fmt1,' saved successfully!'!itime/icheckpoint,'saved successfully!'
          ! write(*,*) 'Elapsed time (s)',real(trestart,4)
          ! write(*,*) 'Approximated writing speed (MB/s)',real(((s3df*16.)*1e-6)/trestart,4)
          write(*,*) 'If necessary restart from:',itime+1
       endif
    end if

  end subroutine restart
!-----------------------------------------------------------------------------! 
  subroutine init_restart_adios2()

    use decomp_2d, only : mytype, phG
    use decomp_2d_io, only : decomp_2d_register_variable, decomp_2d_init_io
    use variables, only : numscalar
    use param, only : ilmn, nrhotime, ntime
    use var, only : itimescheme, iibm
    
    implicit none

    integer :: ierror
    
    integer :: is
    character(len=80) :: varname
    
    call decomp_2d_init_io(io_restart)
    
    call decomp_2d_register_variable(io_restart, "ux", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart, "uy", 1, 0, 0, mytype)
    call decomp_2d_register_variable(io_restart, "uz", 1, 0, 0, mytype)

    call decomp_2d_register_variable(io_restart, "pp", 3, 0, 0, mytype, phG) !! XXX: need some way to handle the different grid here...

    do is = 1, numscalar
       write(varname,*) "phi-", is
       call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
    end do

    if ((itimescheme.eq.2) .or. (itimescheme.eq.3)) then
       call decomp_2d_register_variable(io_restart, "dux-2", 1, 0, 0, mytype)
       call decomp_2d_register_variable(io_restart, "duy-2", 1, 0, 0, mytype)
       call decomp_2d_register_variable(io_restart, "duz-2", 1, 0, 0, mytype)

       do is = 1, numscalar
          write(varname,*) "dphi-", is, "-2"
          call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
       end do

       if (itimescheme.eq.3) then
          call decomp_2d_register_variable(io_restart, "dux-3", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "duy-3", 1, 0, 0, mytype)
          call decomp_2d_register_variable(io_restart, "duz-3", 1, 0, 0, mytype)

          do is = 1, numscalar
             write(varname,*) "dphi-", is, "-3"
             call decomp_2d_register_variable(io_restart, trim(varname), 1, 0, 0, mytype)
          end do
       endif
    endif

    if (iibm .ne. 0) then
       call decomp_2d_register_variable(io_restart, "ep", 1, 0, 0, mytype)
    endif

    if (ilmn) then
       do is = 1, nrhotime
          write(varname, *) "rho-", is
          call decomp_2d_register_variable(io_restart, varname, 1, 0, 0, mytype)
       end do
       do is = 1, ntime
          write(varname, *) "drho-", is
          call decomp_2d_register_variable(io_restart, varname, 1, 0, 0, mytype)
       end do
    end if
    
  end subroutine init_restart_adios2
  !-----------------------------------------------------------------------------!
  subroutine apply_spatial_filter(ux1,uy1,uz1,phi1)

    use decomp_2d
    use param
    use var, only: uxf1,uyf1,uzf1,uxf2,uyf2,uzf2,uxf3,uyf3,uzf3,di1,di2,di3,phif1,phif2,phif3
    use variables
    use ibm_param, only : ubcx,ubcy,ubcz

    implicit none
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),            intent(inout)           :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3), numscalar), intent(inout), optional :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi11
    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2, phi2
    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3, phi3

    integer :: i,j,k,npaire

    !if (iscalar == 1) phi11=phi1(:,:,:,1) !currently only first scalar
    if (ifilter==1.or.ifilter==2) then
      call filx(uxf1,ux1,di1,fisx,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0,ubcx)
      call filx(uyf1,uy1,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcy)
      call filx(uzf1,uz1,di1,fisx,fiffxp,fifsxp,fifwxp,xsize(1),xsize(2),xsize(3),1,ubcz)
    else
      uxf1=ux1
      uyf1=uy1
      uzf1=uz1
      !if (iscalar == 1) phif1=phi11
    end if

    call transpose_x_to_y(uxf1,ux2)
    call transpose_x_to_y(uyf1,uy2)
    call transpose_x_to_y(uzf1,uz2)
    !if (iscalar == 1) call transpose_x_to_y(phif1,phi2)

    if (ifilter==1.or.ifilter==3) then ! all filter or y filter
      call fily(uxf2,ux2,di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcx)
      call fily(uyf2,uy2,di2,fisy,fiffy,fifsy,fifwy,ysize(1),ysize(2),ysize(3),0,ubcy)
      call fily(uzf2,uz2,di2,fisy,fiffyp,fifsyp,fifwyp,ysize(1),ysize(2),ysize(3),1,ubcz)
      !if (iscalar.eq.1) call fily(phif2,phi2,di2,fisy,fiffy,fifsy,fifwy,ysize(1),ysize(2),ysize(3),0)
    else
      uxf2=ux2
      uyf2=uy2
      uzf2=uz2
      !if (iscalar == 1) phif2=phi2
    end if

    call transpose_y_to_z(uxf2,ux3)
    call transpose_y_to_z(uyf2,uy3)
    call transpose_y_to_z(uzf2,uz3)
    !if (iscalar == 1) call transpose_y_to_z(phif2,phi3)

    if (ifilter==1.or.ifilter==2) then
      call filz(uxf3,ux3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcx)
      call filz(uyf3,uy3,di3,fisz,fiffzp,fifszp,fifwzp,zsize(1),zsize(2),zsize(3),1,ubcy)
      call filz(uzf3,uz3,di3,fisz,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),0,ubcz)
      !if (iscalar.eq.1) call filz(phif3,phi3,di3,fisz,fiffz,fifsz,fifwz,zsize(1),zsize(2),zsize(3),0)
    else
      uxf3=ux3
      uyf3=uy3
      uzf3=uz3
      !if (iscalar == 1) phif3=phi3
    end if

    call transpose_z_to_y(uxf3,ux2)
    call transpose_z_to_y(uyf3,uy2)
    call transpose_z_to_y(uzf3,uz2)
    !if (iscalar == 1) call transpose_z_to_y(phif3,phi2)

    call transpose_y_to_x(ux2,ux1)
    call transpose_y_to_x(uy2,uy1)
    call transpose_y_to_x(uz2,uz1)
    !if (iscalar == 1) call transpose_y_to_x(phi2,phi11)

    !if (iscalar == 1) phi1(:,:,:,1)=phi11

  end subroutine apply_spatial_filter
  
  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: compute_cfldiff
  ! DESCRIPTION: Computes Diffusion/Numerical Fourier number (D < 0.5).
  !      AUTHOR: Kay Schäfer
  !-----------------------------------------------------------------------------!
  subroutine compute_cfldiff()
     
     use param, only : xnu,dt,dx,dy,dz,istret
     use param, only : cfl_diff_sum, cfl_diff_x, cfl_diff_y, cfl_diff_z
     use variables, only : dyp
     use decomp_2d, only : nrank

     implicit none

     cfl_diff_x = xnu * dt/ (dx**2)
     cfl_diff_z = xnu * dt/ (dz**2)

     if (istret == 0) then
        cfl_diff_y = xnu * dt / (dy**2)
     else
        cfl_diff_y = xnu * dt / ((minval(dyp))**2)
     end if

     cfl_diff_sum = cfl_diff_x + cfl_diff_y + cfl_diff_z

     if (nrank==0) then
        write(*,*) 'Diffusion number D (or numerical Fourier, Fo)'
        write(*,"(' D,x                    : ',F17.8)") cfl_diff_x
        write(*,"(' D,y                    : ',F17.8)") cfl_diff_y
        write(*,"(' D,z                    : ',F17.8)") cfl_diff_z
        write(*,"(' D,sum                  : ',F17.8)") cfl_diff_sum
        write(*,*) '-----------------------------------------------------------'
     endif

     return
  end subroutine compute_cfldiff

  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: compute_cfl
  ! DESCRIPTION: Computes CFl number for stretched mesh
  !              and adjust time-step if required.
  !      AUTHOR: Kay Schäfer, Filippo Moroni
  !-----------------------------------------------------------------------------!
  subroutine compute_cfl(ux,uy,uz)
    
    use param,       only : dx,dy,dz,dt,istret,cfl_limit,icfllim,t,itimescheme,cflmax
    use decomp_2d,   only : nrank, mytype, xsize, xstart, xend, real_type
    use mpi
    use variables,   only : dyp

    implicit none

    integer      :: code, i,j,k,jloc
    real(mytype) :: value_x, value_y, value_z, value_sum
    real(mytype) :: maxvalue_sum, maxvalue_sum_out, maxvalue_x, maxvalue_y,  maxvalue_z
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(4) :: cflmax_in, cflmax_out
    !
    maxvalue_x  =-1609._mytype
    maxvalue_y  =-1609._mytype
    maxvalue_z  =-1609._mytype
    maxvalue_sum=-1609._mytype
    !
    if (istret == 0) then
       do j = xstart(2), xend(2)
          jloc = j-xstart(2)+1
          value_x    = maxval(abs(ux(:,jloc,:)) / dx)
          value_y    = maxval(abs(uy(:,jloc,:)) / dy)
          value_z    = maxval(abs(uz(:,jloc,:)) / dz)
          value_sum  = maxval(abs(ux(:,jloc,:)) / dx + abs(uy(:,jloc,:)) / dy + abs(uz(:,jloc,:)) / dz)
          !
          maxvalue_x   = maxval((/maxvalue_x,   value_x /))
          maxvalue_y   = maxval((/maxvalue_y,   value_y /))
          maxvalue_z   = maxval((/maxvalue_z,   value_z /))
          maxvalue_sum = maxval((/maxvalue_sum, value_sum /))
       end do
    else
       do j = xstart(2), xend(2)
          jloc = j-xstart(2)+1
          value_x    = maxval(abs(ux(:,jloc,:)) / dx)
          value_y    = maxval(abs(uy(:,jloc,:)) / dyp(j))
          value_z    = maxval(abs(uz(:,jloc,:)) / dz)
          value_sum  = maxval(abs(ux(:,jloc,:)) / dx + abs(uy(:,jloc,:)) / dyp(j) + abs(uz(:,jloc,:)) /dz)
          !
          maxvalue_x   = maxval((/maxvalue_x,   value_x /))
          maxvalue_y   = maxval((/maxvalue_y,   value_y /))
          maxvalue_z   = maxval((/maxvalue_z,   value_z /))
          maxvalue_sum = maxval((/maxvalue_sum, value_sum /))
       end do
    end if

    cflmax_in =  (/maxvalue_x, maxvalue_y, maxvalue_z, maxvalue_sum/)

    call MPI_ALLREDUCE(cflmax_in,cflmax_out,4,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    
    if (nrank == 0) then
      write(*,*) '-----------------------------------------------------------'
      write(*,*) 'CFL Number (or Courant, Co)'
      write(*,"(' CFL,x                  : ',F17.8)") cflmax_out(1) * dt
      write(*,"(' CFL,y                  : ',F17.8)") cflmax_out(2) * dt
      write(*,"(' CFL,z                  : ',F17.8)") cflmax_out(3) * dt
      write(*,"(' CFL,sum                : ',F17.8)") cflmax_out(4) * dt
      write(*,*) '-----------------------------------------------------------'
    end if
    
    ! Store the maximum CFL
    cflmax = cflmax_out(4) * dt
    
    ! Adjust time-step if adjustable time-step option is enabled and if we are exiting the specified interval (cfl_lim - 0.05, cfl_lim)
    ! Valid only for RK3
    if (icfllim == 1 .and. itimescheme .eq. 5 .and. (cflmax_out(4)*dt > cfl_limit .or. cflmax_out(4)*dt < (cfl_limit - 0.05))) then
    
        dt = (cfl_limit / cflmax_out(4))
        
        ! Update coefficients for time integration schemes
        call update_time_int_coeff()
                   
    end if
  end subroutine compute_cfl

  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: compute_reynolds_cell
  ! DESCRIPTION: Computes numerical Péclet number or Reynolds cell number 
  !              for a stretched mesh.
  !      AUTHOR: Filippo Moroni (adapted from compute_cfl subroutine)
  !-----------------------------------------------------------------------------!
  subroutine compute_reynolds_cell(ux,uy,uz)
    
    use param, only : dx,dy,dz,dt,istret,xnu
    use decomp_2d, only : nrank, mytype, xsize, xstart, xend, real_type
    use mpi
    use variables, only : dyp

    implicit none

    integer      :: code,i,j,k,jloc
    real(mytype) :: value_x, value_y, value_z, value_sum
    real(mytype) :: maxvalue_sum, maxvalue_sum_out, maxvalue_x, maxvalue_y,  maxvalue_z
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    
    ! Reynolds cell max in/out
    real(mytype),dimension(4) :: recmax_in, recmax_out
   
    maxvalue_x  =-1609._mytype
    maxvalue_y  =-1609._mytype
    maxvalue_z  =-1609._mytype
    maxvalue_sum=-1609._mytype
    
    ! Uniform mesh along y
    if (istret == 0) then
       do j = xstart(2), xend(2)
          jloc = j-xstart(2)+1
          value_x    = maxval(abs(ux(:,jloc,:)) * dx)
          value_y    = maxval(abs(uy(:,jloc,:)) * dy)
          value_z    = maxval(abs(uz(:,jloc,:)) * dz)
          value_sum  = maxval(abs(ux(:,jloc,:)) * dx + abs(uy(:,jloc,:)) * dy + abs(uz(:,jloc,:)) * dz)
          
          maxvalue_x   = maxval((/maxvalue_x,   value_x /))
          maxvalue_y   = maxval((/maxvalue_y,   value_y /))
          maxvalue_z   = maxval((/maxvalue_z,   value_z /))
          maxvalue_sum = maxval((/maxvalue_sum, value_sum /))
       end do
    ! Stretched mesh along y
    else
       do j = xstart(2), xend(2)
          jloc = j-xstart(2)+1
          value_x    = maxval(abs(ux(:,jloc,:)) * dx)
          value_y    = maxval(abs(uy(:,jloc,:)) * dyp(j))
          value_z    = maxval(abs(uz(:,jloc,:)) * dz)
          value_sum  = maxval(abs(ux(:,jloc,:)) * dx + abs(uy(:,jloc,:)) * dyp(j) + abs(uz(:,jloc,:)) * dz)
          
          maxvalue_x   = maxval((/maxvalue_x,   value_x /))
          maxvalue_y   = maxval((/maxvalue_y,   value_y /))
          maxvalue_z   = maxval((/maxvalue_z,   value_z /))
          maxvalue_sum = maxval((/maxvalue_sum, value_sum /))
       end do
    end if

    recmax_in =  (/maxvalue_x, maxvalue_y, maxvalue_z, maxvalue_sum/)

    call MPI_REDUCE(recmax_in,recmax_out,4,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)

    if (nrank == 0) then
      write(*,*) 'Reynolds cell (or numerical Péclet, Pé)'
      write(*,"(' Pé,x                   : ',F17.8)") recmax_out(1) / xnu
      write(*,"(' Pé,y                   : ',F17.8)") recmax_out(2) / xnu
      write(*,"(' Pé,z                   : ',F17.8)") recmax_out(3) / xnu
      write(*,"(' Pé,sum                 : ',F17.8)") recmax_out(4) / xnu
      write(*,*) '-----------------------------------------------------------'
    end if
  end subroutine compute_reynolds_cell

  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: compute_stab_param
  ! DESCRIPTION: Computes stability parameter S < 1 (Thompson et al. (1985)). 
  !      AUTHOR: Filippo Moroni (adapted from compute_cfl subroutine)
  !-----------------------------------------------------------------------------!
  subroutine compute_stab_param(ux,uy,uz)
    
    use param, only : dt,xnu,two
    use decomp_2d, only : nrank, mytype, xsize, xstart, xend, real_type
    use mpi

    implicit none

    integer      :: code,i,j,k,jloc
    real(mytype) :: value_x, value_y, value_z, value_sum
    real(mytype) :: maxvalue_sum, maxvalue_sum_out, maxvalue_x, maxvalue_y,  maxvalue_z
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    
    ! Stability parameter S max in/out
    real(mytype),dimension(4) :: stparmax_in, stparmax_out
   
    maxvalue_x  =-1609._mytype
    maxvalue_y  =-1609._mytype
    maxvalue_z  =-1609._mytype
    maxvalue_sum=-1609._mytype
    
    ! Unique cycle, since S does not depend on the mesh but just on velocity field
    do j = xstart(2), xend(2)
       jloc = j-xstart(2)+1
       value_x    = maxval((ux(:,jloc,:))**2)
       value_y    = maxval((uy(:,jloc,:))**2)
       value_z    = maxval((uz(:,jloc,:))**2)
       value_sum  = maxval((ux(:,jloc,:))**2 + (uy(:,jloc,:))**2 + (uz(:,jloc,:))**2)
          
       maxvalue_x   = maxval((/maxvalue_x,   value_x /))
       maxvalue_y   = maxval((/maxvalue_y,   value_y /))
       maxvalue_z   = maxval((/maxvalue_z,   value_z /))
       maxvalue_sum = maxval((/maxvalue_sum, value_sum /))
    end do

    stparmax_in =  (/maxvalue_x, maxvalue_y, maxvalue_z, maxvalue_sum/)

    call    MPI_REDUCE(stparmax_in,stparmax_out,4,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)

    if (nrank == 0) then
      write(*,*) 'Stability parameter S (Thompson et al. (1985))'
      write(*,"(' S,x                    : ',F17.8)") stparmax_out(1) * dt / two / xnu
      write(*,"(' S,y                    : ',F17.8)") stparmax_out(2) * dt / two / xnu
      write(*,"(' S,z                    : ',F17.8)") stparmax_out(3) * dt / two / xnu
      write(*,"(' S,sum                  : ',F17.8)") stparmax_out(4) * dt / two / xnu
      write(*,*) '-----------------------------------------------------------'
    end if
  end subroutine compute_stab_param
  
  !---------------------------------------------------------------------------!
  ! Update coefficients for time integration schemes
  ! if adaptive time-step is used (max CFL condition)
  ! (taken from init_variables in 'variables' module).
  !
  ! - Used in subroutine for CFL calculation in tools module.
  ! - At the moment, valid only for RK3 scheme.
  ! - To do: AB schemes.
  !---------------------------------------------------------------------------!
  subroutine update_time_int_coeff()
  
  use param
  
  implicit none
  
    ! Runge-Kutta 3 (RK3)
    adt(1)=(eight/fifteen)*dt
    bdt(1)= zero
    gdt(1)=adt(1)
    adt(2)=(      five/twelve)*dt
    bdt(2)=(-seventeen/ sixty)*dt
    gdt(2)=adt(2)+bdt(2)
    adt(3)=( three/four)*dt
    bdt(3)=(-five/twelve)*dt
    gdt(3)=adt(3)+bdt(3)
           
  end subroutine update_time_int_coeff
  
  !-----------------------------------------------------------------------------!  
  ! Rescale pressure to physical pressure
  ! Written by Kay Schäfer 2019
  !-----------------------------------------------------------------------------!
  elemental subroutine rescale_pressure(pre1)

    use decomp_2d, only : mytype
    use param, only : itimescheme, gdt
    
    implicit none

    real(mytype), intent(inout) :: pre1

    ! Adjust pressure to physical pressure
    ! Multiply pressure by factor of time-scheme
    ! 1/gdt = 1  / (dt * c_k)
    
    ! Explicit Euler, AB2, AB3, AB4, RK3
    if (itimescheme>=1 .and. itimescheme<=5) then
       pre1 = pre1 / gdt(3)
    ! RK4
    elseif (itimescheme==6) then
       pre1 = pre1 / gdt(5)
    endif

  end subroutine
  
  !-----------------------------------------------------------------------------!
  subroutine mean_plane_x (f1,nx,ny,nz,fm1)

    use param, only : mytype, zero

    implicit none

    integer,intent(in) :: nx, ny, nz
    real(mytype),intent(in),dimension(nx,ny,nz) :: f1
    real(mytype),intent(out),dimension(ny,nz) :: fm1
    integer :: i,j,k

    fm1 = sum(f1, DIM=1) / real(nx, mytype)
    return

  end subroutine mean_plane_x
  
  !-----------------------------------------------------------------------------!
  subroutine mean_plane_y (f2,nx,ny,nz,fm2)

    use param, only : mytype, zero

    implicit none

    integer,intent(in) :: nx, ny, nz
    real(mytype),intent(in),dimension(nx,ny,nz) :: f2
    real(mytype),intent(out),dimension(nx,nz) :: fm2
    integer :: i,j,k

    fm2 = sum(f2, DIM=2) / real(ny, mytype)
    return

  end subroutine mean_plane_y
  
  !-----------------------------------------------------------------------------!
  subroutine mean_plane_z (f3,nx,ny,nz,fm3)

    use param, only : mytype, zero

    implicit none

    integer,intent(in) :: nx, ny, nz
    real(mytype),intent(in),dimension(nx,ny,nz) :: f3
    real(mytype),intent(out),dimension(nx,ny) :: fm3
    integer :: i,j,k

    fm3 = sum(f3, DIM=3) / real(nz,mytype)
    return

  end subroutine mean_plane_z

  !-----------------------------------------------------------------------------!
  !  SUBROUTINE: avg3d
  !      AUTHOR: Stefano Rolfo
  ! DESCRIPTION: Compute the total sum of a a 3d field
  !-----------------------------------------------------------------------------!
  subroutine avg3d (var, avg)

    use decomp_2d, only: real_type, xsize, xend
    use param
    use variables, only: nx,ny,nz,nxm,nym,nzm
    use mpi

    implicit none

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: var
    real(mytype), intent(out) :: avg
    real(mytype)              :: dep

    integer :: i,j,k, code
    integer :: nxc, nyc, nzc, xsize1, xsize2, xsize3

    if (nclx1==1.and.xend(1)==nx) then
       xsize1=xsize(1)-1
    else
       xsize1=xsize(1)
    endif
    if (ncly1==1.and.xend(2)==ny) then
       xsize2=xsize(2)-1
    else
       xsize2=xsize(2)
    endif
    if (nclz1==1.and.xend(3)==nz) then
       xsize3=xsize(3)-1
    else
       xsize3=xsize(3)
    endif
    if (nclx1==1) then
       nxc=nxm
    else
       nxc=nx
    endif
    if (ncly1==1) then
       nyc=nym
    else
       nyc=ny
    endif
    if (nclz1==1) then
       nzc=nzm
    else
       nzc=nz
    endif

    dep=zero
    do k=1,xsize3
       do j=1,xsize2
          do i=1,xsize1
             !dep=dep+var(i,j,k)**2
             dep=dep+var(i,j,k)
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(dep,avg,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    avg=avg/(nxc*nyc*nzc)

    return

  end subroutine avg3d
end module tools

!-----------------------------------------------------------------------------!
! Subroutine for computing the local and global CFL
! number, according to Lele 1992.
!-----------------------------------------------------------------------------!
subroutine cfl_compute(uxmax,uymax,uzmax)

  use param
  use variables
  use var

  implicit none

  real(mytype),intent(in) :: uxmax,uymax,uzmax
  real(mytype) :: cfl_x_adv,cfl_x_diff,cfl_y_adv,cfl_y_diff,cfl_z_adv,cfl_z_diff
  real(mytype) :: cfl_conv_lim, cfl_diff_lim
  real(mytype) :: sigma_conv(3), sigma_diff(3)
  real(mytype) :: visc

  ! Set the constants (this is true for periodic boundaries)
  sigma_conv=[zero, sqrt(three), 2.85_mytype]
  sigma_diff=[two, 2.5_mytype, 2.9_mytype]

  if(jles==0) then
     visc=xnu
  elseif (jles==1) then
     visc=xnu
  endif

  ! This is considering 1D peridic boundaries
  ! Do x-direction
  cfl_x_adv =abs(uxmax) * dt / dx
  cfl_x_diff = visc * dt / dx**2
  ! Do y-direction
  cfl_y_adv = abs(uymax) * dt / dy
  cfl_y_diff = visc * dt / dy**2
  ! Do z-direction
  cfl_z_adv = abs(uzmax) * dt / dz
  cfl_z_diff = visc * dt / dz**2

  ! So far we will focus on uniform grids
  if(nrank == 0) then
     write(*,*) ' '
     write(*,1002) cfl_x_adv, cfl_x_diff
1002 format('CFL x-direction (Adv and Diff) =',F9.4,',',F9.4)
     write(*,1003) cfl_y_adv, cfl_y_diff
1003 format('CFL y-direction (Adv and Diff) =',F9.4,',',F9.4)
     write(*,1004) cfl_z_adv, cfl_z_diff
1004 format('CFL z-direction (Adv and Diff) =',F9.4,',',F9.4)
     cfl_conv_lim = sigma_conv(itimescheme) / sqrt(three)
     cfl_diff_lim = sigma_diff(itimescheme) / six
     write(*,1005) cfl_conv_lim, cfl_diff_lim
     write(*,*) ' '
1005 format('CFL limits (Adv and Diff) : ',F9.4,',',F9.4)
  endif

end subroutine cfl_compute

!-----------------------------------------------------------------------------!
subroutine stretching()

  use decomp_2d
  use variables
  use param
  use var
  use mpi

  implicit none

  real(mytype) :: yinf,den,xnum,xcx,den1,den2,den3,den4,xnum1,cst
  integer :: j

  yinf=-yly/two
  den=two*beta*yinf
  xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
  alpha=abs(xnum/den)
  xcx=one/beta/alpha
  if (alpha.ne.0.) then
     if (istret.eq.1) yp(1)=zero
     if (istret.eq.2) yp(1)=zero
     if (istret.eq.1) yeta(1)=zero
     if (istret.eq.2) yeta(1)=-half
     if (istret.eq.3) yp(1)=zero
     if (istret.eq.3) yeta(1)=-half
     do j=2,ny
        if (istret==1) yeta(j)=real(j-1,mytype)*(one/nym)
        if (istret==2) yeta(j)=real(j-1,mytype)*(one/nym)-half
        if (istret==3) yeta(j)=real(j-1,mytype)*(half/nym)-half
        den1=sqrt(alpha*beta+one)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
        den=two*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
        den4=two*alpha*beta-cos(two*pi*yeta(j))+one
        xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(two*sqrt(alpha)*sqrt(alpha*beta+one))
        if (istret==1) then
           if (yeta(j).lt.half) yp(j)=xnum1-cst-yinf
           if (yeta(j).eq.half) yp(j)=zero-yinf
           if (yeta(j).gt.half) yp(j)=xnum1+cst-yinf
        endif
        if (istret==2) then
           if (yeta(j).lt.half) yp(j)=xnum1-cst+yly
           if (yeta(j).eq.half) yp(j)=zero+yly
           if (yeta(j).gt.half) yp(j)=xnum1+cst+yly
        endif
        if (istret==3) then
           if (yeta(j).lt.half) yp(j)=(xnum1-cst+yly)*two
           if (yeta(j).eq.half) yp(j)=(zero+yly)*two
           if (yeta(j).gt.half) yp(j)=(xnum1+cst+yly)*two
        endif
     enddo
  endif
  if (alpha.eq.0.) then
     yp(1)=-1.e10
     do j=2,ny
        yeta(j)=real(j-1,mytype)*(one/ny)
        yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
     enddo
  endif
  if (alpha.ne.0.) then
     do j=1,ny
        if (istret==1) yetai(j)=(real(j,mytype)-half)*(one/nym)
        if (istret==2) yetai(j)=(real(j,mytype)-half)*(one/nym)-half
        if (istret==3) yetai(j)=(real(j,mytype)-half)*(half/nym)-half
        den1=sqrt(alpha*beta+one)
        xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
        den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
        den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
        den4=two*alpha*beta-cos(two*pi*yetai(j))+one
        xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
        cst=sqrt(beta)*pi/(two*sqrt(alpha)*sqrt(alpha*beta+one))
        if (istret==1) then
           if (yetai(j).lt.half) ypi(j)=xnum1-cst-yinf
           if (yetai(j).eq.half) ypi(j)=zero-yinf
           if (yetai(j).gt.half) ypi(j)=xnum1+cst-yinf
        endif
        if (istret==2) then
           if (yetai(j).lt.half) ypi(j)=xnum1-cst+yly
           if (yetai(j).eq.half) ypi(j)=zero+yly
           if (yetai(j).gt.half) ypi(j)=xnum1+cst+yly
        endif
        if (istret==3) then
           if (yetai(j).lt.half) ypi(j)=(xnum1-cst+yly)*two
           if (yetai(j).eq.half) ypi(j)=(zero+yly)*two
           if (yetai(j).gt.half) ypi(j)=(xnum1+cst+yly)*two
        endif
     enddo
  endif
  if (alpha.eq.0.) then
     ypi(1)=-1.e10
     do j=2,ny
        yetai(j)=real(j-1,mytype)*(one/ny)
        ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
     enddo
  endif

  ! Mapping, metric terms
  if (istret .ne. 3) then
     do j=1,ny
        ppy(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yeta(j))*sin(pi*yeta(j)))
        pp2y(j)=ppy(j)*ppy(j)
        pp4y(j)=(-two/beta*cos(pi*yeta(j))*sin(pi*yeta(j)))
     enddo
     do j=1,ny
        ppyi(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yetai(j))*sin(pi*yetai(j)))
        pp2yi(j)=ppyi(j)*ppyi(j)
        pp4yi(j)=(-two/beta*cos(pi*yetai(j))*sin(pi*yetai(j)))
     enddo
  endif

  if (istret .eq. 3) then
     do j=1,ny
        ppy(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yeta(j))*sin(pi*yeta(j)))
        pp2y(j)=ppy(j)*ppy(j)
        pp4y(j)=(-two/beta*cos(pi*yeta(j))*sin(pi*yeta(j)))/two
     enddo
     do j=1,ny
        ppyi(j)=yly*(alpha/pi+(one/pi/beta)*sin(pi*yetai(j))*sin(pi*yetai(j)))
        pp2yi(j)=ppyi(j)*ppyi(j)
        pp4yi(j)=(-two/beta*cos(pi*yetai(j))*sin(pi*yetai(j)))/two
     enddo
  endif

  if (nrank == 0) then
     open(10,file='yp.dat', form='formatted')
     do j=1,ny
        write(10,*)yp(j)
     enddo
     close(10)
     open(10,file='ypi.dat', form='formatted')
     do j=1,nym
        write(10,*)ypi(j)
     enddo
     close(10)
  endif

end subroutine stretching

!-----------------------------------------------------------------------------!
subroutine inversion5_v1(aaa_in,eee,spI)

  use decomp_2d
  use variables
  use param
  use var
  use mpi

  implicit none

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

  complex(mytype),dimension(spI%yst(1):spI%yen(1),ny/2,spI%yst(3):spI%yen(3),5) :: aaa, aaa_in
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(2):spI%yen(2),spI%yst(3):spI%yen(3)) :: eee
  integer :: i,j,k,m,mi,jc
  integer,dimension(2) :: ja,jb
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

  real(mytype) :: tmp1,tmp2,tmp3,tmp4

  complex(mytype) :: cx
  real(mytype) :: rl, iy
  external cx, rl, iy

  aaa = aaa_in

  do i = 1, 2
     ja(i) = 4 - i
     jb(i) = 5 - i
  enddo
  do m = 1, ny/2 - 2
     do i = 1, 2
        mi = m + i
        do k = spI%yst(3), spI%yen(3)
           do j = spI%yst(1), spI%yen(1)
              if (rl(aaa(j,m,k,3)) /= zero) tmp1 = rl(aaa(j,mi,k,3-i)) / rl(aaa(j,m,k,3))
              if (iy(aaa(j,m,k,3)) /= zero) tmp2 = iy(aaa(j,mi,k,3-i)) / iy(aaa(j,m,k,3))
              sr(j,k)=cx(tmp1,tmp2)
              eee(j,mi,k)=cx(rl(eee(j,mi,k)) - tmp1 * rl(eee(j,m,k)),&
                             iy(eee(j,mi,k)) - tmp2 * iy(eee(j,m,k)))
           enddo
        enddo
        do jc = ja(i), jb(i)
           do k = spI%yst(3), spI%yen(3)
              do j = spI%yst(1), spI%yen(1)
                 aaa(j,mi,k,jc) = cx(rl(aaa(j,mi,k,jc)) - rl(sr(j,k)) * rl(aaa(j,m,k,jc+i)),&
                                     iy(aaa(j,mi,k,jc)) - iy(sr(j,k)) * iy(aaa(j,m,k,jc+i)))
              enddo
           enddo
        enddo
     enddo
  enddo

  do k = spI%yst(3), spI%yen(3)
     do j = spI%yst(1), spI%yen(1)
        if (abs(rl(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp1 = rl(aaa(j,ny/2,k,2)) / rl(aaa(j,ny/2-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs(iy(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp2 = iy(aaa(j,ny/2,k,2)) / iy(aaa(j,ny/2-1,k,3))
        else
           tmp2 = zero
        endif
        sr(j,k) = cx(tmp1,tmp2)
        b1(j,k) = cx(rl(aaa(j,ny/2,k,3)) - tmp1 * rl(aaa(j,ny/2-1,k,4)),&
                     iy(aaa(j,ny/2,k,3)) - tmp2 * iy(aaa(j,ny/2-1,k,4)))

        if (abs(rl(b1(j,k))) > epsilon) then
           tmp1 = rl(sr(j,k)) / rl(b1(j,k))
           tmp3 = rl(eee(j,ny/2,k)) / rl(b1(j,k)) - tmp1 * rl(eee(j,ny/2-1,k))
        else
           tmp1 = zero
           tmp3 = zero
        endif
        if (abs(iy(b1(j,k))) > epsilon) then
           tmp2 = iy(sr(j,k)) / iy(b1(j,k))
           tmp4 = iy(eee(j,ny/2,k)) / iy(b1(j,k)) - tmp2 * iy(eee(j,ny/2-1,k))
        else
           tmp2 = zero
           tmp4 = zero
        endif
        a1(j,k) = cx(tmp1,tmp2)
        eee(j,ny/2,k) = cx(tmp3,tmp4)

        if (abs(rl(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp1 = one / rl(aaa(j,ny/2-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs(iy(aaa(j,ny/2-1,k,3))) > epsilon) then
           tmp2 = one / iy(aaa(j,ny/2-1,k,3))
        else
           tmp2 = zero
        endif
        b1(j,k) = cx(tmp1, tmp2)
        a1(j,k) = cx(rl(aaa(j,ny/2-1,k,4)) * rl(b1(j,k)),&
                     iy(aaa(j,ny/2-1,k,4)) * iy(b1(j,k)))
        eee(j,ny/2-1,k) = cx(rl(eee(j,ny/2-1,k)) * rl(b1(j,k)) - rl(a1(j,k)) * rl(eee(j,ny/2,k)),&
                             iy(eee(j,ny/2-1,k)) * iy(b1(j,k)) - iy(a1(j,k)) * iy(eee(j,ny/2,k)))
     enddo
  enddo

  do i = ny/2 - 2, 1, -1
     do k = spI%yst(3), spI%yen(3)
        do j = spI%yst(1), spI%yen(1)
           if (abs(rl(aaa(j,i,k,3))) > epsilon) then
              tmp1 = one / rl(aaa(j,i,k,3))
           else
              tmp1 = zero
           endif
           if (abs(iy(aaa(j,i,k,3))) > epsilon) then
              tmp2 = one/iy(aaa(j,i,k,3))
           else
              tmp2 = zero
           endif
           sr(j,k) = cx(tmp1,tmp2)
           a1(j,k) = cx(rl(aaa(j,i,k,4)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,4)) * iy(sr(j,k)))
           b1(j,k) = cx(rl(aaa(j,i,k,5)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,5)) * iy(sr(j,k)))
           eee(j,i,k) = cx(rl(eee(j,i,k)) * rl(sr(j,k)) - rl(a1(j,k)) * rl(eee(j,i+1,k)) - rl(b1(j,k)) * rl(eee(j,i+2,k)),&
                           iy(eee(j,i,k)) * iy(sr(j,k)) - iy(a1(j,k)) * iy(eee(j,i+1,k)) - iy(b1(j,k)) * iy(eee(j,i+2,k)))
        enddo
     enddo
  enddo

  return

end subroutine inversion5_v1

!-----------------------------------------------------------------------------!
subroutine inversion5_v2(aaa,eee,spI)

  use decomp_2d
  use variables
  use param
  use var
  use MPI

  implicit none

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

  complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3),5) :: aaa
  complex(mytype),dimension(spI%yst(1):spI%yen(1),nym,spI%yst(3):spI%yen(3)) :: eee
  integer :: i,j,k,m,mi,jc
  integer,dimension(2) :: ja,jb
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
  complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

  real(mytype) :: tmp1,tmp2,tmp3,tmp4

  complex(mytype) :: cx
  real(mytype) :: rl, iy
  external cx, rl, iy

  do i = 1, 2
     ja(i) = 4 - i
     jb(i) = 5 - i
  enddo
  do m = 1, nym - 2
     do i = 1, 2
        mi = m + i
        do k = spI%yst(3), spI%yen(3)
           do j = spI%yst(1), spI%yen(1)
              if (rl(aaa(j,m,k,3)) /= zero) tmp1 = rl(aaa(j,mi,k,3-i)) / rl(aaa(j,m,k,3))
              if (iy(aaa(j,m,k,3)) /= zero) tmp2 = iy(aaa(j,mi,k,3-i)) / iy(aaa(j,m,k,3))
              sr(j,k) = cx(tmp1, tmp2)
              eee(j,mi,k) = cx(rl(eee(j,mi,k)) - tmp1 * rl(eee(j,m,k)),&
                               iy(eee(j,mi,k)) - tmp2 * iy(eee(j,m,k)))
           enddo
        enddo
        do jc = ja(i), jb(i)
           do k = spI%yst(3), spI%yen(3)
              do j = spI%yst(1), spI%yen(1)
                 aaa(j,mi,k,jc) = cx(rl(aaa(j,mi,k,jc)) - rl(sr(j,k)) * rl(aaa(j,m,k,jc+i)),&
                                     iy(aaa(j,mi,k,jc)) - iy(sr(j,k)) * iy(aaa(j,m,k,jc+i)))
              enddo
           enddo
        enddo
     enddo
  enddo
  do k = spI%yst(3), spI%yen(3)
     do j = spI%yst(1), spI%yen(1)
        if (abs(rl(aaa(j,nym-1,k,3))) > epsilon) then
           tmp1 = rl(aaa(j,nym,k,2)) / rl(aaa(j,nym-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs(iy(aaa(j,nym-1,k,3))) > epsilon) then
           tmp2 = iy(aaa(j,nym,k,2)) / iy(aaa(j,nym-1,k,3))
        else
           tmp2 = zero
        endif
        sr(j,k) = cx(tmp1,tmp2)
        b1(j,k) = cx(rl(aaa(j,nym,k,3)) - tmp1 * rl(aaa(j,nym-1,k,4)),&
                     iy(aaa(j,nym,k,3)) - tmp2 * iy(aaa(j,nym-1,k,4)))
        if (abs(rl(b1(j,k))) > epsilon) then
           tmp1 = rl(sr(j,k)) / rl(b1(j,k))
           tmp3 = rl(eee(j,nym,k)) / rl(b1(j,k)) - tmp1 * rl(eee(j,nym-1,k))
        else
           tmp1 = zero
           tmp3 = zero
        endif
        if (abs(iy(b1(j,k))) > epsilon) then
           tmp2 = iy(sr(j,k)) / iy(b1(j,k))
           tmp4 = iy(eee(j,nym,k)) / iy(b1(j,k)) - tmp2 * iy(eee(j,nym-1,k))
        else
           tmp2 = zero
           tmp4 = zero
        endif
        a1(j,k) = cx(tmp1, tmp2)
        eee(j,nym,k) = cx(tmp3, tmp4)

        if (abs(rl(aaa(j,nym-1,k,3))) > epsilon) then
           tmp1 = one / rl(aaa(j,nym-1,k,3))
        else
           tmp1 = zero
        endif
        if (abs(iy(aaa(j,nym-1,k,3))) > epsilon) then
           tmp2 = one / iy(aaa(j,nym-1,k,3))
        else
           tmp2 = zero
        endif
        b1(j,k) = cx(tmp1,tmp2)
        a1(j,k) = cx(rl(aaa(j,nym-1,k,4)) * rl(b1(j,k)),&
                     iy(aaa(j,nym-1,k,4)) * iy(b1(j,k)))
        eee(j,nym-1,k) = cx(rl(eee(j,nym-1,k)) * rl(b1(j,k)) - rl(a1(j,k)) * rl(eee(j,nym,k)),&
                            iy(eee(j,nym-1,k)) * iy(b1(j,k)) - iy(a1(j,k)) * iy(eee(j,nym,k)))
     enddo
  enddo

  do i = nym - 2, 1, -1
     do k = spI%yst(3), spI%yen(3)
        do j = spI%yst(1), spI%yen(1)
           if (abs(rl(aaa(j,i,k,3))) > epsilon) then
              tmp1 = one / rl(aaa(j,i,k,3))
           else
              tmp1 = zero
           endif
           if (abs(iy(aaa(j,i,k,3))) > epsilon) then
              tmp2 = one / iy(aaa(j,i,k,3))
           else
              tmp2 = zero
           endif
           sr(j,k) = cx(tmp1,tmp2)
           a1(j,k) = cx(rl(aaa(j,i,k,4)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,4)) * iy(sr(j,k)))
           b1(j,k) = cx(rl(aaa(j,i,k,5)) * rl(sr(j,k)),&
                        iy(aaa(j,i,k,5)) * iy(sr(j,k)))
           eee(j,i,k) = cx(rl(eee(j,i,k)) * rl(sr(j,k)) - rl(a1(j,k)) * rl(eee(j,i+1,k)) -rl(b1(j,k)) * rl(eee(j,i+2,k)),&
                           iy(eee(j,i,k)) * iy(sr(j,k)) - iy(a1(j,k)) * iy(eee(j,i+1,k)) -iy(b1(j,k)) * iy(eee(j,i+2,k)))
        enddo
     enddo
  enddo

  return

end subroutine inversion5_v2
!-----------------------------------------------------------------------------!
function rl(complexnumber)

  use param

  implicit none

  real(mytype) :: rl
  complex(mytype) :: complexnumber

  rl = real(complexnumber, kind=mytype)

end function rl
!-----------------------------------------------------------------------------!
function iy(complexnumber)

  use param

  implicit none

  real(mytype) :: iy
  complex(mytype) :: complexnumber

  iy = aimag(complexnumber)

end function iy
!-----------------------------------------------------------------------------!
function cx(realpart,imaginarypart)

  use param

  implicit none

  complex(mytype) :: cx
  real(mytype) :: realpart, imaginarypart

  cx = cmplx(realpart, imaginarypart, kind=mytype)

end function cx
!-----------------------------------------------------------------------------!
subroutine calc_temp_eos(temp, rho, phi, mweight, xlen, ylen, zlen)

  use decomp_2d
  use param, only : pressure0, imultispecies
  use var, only : numscalar

  implicit none

  !! inputs
  integer, intent(in) :: xlen, ylen, zlen
  real(mytype), intent(in), dimension(xlen, ylen, zlen) :: rho
  real(mytype), intent(in), dimension(xlen, ylen, zlen, numscalar) :: phi

  !! outputs
  real(mytype), intent(out), dimension(xlen, ylen, zlen) :: temp

  !! locals
  real(mytype), dimension(xlen, ylen, zlen) :: mweight

  temp(:,:,:) = pressure0 / rho(:,:,:)
  if (imultispecies) then
     call calc_mweight(mweight, phi, xlen, ylen, zlen)
     temp(:,:,:) = temp(:,:,:) * mweight(:,:,:)
  endif

endsubroutine calc_temp_eos
!-----------------------------------------------------------------------------!
subroutine calc_rho_eos(rho, temp, phi, mweight, xlen, ylen, zlen)

  use decomp_2d
  use param, only : pressure0, imultispecies
  use var, only : numscalar

  implicit none

  !! INPUTS
  integer, intent(in) :: xlen, ylen, zlen
  real(mytype), intent(in), dimension(xlen, ylen, zlen) :: temp
  real(mytype), intent(in), dimension(xlen, ylen, zlen, numscalar) :: phi

  !! OUTPUTS
  real(mytype), intent(out), dimension(xlen, ylen, zlen) :: rho

  !! LOCALS
  real(mytype), dimension(xlen, ylen, zlen) :: mweight

  rho(:,:,:) = pressure0 / temp(:,:,:)
  if (imultispecies) then
     call calc_mweight(mweight, phi, xlen, ylen, zlen)
     rho(:,:,:) = rho(:,:,:) * mweight(:,:,:)
  endif

endsubroutine calc_rho_eos
!-----------------------------------------------------------------------------!
subroutine calc_mweight(mweight, phi, xlen, ylen, zlen)

  use decomp_2d
  use param, only : zero, one
  use param, only : massfrac, mol_weight
  use var, only : numscalar

  implicit none

  integer, intent(in) :: xlen, ylen, zlen
  real(mytype), intent(in), dimension(xlen, ylen, zlen, numscalar) :: phi

  !! LOCALS
  real(mytype), dimension(xlen, ylen, zlen) :: mweight
  integer :: is

  mweight(:,:,:) = zero
  do is = 1, numscalar
     if (massfrac(is)) then
        mweight(:,:,:) = mweight(:,:,:) + phi(:,:,:,is) / mol_weight(is)
     endif
  enddo
  mweight(:,:,:) = one / mweight(:,:,:)

endsubroutine calc_mweight
!-----------------------------------------------------------------------------!
!  R8_RANDOM returns a pseudorandom number between 0 and 1.
!
!  Discussion:
!
!    This function returns a pseudo-random number rectangularly distributed
!    between 0 and 1.   The cycle length is 6.95E+12.  (See page 123
!    of Applied Statistics (1984) volume 33), not as claimed in the
!    original article.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2008
!
!  Author:
!
!    FORTRAN77 original version by Brian Wichman, David Hill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Brian Wichman, David Hill,
!    Algorithm AS 183: An Efficient and Portable Pseudo-Random
!    Number Generator,
!    Applied Statistics,
!    Volume 31, Number 2, 1982, pages 188-190.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) S1, S2, S3, three values used as the
!    seed for the sequence.  These values should be positive
!    integers between 1 and 30,000.
!
!    Output, real ( kind = 8 ) R8_RANDOM, the next value in the sequence.
!-----------------------------------------------------------------------------!
function r8_random ( s1, s2, s3 )

  implicit none

  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) s3
  real ( kind = 8 ) r8_random

  s1 = mod ( 171 * s1, 30269 )
  s2 = mod ( 172 * s2, 30307 )
  s3 = mod ( 170 * s3, 30323 )

  r8_random = mod ( real ( s1, kind = 8 ) / 30269.0D+00 &
                  + real ( s2, kind = 8 ) / 30307.0D+00 &
                  + real ( s3, kind = 8 ) / 30323.0D+00, 1.0D+00 )

  return
end
!-----------------------------------------------------------------------------!
function return_30k(x) result(y)

  integer ( kind = 4 ), intent(in) :: x
  integer ( kind = 4 )             :: y
  integer ( kind = 4 ), parameter  :: xmax = 30000

  y = iabs(x) - int(iabs(x)/xmax)*xmax
end function return_30k
!-----------------------------------------------------------------------------!
!  R8_UNI returns a pseudorandom number between 0 and 1.
!
!  Discussion:
!
!    This function generates uniformly distributed pseudorandom numbers
!    between 0 and 1, using the 32-bit generator from figure 3 of
!    the article by L'Ecuyer.
!
!    The cycle length is claimed to be 2.30584E+18.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 July 2008
!
!  Author:
!
!    Original Pascal original version by Pierre L'Ecuyer
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    Pierre LEcuyer,
!    Efficient and Portable Combined Random Number Generators,
!    Communications of the ACM,
!    Volume 31, Number 6, June 1988, pages 742-751.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) S1, S2, two values used as the
!    seed for the sequence.  On first call, the user should initialize
!    S1 to a value between 1 and 2147483562;  S2 should be initialized
!    to a value between 1 and 2147483398.
!
!    Output, real ( kind = 8 ) R8_UNI, the next value in the sequence.
!-----------------------------------------------------------------------------!
function r8_uni ( s1, s2 )

  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uni
  integer ( kind = 4 ) s1
  integer ( kind = 4 ) s2
  integer ( kind = 4 ) z

  k = s1 / 53668
  s1 = 40014 * ( s1 - k * 53668 ) - k * 12211
  if ( s1 < 0 ) then
    s1 = s1 + 2147483563
  end if

  k = s2 / 52774
  s2 = 40692 * ( s2 - k * 52774 ) - k * 3791
  if ( s2 < 0 ) then
    s2 = s2 + 2147483399
  end if

  z = s1 - s2
  if ( z < 1 ) then
    z = z + 2147483562
  end if

  r8_uni = real ( z, kind = 8 ) / 2147483563.0D+00

  return
end
!-----------------------------------------------------------------------------!
subroutine test_min_max(name,text,array_tmp,i_size_array_tmp)

  use param
  use variables
  use decomp_2d
  use MPI

  implicit none

  integer :: ierror, i, i_size_array_tmp
  real(mytype) :: max_tmp, min_tmp, tot_tmp, max_tot, min_tot, tot_tot
  real(mytype), dimension(i_size_array_tmp) :: array_tmp
  character(len=5) :: name
  character(len=15) :: text

  max_tmp=-0.000000000000000001_mytype
  tot_tmp=0._mytype
  min_tmp=+1000000000000000000._mytype
  do i=1,size(array_tmp)
    max_tmp=max(max_tmp,array_tmp(i))
    tot_tmp=tot_tmp + array_tmp(i)
    min_tmp=min(min_tmp,array_tmp(i))
  enddo
  call MPI_ALLREDUCE(max_tmp,max_tot,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
  call MPI_ALLREDUCE(min_tmp,min_tot,1,real_type,MPI_MIN,MPI_COMM_WORLD,ierror)
  call MPI_ALLREDUCE(tot_tmp,tot_tot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
  if (nrank == 0) then
     write(*,*) " "
     write(*,*) trim(text)//' Max ',name,max_tot
     write(*,*) trim(text)//' Tot ',name,tot_tot
     write(*,*) trim(text)//' Min ',name,min_tot
     write(*,*) " "
     flush(6)
  endif

  return
end subroutine test_min_max
!-----------------------------------------------------------------------------!
