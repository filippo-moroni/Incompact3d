
!This file is not part of standard Xcompact3d releases (xcompact3d.com).

!-----------------------------------------------------------------------------!
! DESCRIPTION: This file is used to store useful subroutines for 
!              post-processing program 'post_incompact3d'.   
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
! DESCRIPTION: This subroutine is used to read the post-processing file
!              'post.prm' and to continue the initialization of 
!              post-processing work variables.   
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine read_post_file()

  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
  use decomp_2d_io
  
  use variables
  use param
  use var

  use post_processing

  implicit none

  ! Index for the number of post-processing selections employed (selector index)
  integer,dimension(4) :: sel
  integer              :: i
                    
  ! Reading of the input file of post-processing ('post.prm')
  open(newunit=iunit,file='post.prm',status='old',form='formatted')
    
     read (iunit,'(A1)') a
     read (iunit,'(A1)') a
     read (iunit,'(A1)') a 
     read (iunit,'(A1)') a
     read (iunit,'(A1)') a
     read (iunit,*) file1
     read (iunit,*) filen
     read (iunit,*) icrfile
     read (iunit,*) nr 
     read (iunit,'(A1)') a
     read (iunit,'(A1)') a
     read (iunit,'(A1)') a
     
     do i=1,size(sel)
        read (iunit,*) sel(i)
     enddo
     
  close(iunit)
  
  ! Set to 'true' logicals if they have been set to 1 in 'post.prm' file
  if (sel(1)==1) post_mean   = .true.
  if (sel(2)==1) post_grad   = .true.
  if (sel(3)==1) post_corz   = .true.
  if (sel(4)==1) post_tke_eq = .true.
  
  ! Return an error if all switchers are set to zero
  if (nrank==0) then
     if ((.not.post_mean   ) .and. &
         (.not.post_grad   ) .and. &
         (.not.post_corz   ) .and. &
         (.not.post_tke_eq))       &
        call decomp_2d_abort(iunit,'Invalid post-processing switchers specified, no work to be done here!')
  endif
  
  ! Logicals for reading snapshots
  if (post_mean) then
     read_vel=.true.
     read_pre=.true.
  endif
  
  ! Reading of velocity only if necessary
  if (post_grad .or. post_corz .or. post_tke_eq) read_vel=.true.
  
  ! Read of scalar field only if necessary
  if (iscalar==1) read_phi=.true. 
  
  ! Total number of snapshots in time
  nt = (filen-file1)/icrfile+1                         

end subroutine read_post_file

!-----------------------------------------------------------------------------!
! DESCRIPTION: This subroutine is used to read the 'time' attribute in a
!              .xdmf file.   
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine read_xdmf_time(filename, time_value, time_found)

  implicit none

  character(len=*),  intent(in)  :: filename    ! Input filename of the .xdmf file
  character(len=20), intent(out) :: time_value  ! Output variable for the extracted time value
  logical, intent(out) :: time_found            ! Flag to indicate if the time value was found
    
  character(len=256) :: line                    ! Buffer for reading lines from the file
  real    :: time_real                          ! Time value directly as real
  real    :: rounded_time                       ! To store the rounded time as real for formatting
  integer :: start_pos, end_pos                 ! To find the position of the value in the string
  integer :: ios                                ! I/O status variable for reading the file
  integer :: iunit                              ! Integer for fileunit
  
  ! Initialize output variables
  time_found = .false.
  time_value = ''

  ! Open the file (assuming newunit is supported by your compiler)
  open(newunit=iunit, file=filename, status='old', action='read')

  ! Loop through the file to find the <Time Value="..."> tag
  do while (.true.)
  
    ! Read one line at a time
    read(iunit, '(A)', iostat=ios) line  

    ! Exit loop if we reach the end of the file
    if (ios /= 0) exit                   

    if (index(line, '<Time Value=') > 0) then
        
        ! Find the position of the first and second quotes around the time value
        start_pos = index(line, '"') + 1
        end_pos = index(line(start_pos:), '"') + start_pos - 1

        ! Extract the time value directly and convert to real
        read(line(start_pos:end_pos), *) time_real

        ! Round the real number without converting to integer
        rounded_time = nint(time_real * 100.0) / 100.0  ! Keep 2 decimal places

        ! Format the rounded time to 2 decimal places as a string
        write(time_value, '(F6.2)') rounded_time        ! Format with 2 decimal places

        ! Set flag indicating that the time value was found
        time_found = .true.  
        
        ! Exit the loop once time is found
        exit  
    endif

  end do

  ! Close the file
  close(iunit)

end subroutine read_xdmf_time

!-----------------------------------------------------------------------------!
! DESCRIPTION: This subroutine is used to write the header of statistics
!              generated by 'post_incompact3d', depending on the flow case
!              and type of statistics we are saving.   
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine stats_header(iunit)

  use param
  use variables
  use post_processing

  implicit none
  
  integer, intent(in) :: iunit
  
  ! Start to write to .txt file with an empty row
  write(iunit, '(A)') ' '
  
  if (post_mean) then
  
        write(iunit, '(A)') ' Mean statistics for velocity, pressure and scalar fields.'
        write(iunit, '(A)') ' '       
  
  end if
  
  if (post_grad) then
  
        write(iunit, '(A)') ' Mean statistics for gradients-related quantities (vorticity, mean gradients, total dissipation).'
        write(iunit, '(A)') ' '       
  
  end if       
    
  ! TTBL
  if (itype .eq. itype_ttbl) then
  
      write(iunit, '(A)')      ' Temporal Turbulent Boundary Layer (TTBL),'
      write(iunit, '(A)')      ' with initialisation as Kozul et al. (2016).'
      write(iunit, '(A)')      ' '
      write(iunit, '(A)')      ' Statistics are adimensionalised with wall velocity Uw and trip wire diameter D, both unitary.'
      write(iunit, '(A)')      ' '       
      write(iunit, '(A)')      ' Simulation & post-processing details:'      
      write(iunit, '(A,F8.2)') '  - Trip Reynolds number,      Re_D = ', re
      write(iunit, '(A,I4)')   '  - Number of flow realizations, nr = ', nr 
  
  ! Channel
  else if (itype .eq. itype_channel) then
      
      write(iunit, '(A)')      ' Turbulent Channel.'
      write(iunit, '(A)')      ' '
      write(iunit, '(A)')      ' Statistics are adimensionalised with centerline velocity of the related laminar Poiseuille flow'
      write(iunit, '(A)')      ' and channel half-height, both unitary.'
      write(iunit, '(A)')      ' '       
      write(iunit, '(A)')      ' Simulation & post-processing details:'      
      write(iunit, '(A,F8.2)') '  - Centerline Reynolds number of the related laminar Poiseuille flow, Re_0 = ', re
      write(iunit, '(A,I4)')   '  - Number of snapshots in time, nt = ', nt  
  
  end if    
  
  write(iunit, '(A)')          ' '      
  write(iunit, '(A)')          ' Domain dimensions:'      
  write(iunit, '(A,F6.2)')     '  - Lx = ', xlx
  write(iunit, '(A,F6.2)')     '  - Ly = ', yly                  
  write(iunit, '(A,F6.2)')     '  - Lz = ', zlz
  write(iunit, '(A)')          ' '               
  write(iunit, '(A)')          ' Number of points:'      
  write(iunit, '(A,I4)')       '  - nx = ', nx
  write(iunit, '(A,I4)')       '  - ny = ', ny                  
  write(iunit, '(A,I4)')       '  - nz = ', nz
  write(iunit, '(A)')          ' '                      
  write(iunit, '(A)')          ' Nomenclature:'
  write(iunit, '(A)')          ' x          : streamwise  direction;'
  write(iunit, '(A)')          ' y          : wall-normal  direction;'  
  write(iunit, '(A)')          ' z          : spanwise  direction;'

  if (post_mean) then

      write(iunit, '(A)')      ' u          : streamwise  velocity;'
      write(iunit, '(A)')      ' v          : wall-normal velocity;'      
      write(iunit, '(A)')      ' w          : spanwise    velocity;'
      write(iunit, '(A)')      ' p          : pressure;'
      write(iunit, '(A)')      ' phi        : scalar field;'
      write(iunit, '(A)')      ' < >        : average operator;'
      write(iunit, '(A)')      " '          : fluctuating quantity;"      
      write(iunit, '(A)')      ' mean       : mean/average;'
      write(iunit, '(A)')      ' var        : variance;'
      write(iunit, '(A)')      ' skew       : skewness;'  
      write(iunit, '(A)')      ' kurt       : kurtosis.'
      
  end if
  
  if (post_grad) then

      write(iunit, '(A)')      ' omega_x    : streamwise vorticity;'
      write(iunit, '(A)')      ' omega_y    : wall-normal vorticity;'
      write(iunit, '(A)')      ' omega_z    : spanwise vorticity;'
      write(iunit, '(A)')      ' dU/dy      : mean streamwise velocity gradient;'
      write(iunit, '(A)')      ' dW/dy      : mean spanwise velocity gradient;'
      write(iunit, '(A)')      ' dPhi/dy    : mean scalar gradient;'
      write(iunit, '(A)')      ' eps        : total dissipation rate of kinetic energy;'                  
      write(iunit, '(A)')      ' mean       : mean/average;'

  end if

  if (post_tke_eq) then

      write(iunit, '(A)')      ' tke_turbt  : wall-normal turbulent transport of TKE;'
      write(iunit, '(A)')      ' tke_presst : wall-normal pressure transport of TKE;'
      write(iunit, '(A)')      ' tke_difft  : viscous diffusion of TKE;'
      write(iunit, '(A)')      ' tke_prod   : production of TKE;'
      write(iunit, '(A)')      ' tke_pseps  : pseudo-dissipation of TKE;'

end if
  
  ! Empty row
  write(iunit, '(A)') ' '  

  ! Write time and date to the .txt file
  call write_time_and_date(iunit) 
  
#ifdef TTBL_MODE  
        ! Add the time unit to the header
        write(iunit, '(A, A)') 'Time unit, t = ', time_value
        write(iunit, *) ' '
#endif 
  
end subroutine stats_header

!-----------------------------------------------------------------------------!
! DESCRIPTION: This subroutine is used to write to .txt files time and date.   
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine write_time_and_date(iunit)

  implicit none

  integer, intent(in) :: iunit
  
  ! Work variables to write the date and time to .txt files
  character(len=8)  :: date_str
  character(len=10) :: time_str
  integer :: year, month, day, hour, minute, second

  ! Get the current date and time with intrinsic subroutine of Fortran 90 and later versions
  call date_and_time(date=date_str, time=time_str)

  ! Extract year, month, day from the date string
  read(date_str, '(I4, I2, I2)') year, month, day

  ! Extract hour, minute, second from the time string
  read(time_str, '(I2, I2, I2)') hour, minute, second
  
  ! Write to .txt file
  write(iunit, '(A)')      ' Date and time of creation:'
  write(iunit, '(A, I2.2, A, I2.2, A, I4  )') " - Date (day/month/year) : ", day, "/", month, "/", year
  write(iunit, '(A, I2.2, A, I2.2, A, I2.2)') " - Time (hour/min/sec)   : ", hour, ":", minute, ":", second  
  write(iunit, '(A)')      ' '
  write(iunit, '(A)')      ' ' 

end subroutine write_time_and_date











