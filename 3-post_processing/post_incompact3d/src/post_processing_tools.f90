
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
  integer,dimension(5) :: sel
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
  if (sel(2)==1) post_vort   = .true.
  if (sel(3)==1) post_diss   = .true.
  if (sel(4)==1) post_corz   = .true.
  if (sel(5)==1) post_tke_eq = .true.
  
  ! Return an error if all switchers are set to zero
  if (nrank==0) then
     if ((.not.post_mean   ) .and. &
         (.not.post_vort   ) .and. &
         (.not.post_diss   ) .and. &
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
  if (post_vort .or. post_diss .or. post_corz .or. post_tke_eq) read_vel=.true.
  
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
! DESCRIPTION: This subroutine is used to write the header of mean statistics
!              file 'mean_stats.txt', generated by 'post_incompact3d'.   
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine mean_stats_header(iunit)

  use param

  implicit none
  
  integer, intent(in) :: iunit
  
  ! TTBL
  if (itype .eq. itype_ttbl) then
  
      write(iunit, '(A)')      'Mean statistics for a Temporal Turbulent Boundary Layer (TTBL),'
      write(iunit, '(A)')      'with initialisation as Kozul et al. (2016).'
      write(iunit, '(A)')      ' '      
      write(iunit, '(A)')      'Simulation details:'      
      write(iunit, '(A,F5.2)') ' - Trip Reynolds number, Re_D = ', re      
      write(iunit, '(A)')      ' '      
      write(iunit, '(A)')      'Domain dimensions:'      
      write(iunit, '(A,F5.2)') ' - Lx = ', xlx         
  end if
  
end subroutine mean_stats_header











