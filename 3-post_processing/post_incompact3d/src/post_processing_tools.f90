
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
    
  character(len=*),  intent(in)  :: filename    ! Input filename of the .xdmf file, len=* is a dummy variable
  character(len=20), intent(out) :: time_value  ! Output variable for the extracted time value
  logical, intent(out) :: time_found            ! Flag to indicate if the time value was found
    
  character(len=256) :: line                    ! Buffer for reading lines from the file
    
  integer :: ios    ! I/O status variable
  integer :: iunit  ! Unit number for file handling
  integer :: start_pos, end_pos  ! To handle the position of quotes for extracting the time value
  
  ! Declare a real variable to hold the time value for rounding
  real :: time_real

  ! Initialize output variables
  time_found = .false.  ! Set time_found to false initially
  time_value = ''       ! Clear the time_value string

  ! Open the .xdmf file for reading using newunit
  open(newunit=iunit, file=filename, status='old', action='read')

  ! Read through the file line by line
  do while (.true.)

      ! Read a line from the file
      read(iunit, '(A)', iostat=ios) line
      
      ! Exit the loop if end of file
      if (ios /= 0) exit

      ! Check if the line contains the <Time Value="..."> tag
      if (index(line, '<Time Value=') > 0) then
          ! Find the position of the first quote after the '=' sign
          start_pos = index(line, '"') + 1
          end_pos = index(line(start_pos:), '"') + start_pos - 1
          
          ! Extract the time value
          time_value = line(start_pos:end_pos)
          
          ! Round the time value
          time_real = round(time_real)

          ! Convert the rounded value back to a string
          write(time_value, '(F6.2)') time_real  ! Example format with 2 decimal places
    
          time_found = .true.
          exit  ! Exit the loop once the time is found
      endif

  end do

  ! Close the file
  close(iunit)

end subroutine read_xdmf_time


