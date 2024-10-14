
!This file is not part of standard Xcompact3d releases (xcompact3d.com).

!-----------------------------------------------------------------------------!
! DESCRIPTION: This file is used to store useful subroutines for 
!              post-processing (pp) program 'post_incompact3d'.   
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!

!-----------------------------------------------------------------------------!
! DESCRIPTION: This subroutine is used to read the post-processing file
!              'post.prm'.   
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!
subroutine read_post_file

  


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
    
  integer :: ios    ! I/O status variable
  integer :: iunit  ! Unit number for file handling

  ! Initialize output variables
  time_found = .false.  ! Set time_found to false initially
  time_value = ''       ! Clear the time_value string

  ! Open the .xdmf file for reading using newunit
  open(newunit=iunit, file=filename, status='old', action='read')

  ! Read through the file line by line
  do while (.true.)

      ! Read a line from the file
      read(iunit, '(A)', iostat=ios) line
        
      ! Exit loop if end of file is reached
      if (ios /= 0) exit                  

      ! Check if the line contains the <Time Value="> tag
      if (index(line, '<Time Value="') > 0) then
            
          ! Extract the time value from the line
          read(line, '(A)', iostat=ios) time_value  ! Read the line into time_value
            
          ! Clean up the time_value to extract only the numerical value
          time_value = trim(adjustl(time_value))                       ! Trim leading spaces
          time_value = trim(replace(time_value, '<Time Value="', ''))  ! Remove opening tag
          time_value = trim(replace(time_value, '" />', ''))           ! Remove closing tag
            
          ! Set flag to true indicating time was found
          time_found = .true.  
            
          ! Exit the loop after finding the time
          exit  
      end if
  end do

  ! Close the file
  close(iunit)  

contains

  ! Function to replace a substring in a string
  function replace(str, old, new) result(res)
        
    character(len=*), intent(in) :: str, old, new  ! Input strings
    character(len=len(str)) :: res                 ! Resulting string
    integer :: pos                                 ! Position of the old substring

    ! Initialize result to the original string
    res = str

    ! Find position of the old substring
    pos = index(res, old)

    if (pos > 0) then

        ! Concatenate parts of the string to perform replacement
        res = trim(adjustl(res(1:pos-1))) // new // trim(adjustl(res(pos + len(old):)))
        
    end if
    
  end function replace

end subroutine read_xdmf_time

