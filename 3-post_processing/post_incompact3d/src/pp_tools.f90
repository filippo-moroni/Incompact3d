
!This file is not part of standard Xcompact3d releases (xcompact3d.com).

!-----------------------------------------------------------------------------!
! DESCRIPTION: This file is used to store useful subroutines for 
!              post-processing (pp) program 'post_incompact3d'.   
!   AUTHOR(s): Filippo Moroni <filippo.moroni@unimore.it> 
!-----------------------------------------------------------------------------!


subroutine read_xdmf_time(filename, time_value, time_found)
    implicit none
    character(len=*), intent(in) :: filename
    character(len=20), intent(out) :: time_value
    logical, intent(out) :: time_found
    character(len=256) :: line
    integer :: ios

    ! Initialize the output
    time_found = .false.
    time_value = ''

    ! Open the XDMF file
    open(unit=10, file=filename, status='old', action='read')

    ! Read through the file until we find the <Time> line
    do while (.true.)
        read(10, '(A)', iostat=ios) line
        if (ios /= 0) exit  ! Exit loop on end of file

        ! Check if the line contains the <Time Value="> element
        if (index(line, '<Time Value="') > 0) then
            ! Extract the time value
            read(line, '(A)', iostat=ios) time_value
            if (ios == 0) then
                ! Remove the <Time Value=" and closing quote
                time_value = trim(adjustl(time_value))
                time_value = trim(replace(time_value, '<Time Value="', ''))
                time_value = trim(replace(time_value, '" />', ''))
                time_found = .true.
                exit
            end if
        end if
    end do

    close(10)

contains

    function replace(str, old, new) result(res)
        character(len=*), intent(in) :: str, old, new
        character(len=len(str)) :: res
        integer :: pos

        res = str
        pos = index(res, old)
        if (pos > 0) then
            res = trim(adjustl(res(1:pos-1))) // new // trim(adjustl(res(pos + len(old):)))
        end if
    end function replace

end subroutine read_xdmf_time
