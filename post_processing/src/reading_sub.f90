subroutine reading_input_file()

  USE decomp_2d
  USE decomp_2d_io
  USE variables
  USE param
  USE var
  USE MPI
  
  implicit none
  
  ! Variables to read the input.i3d file
  integer :: nargin, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase
  
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
  
end subroutine reading_input_file
