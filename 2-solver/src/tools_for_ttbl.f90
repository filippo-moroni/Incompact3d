  
  
module tools_for_ttbl

  implicit none
  
  private 
  
  public :: calculate_friction_coefficient

contains

  !---------------------------------------------------------------------------!
  ! Calculate skin friction coefficient at the bottom wall and shear velocity
  ! Adapted from visu_ttbl subroutine.
  !---------------------------------------------------------------------------!
  subroutine calculate_friction_coefficient(ux1,uz1)
      
    use param
    use variables
    use var
    
    use ibm_param,   only : ubcx,ubcy,ubcz
    use dbg_schemes, only : sqrt_prec
    
    use MPI
    use decomp_2d
    use decomp_2d_io
    
    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uz1
    
    real(mytype) :: mean_gw      ! mean gradient at the wall at each processor
    real(mytype) :: mean_gw_tot  ! mean gradient at the wall total
    
    integer :: ierr  ! for MPI (initialized in init_xcompact3d subroutine)
    integer :: i,k
    
    ! Set again variables to zero
    mean_gw     = zero
    mean_gw_tot = zero
    fric_coeff  = zero
    sh_vel      = zero
    
    ! Perform communications if needed
    if (sync_vel_needed) then
      call transpose_x_to_y(ux1,ux2)
      call transpose_x_to_y(uz1,uz2)
      sync_vel_needed = .false.
    endif

    ! y-derivatives
    call dery (ta2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcx)
    call dery (tc2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1,ubcz)
    
    ! du/dy=ta2   
    ! dw/dy=tc2
    
    ! Mean velocity gradient at the wall, sqrt(du/dy**2 + dw/dy**2) and summation over all points
    do k=ystart(3),yend(3)
       do i=ystart(1),yend(1)
           
             ! Index for j is 1, since we are in global coordinates (y-pencils)
             mean_gw = mean_gw + (sqrt_prec(ta2(i,1,k)**2 + tc2(i,1,k)**2)) / real(nx*nz,mytype)                  
       enddo
    enddo
    
    ! Summation over all MPI processes and broadcast the result
    call MPI_ALLREDUCE(mean_gw,mean_gw_tot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    
    ! Calculate cf and shear velocity from the mean gradient at the wall    
    fric_coeff = mean_gw_tot * two * xnu / (uwall**2)
    sh_vel     = sqrt_prec(xnu * mean_gw_tot)
        
    return 
  
  end subroutine calculate_friction_coefficient

end module tools_for_ttbl

