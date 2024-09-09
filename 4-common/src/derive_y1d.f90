
!This file is not part of standard Xcompact3d releases (xcompact3d.com).

!---------------------------------------------------------------------------!
! DESCRIPTION: This file collects subroutines for derivatives in y-direction        
!              of arrays of size 'ysize(2)'. Useful for derivatives of 
!              averaged quantities, e.g. U = U (y) (mean velocity profile). 
!              Adapted from default derivative subroutines of 'derive.f90' 
!              of standard Incompact3d. 
!   AUTHOR(s): Roberto Corsini <roberto.corsini@unimore.it> 
!---------------------------------------------------------------------------!

subroutine dery1D_00(ty,uy,ry,sy,ffy,fsy,fwy,ppy,ny,npaire)

  use param
  use derivY

  implicit none

  integer :: ny,npaire,j
  real(mytype), dimension(ny) :: ty,uy
  real(mytype), dimension(ny) :: ry
  real(mytype) :: sy
  real(mytype), dimension(ny) :: ffy,fsy,fwy,ppy

  ty(1)=afjy*(uy(2)-uy(ny))&
       +bfjy*(uy(3)-uy(ny-1))
  ry(1)=-one
  ty(2)=afjy*(uy(3)-uy(1))&
       +bfjy*(uy(4)-uy(ny))
  ry(2)=zero
  
  do j=3,ny-2
     ty(j)=afjy*(uy(j+1)-uy(j-1))&
          +bfjy*(uy(j+2)-uy(j-2))
     ry(j)=zero
  enddo
  
  ty(ny-1)=afjy*(uy(ny)-uy(ny-2))&
       +bfjy*(uy(1)-uy(ny-3))
  ry(ny-1)=zero
  ty(ny)=afjy*(uy(1)-uy(ny-1))&
       +bfjy*(uy(2)-uy(ny-2))
  ry(ny)=alfajy
  
  do j=2,ny
     ty(j)=ty(j)-ty(j-1)*fsy(j)
     ry(j)=ry(j)-ry(j-1)*fsy(j)
  enddo
  
  ty(ny)=ty(ny)*fwy(ny)
  ry(ny)=ry(ny)*fwy(ny)
  
  do j=ny-1,1,-1
     ty(j)=(ty(j)-ffy(j)*ty(j+1))*fwy(j)
     ry(j)=(ry(j)-ffy(j)*ry(j+1))*fwy(j)
  enddo
  
  sy=(ty(1)-alfajy*ty(ny))&
         /(one+ry(1)-alfajy*ry(ny))
  
  do j=1,ny
     ty(j)=ty(j)-sy*ry(j)
  enddo
  if (istret.ne.0) then
     do j=1,ny
        ty(j)=ty(j)*ppy(j)
     enddo
  endif
  return
end subroutine dery1D_00

!---------------------------------------------------------------------------!
subroutine dery1D_11(ty,uy,ry,sy,ffy,fsy,fwy,ppy,ny,npaire) 

  use param
  use derivY

  implicit none

  integer :: ny,j,npaire
  real(mytype), dimension(ny) :: ty,uy
  real(mytype), dimension(ny) :: ry
  real(mytype) :: sy
  real(mytype), dimension(ny) :: ffy,fsy,fwy,ppy

  if (npaire==1) then
        
     ty(1)=zero
     ty(2)=afjy*(uy(3)-uy(1))&
       +bfjy*(uy(4)-uy(2))
     
     do j=3,ny-2
        ty(j)=afjy*(uy(j+1)-uy(j-1))&
             +bfjy*(uy(j+2)-uy(j-2))
     enddo
        
     ty(ny-1)=afjy*(uy(ny)-uy(ny-2))&
          +bfjy*(uy(ny-1)-uy(ny-3))
     ty(ny)=zero
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*fsy(j)
     enddo
     
     ty(ny)=ty(ny)*fwy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-ffy(j)*ty(j+1))*fwy(j)
     enddo
  endif
  if (npaire==0) then
        
     ty(1)=afjy*(uy(2)+uy(2))&
          +bfjy*(uy(3)+uy(3))
     ty(2)=afjy*(uy(3)-uy(1))&
          +bfjy*(uy(4)+uy(2))
     
     do j=3,ny-2
        ty(j)=afjy*(uy(j+1)-uy(j-1))&
             +bfjy*(uy(j+2)-uy(j-2))
     enddo
     
     ty(ny-1)=afjy*(uy(ny)-uy(ny-2))&
          +bfjy*((-uy(ny-1))-uy(ny-3))
     ty(ny)=afjy*((-uy(ny-1))-uy(ny-1))&
          +bfjy*((-uy(ny-2))-uy(ny-2))
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*fsy(j)
     enddo
     
     ty(ny)=ty(ny)*fwy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-ffy(j)*ty(j+1))*fwy(j)
     enddo
  endif
  if (istret.ne.0) then
     
     do j=1,ny
        ty(j)=ty(j)*ppy(j)
     enddo
  endif
  return
end subroutine dery1D_11

!---------------------------------------------------------------------------!
subroutine dery1D_12(ty,uy,ry,sy,ffy,fsy,fwy,ppy,ny,npaire) 

  use param
  use derivY

  implicit none

  integer :: ny,j,npaire
  real(mytype), dimension(ny) :: ty,uy
  real(mytype), dimension(ny) :: ry
  real(mytype) :: sy
  real(mytype), dimension(ny) :: ffy,fsy,fwy,ppy

  if (npaire==1) then
     
     ty(1)=zero
     ty(2)=afjy*(uy(3)-uy(1))&
          +bfjy*(uy(4)-uy(2))
     
     do j=3,ny-2
        ty(j)=afjy*(uy(j+1)-uy(j-1))&
             +bfjy*(uy(j+2)-uy(j-2))
     enddo
        
     ty(ny-1)=afmy*(uy(ny)-uy(ny-2))
     ty(ny)=-afny*uy(ny)-bfny*uy(ny-1)-cfny*uy(ny-2)
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*fsy(j)
     enddo
     
     ty(ny)=ty(ny)*fwy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-ffy(j)*ty(j+1))*fwy(j)
     enddo
  endif
  if (npaire==0) then
        
     ty(1)=afjy*(uy(2)+uy(2))&
          +bfjy*(uy(3)+uy(3))
     ty(2)=afjy*(uy(3)-uy(1))&
          +bfjy*(uy(4)+uy(2))
     
     do j=3,ny-2
        ty(j)=afjy*(uy(j+1)-uy(j-1))&
             +bfjy*(uy(j+2)-uy(j-2))
     enddo
     
     ty(ny-1)=afmy*(uy(ny)-uy(ny-2))
     ty(ny)=-afny*uy(ny)-bfny*uy(ny-1)-cfny*uy(ny-2)
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*fsy(j)
     enddo
     
     ty(ny)=ty(ny)*fwy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-ffy(j)*ty(j+1))*fwy(j)
     enddo
  endif
  if (istret.ne.0) then
     do j=1,ny
        ty(j)=ty(j)*ppy(j)
     enddo
  endif
  return
end subroutine dery1D_12

!---------------------------------------------------------------------------!
subroutine dery1D_21(ty,uy,ry,sy,ffy,fsy,fwy,ppy,ny,npaire) 

  use param
  use derivY

  implicit none

  integer :: ny,j,npaire
  real(mytype), dimension(ny) :: ty,uy
  real(mytype), dimension(ny) :: ry
  real(mytype) :: sy
  real(mytype), dimension(ny) :: ffy,fsy,fwy,ppy

  if (npaire==1) then
     
     ty(1)=af1y*uy(1)+bf1y*uy(2)+cf1y*uy(3)
     ty(2)=af2y*(uy(3)-uy(1))
     
     do j=3,ny-2
        ty(j)=afjy*(uy(j+1)-uy(j-1))&
             +bfjy*(uy(j+2)-uy(j-2))
     enddo
     
     ty(ny-1)=afjy*(uy(ny)-uy(ny-2))&
          +bfjy*(uy(ny-1)-uy(ny-3))
     ty(ny)=zero
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*fsy(j)
     enddo
     
     ty(ny)=ty(ny)*fwy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-ffy(j)*ty(j+1))*fwy(j)
     enddo
  endif
  if (npaire==0) then
        
     ty(1)=af1y*uy(1)+bf1y*uy(2)+cf1y*uy(3)
     ty(2)=af2y*(uy(3)-uy(1))
     
     do j=3,ny-2
        ty(j)=afjy*(uy(j+1)-uy(j-1))&
             +bfjy*(uy(j+2)-uy(j-2))
     enddo
     
     ty(ny-1)=afjy*(uy(ny)-uy(ny-2))&
          +bfjy*((-uy(ny-1))-uy(ny-3))
     ty(ny)=afjy*((-uy(ny-1))-uy(ny-1))&
          +bfjy*((-uy(ny-2))-uy(ny-2))
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*fsy(j)
     enddo
        
     ty(ny)=ty(ny)*fwy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-ffy(j)*ty(j+1))*fwy(j)
     enddo
  endif
  if (istret.ne.0) then
     do j=1,ny
        ty(j)=ty(j)*ppy(j)
     enddo
  endif
  return
end subroutine dery1D_21

!---------------------------------------------------------------------------!
subroutine dery1D_22(ty,uy,ry,sy,ffy,fsy,fwy,ppy,ny,npaire) 

  use param
  use derivY

  implicit none

  integer :: ny,j,npaire
  real(mytype), dimension(ny) :: ty,uy
  real(mytype), dimension(ny) :: ry
  real(mytype) :: sy
  real(mytype), dimension(ny) :: ffy,fsy,fwy,ppy
  
     
     ty(1)=af1y*uy(1)+bf1y*uy(2)+cf1y*uy(3)
     ty(2)=af2y*(uy(3)-uy(1))
  
     do j=3,ny-2
        ty(j)=afjy*(uy(j+1)-uy(j-1))&
             +bfjy*(uy(j+2)-uy(j-2))
     enddo
  
     ty(ny-1)=afmy*(uy(ny)-uy(ny-2))
     ty(ny)=-afny*uy(ny)-bfny*uy(ny-1)-cfny*uy(ny-2)
  
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*fsy(j)
     enddo
     
     ty(ny)=ty(ny)*fwy(ny)
  
     do j=ny-1,1,-1
        ty(j)=(ty(j)-ffy(j)*ty(j+1))*fwy(j)
     enddo

  if (istret.ne.0) then
     do j=1,ny
        ty(j)=ty(j)*ppy(j)
     enddo
  endif

  return
end subroutine dery1D_22

!---------------------------------------------------------------------------!
subroutine deryy1D_00(ty,uy,ry,sy,sfy,ssy,swy,ny,npaire) 

  use param
  use derivY

  implicit none

  integer :: ny,npaire,j
  real(mytype), dimension(ny) :: ty,uy,ry
  real(mytype) :: sy
  real(mytype), dimension(ny) :: sfy,ssy,swy
     
  ty(1)=asjy*(uy(2)-uy(1)&
       -uy(1)+uy(ny))&
       +bsjy*(uy(3)-uy(1)&
       -uy(1)+uy(ny-1))&
       +csjy*(uy(4)-uy(1)&
       -uy(1)+uy(ny-2))&
       +dsjy*(uy(5)-uy(1)&
       -uy(1)+uy(ny-3))
  ry(1)=-one
  ty(2)=asjy*(uy(3)-uy(2)&
       -uy(2)+uy(1))&
       +bsjy*(uy(4)-uy(2)&
       -uy(2)+uy(ny))&
       +csjy*(uy(5)-uy(2)&
       -uy(2)+uy(ny-1))&
       +dsjy*(uy(6)-uy(2)&
       -uy(2)+uy(ny-2))
  ry(2)=zero
  ty(3)=asjy*(uy(4)-uy(3)&
       -uy(3)+uy(2))&
       +bsjy*(uy(5)-uy(3)&
       -uy(3)+uy(1))&
       +csjy*(uy(6)-uy(3)&
       -uy(3)+uy(ny))&
       +dsjy*(uy(7)-uy(3)&
       -uy(3)+uy(ny-1))
  ry(3)=zero
  ty(4)=asjy*(uy(5)-uy(4)&
       -uy(4)+uy(3))&
       +bsjy*(uy(6)-uy(4)&
       -uy(4)+uy(2))&
       +csjy*(uy(7)-uy(4)&
       -uy(4)+uy(1))&
       +dsjy*(uy(8)-uy(4)&
       -uy(4)+uy(ny))
  ry(4)=zero
  
  do j=5,ny-4
     ty(j)=asjy*(uy(j+1)-uy(j)&
          -uy(j)+uy(j-1))&
          +bsjy*(uy(j+2)-uy(j)&
          -uy(j)+uy(j-2))&
          +csjy*(uy(j+3)-uy(j)&
          -uy(j)+uy(j-3))&
          +dsjy*(uy(j+4)-uy(j)&
          -uy(j)+uy(j-4))
     ry(j)=zero
  enddo
  
  ty(ny-3)=asjy*(uy(ny-2)-uy(ny-3)&
       -uy(ny-3)+uy(ny-4))&
       +bsjy*(uy(ny-1)-uy(ny-3)&
       -uy(ny-3)+uy(ny-5))&
       +csjy*(uy(ny)-uy(ny-3)&
       -uy(ny-3)+uy(ny-6))&
       +dsjy*(uy(1)-uy(ny-3)&
       -uy(ny-3)+uy(ny-7))
  ry(ny-3)=zero
  ty(ny-2)=asjy*(uy(ny-1)-uy(ny-2)&
       -uy(ny-2)+uy(ny-3))&
       +bsjy*(uy(ny)-uy(ny-2)&
       -uy(ny-2)+uy(ny-4))&
       +csjy*(uy(1)-uy(ny-2)&
       -uy(ny-2)+uy(ny-5))&
       +dsjy*(uy(2)-uy(ny-2)&
       -uy(ny-2)+uy(ny-6))
  ry(ny-2)=zero
  ty(ny-1)=asjy*(uy(ny)-uy(ny-1)&
       -uy(ny-1)+uy(ny-2))&
       +bsjy*(uy(1)-uy(ny-1)&
       -uy(ny-1)+uy(ny-3))&
       +csjy*(uy(2)-uy(ny-1)&
       -uy(ny-1)+uy(ny-4))&
       +dsjy*(uy(3)-uy(ny-1)&
       -uy(ny-1)+uy(ny-5))
  ry(ny-1)=zero
  ty(ny)=asjy*(uy(1)-uy(ny)&
       -uy(ny)+uy(ny-1))&
       +bsjy*(uy(2)-uy(ny)&
    -uy(ny)+uy(ny-2))&
          +csjy*(uy(3)-uy(ny)&
    -uy(ny)+uy(ny-3))&
       +dsjy*(uy(4)-uy(ny)&
       -uy(ny)+uy(ny-4))
  ry(ny)=alsajy
  if (iimplicit.ge.1) return
  
  do j=2,ny
     ty(j)=ty(j)-ty(j-1)*ssy(j)
     ry(j)=ry(j)-ry(j-1)*ssy(j)
  enddo
  
  ty(ny)=ty(ny)*swy(ny)
  ry(ny)=ry(ny)*swy(ny)
  
  do j=ny-1,1,-1
     ty(j)=(ty(j)-sfy(j)*ty(j+1))*swy(j)
     ry(j)=(ry(j)-sfy(j)*ry(j+1))*swy(j)
  enddo
  
  sy=(ty(1)-alsajy*ty(ny))/&
        (one+ry(1)-alsajy*ry(ny))
  
  do j=1,ny
     ty(j)=ty(j)-sy*ry(j)
  enddo
  return
end subroutine deryy1D_00

!---------------------------------------------------------------------------!
subroutine deryy1D_11(ty,uy,ry,sy,sfy,ssy,swy,ny,npaire) 

  use param
  use derivY

  implicit none

  integer :: ny,npaire,j
  real(mytype), dimension(ny) :: ty,uy,ry
  real(mytype) :: sy
  real(mytype), dimension(ny) :: sfy,ssy,swy

  if (npaire==1) then
     
     ty(1)=asjy*(uy(2)-uy(1)&
          -uy(1)+uy(2))&
          +bsjy*(uy(3)-uy(1)&
          -uy(1)+uy(3))&
          +csjy*(uy(4)-uy(1)&
          -uy(1)+uy(4))&
          +dsjy*(uy(5)-uy(1)&
          -uy(1)+uy(5))
     ty(2)=asjy*(uy(3)-uy(2)&
          -uy(2)+uy(1))&
          +bsjy*(uy(4)-uy(2)&
          -uy(2)+uy(2))&
          +csjy*(uy(5)-uy(2)&
          -uy(2)+uy(3))&
          +dsjy*(uy(6)-uy(2)&
          -uy(2)+uy(4))
     ty(3)=asjy*(uy(4)-uy(3)&
          -uy(3)+uy(2))&
          +bsjy*(uy(5)-uy(3)&
          -uy(3)+uy(1))&
          +csjy*(uy(6)-uy(3)&
          -uy(3)+uy(2))&
          +dsjy*(uy(7)-uy(3)&
          -uy(3)+uy(3))
     ty(4)=asjy*(uy(5)-uy(4)&
          -uy(4)+uy(3))&
          +bsjy*(uy(6)-uy(4)&
          -uy(4)+uy(2))&
          +csjy*(uy(7)-uy(4)&
          -uy(4)+uy(1))&
          +dsjy*(uy(8)-uy(4)&
          -uy(4)+uy(2))
     
     do j=5,ny-4
        ty(j)=asjy*(uy(j+1)-uy(j)&
             -uy(j)+uy(j-1))&
             +bsjy*(uy(j+2)-uy(j)&
             -uy(j)+uy(j-2))&
             +csjy*(uy(j+3)-uy(j)&
             -uy(j)+uy(j-3))&
             +dsjy*(uy(j+4)-uy(j)&
             -uy(j)+uy(j-4))
     enddo
     
     ty(ny-3)=asjy*(uy(ny-2)-uy(ny-3)&
          -uy(ny-3)+uy(ny-4))&
          +bsjy*(uy(ny-1)-uy(ny-3)&
          -uy(ny-3)+uy(ny-5))&
          +csjy*(uy(ny)-uy(ny-3)&
          -uy(ny-3)+uy(ny-6))&
          +dsjy*(uy(ny-1)-uy(ny-3)&
          -uy(ny-3)+uy(ny-7))
     ty(ny-2)=asjy*(uy(ny-1)-uy(ny-2)&
          -uy(ny-2)+uy(ny-3))&
          +bsjy*(uy(ny)-uy(ny-2)&
          -uy(ny-2)+uy(ny-4))&
          +csjy*(uy(ny-1)-uy(ny-2)&
          -uy(ny-2)+uy(ny-5))&
          +dsjy*(uy(ny-2)-uy(ny-2)&
          -uy(ny-2)+uy(ny-6))
     ty(ny-1)=asjy*(uy(ny)-uy(ny-1)&
          -uy(ny-1)+uy(ny-2))&
          +bsjy*(uy(ny-1)-uy(ny-1)&
          -uy(ny-1)+uy(ny-3))&
          +csjy*(uy(ny-2)-uy(ny-1)&
          -uy(ny-1)+uy(ny-4))&
          +dsjy*(uy(ny-3)-uy(ny-1)&
          -uy(ny-1)+uy(ny-5))
     ty(ny)=asjy*(uy(ny-1)-uy(ny)&
          -uy(ny)+uy(ny-1))&
          +bsjy*(uy(ny-2)-uy(ny)&
          -uy(ny)+uy(ny-2))&
          +csjy*(uy(ny-3)-uy(ny)&
          -uy(ny)+uy(ny-3))&
          +dsjy*(uy(ny-4)-uy(ny)&
          -uy(ny)+uy(ny-4))
     if (iimplicit.ge.1) return
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*ssy(j)
     enddo
     
     ty(ny)=ty(ny)*swy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-sfy(j)*ty(j+1))*swy(j)
     enddo
  endif
  if (npaire==0) then
     
     ty(1)=zero
     ty(2)=asjy*(uy(3)-uy(2)&
          -uy(2)+uy(1))&
          +bsjy*(uy(4)-uy(2)&
          -uy(2)-uy(2))&
          +csjy*(uy(5)-uy(2)&
          -uy(2)-uy(3))&
          +dsjy*(uy(6)-uy(2)&
          -uy(2)-uy(4))
     ty(3)=asjy*(uy(4)-uy(3)&
          -uy(3)+uy(2))&
          +bsjy*(uy(5)-uy(3)&
          -uy(3)+uy(1))&
          +csjy*(uy(6)-uy(3)&
          -uy(3)-uy(2))&
          +dsjy*(uy(7)-uy(3)&
          -uy(3)-uy(3))
     ty(4)=asjy*(uy(5)-uy(4)&
          -uy(4)+uy(3))&
          +bsjy*(uy(6)-uy(4)&
          -uy(4)+uy(2))&
          +csjy*(uy(7)-uy(4)&
          -uy(4)-uy(1))&
          +dsjy*(uy(8)-uy(4)&
          -uy(4)-uy(2))
 
     do j=5,ny-4
        ty(j)=asjy*(uy(j+1)-uy(j)&
          -uy(j)+uy(j-1))&
          +bsjy*(uy(j+2)-uy(j)&
          -uy(j)+uy(j-2))&
          +csjy*(uy(j+3)-uy(j)&
          -uy(j)+uy(j-3))&
          +dsjy*(uy(j+4)-uy(j)&
          -uy(j)+uy(j-4))
     enddo
           
     ty(ny-3)=asjy*( uy(ny-2)-uy(ny-3)&
          -uy(ny-3)+uy(ny-4))&
          +bsjy*( uy(ny-1)-uy(ny-3)&
          -uy(ny-3)+uy(ny-5))&
          +csjy*(-uy(ny)-uy(ny-3)&
          -uy(ny-3)+uy(ny-6))&
          +dsjy*(-uy(ny-1)-uy(ny-3)&
          -uy(ny-3)+uy(ny-7))
     ty(ny-2)=asjy*( uy(ny-1)-uy(ny-2)&
          -uy(ny-2)+uy(ny-3))&
          +bsjy*( uy(ny)-uy(ny-2)&
          -uy(ny-2)+uy(ny-4))&
          +csjy*(-uy(ny-1)-uy(ny-2)&
          -uy(ny-2)+uy(ny-5))&
          +dsjy*(-uy(ny-2)-uy(ny-2)&
          -uy(ny-2)+uy(ny-6))
     ty(ny-1)=asjy*( uy(ny)-uy(ny-1)&
          -uy(ny-1)+uy(ny-2))&
          +bsjy*(-uy(ny-1)-uy(ny-1)&
          -uy(ny-1)+uy(ny-3))&
          +csjy*(-uy(ny-2)-uy(ny-1)&
          -uy(ny-1)+uy(ny-4))&
          +dsjy*(-uy(ny-3)-uy(ny-1)&
          -uy(ny-1)+uy(ny-5))
     ty(ny)=zero
     if (iimplicit.ge.1) return
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*ssy(j)
     enddo
     
     ty(ny)=ty(ny)*swy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-sfy(j)*ty(j+1))*swy(j)
     enddo
  endif
  return
end subroutine deryy1D_11

!---------------------------------------------------------------------------!
subroutine deryy1D_12(ty,uy,ry,sy,sfy,ssy,swy,ny,npaire) 

  use param
  use derivY

  implicit none

  integer :: ny,npaire,j
  real(mytype), dimension(ny) :: ty,uy,ry
  real(mytype) :: sy
  real(mytype), dimension(ny) :: sfy,ssy,swy

  if (npaire==1) then
        
     ty(1)=asjy*(uy(2)-uy(1)&
          -uy(1)+uy(2))&
          +bsjy*(uy(3)-uy(1)&
          -uy(1)+uy(3))&
          +csjy*(uy(4)-uy(1)&
          -uy(1)+uy(4))&
          +dsjy*(uy(5)-uy(1)&
          -uy(1)+uy(5))
     ty(2)=asjy*(uy(3)-uy(2)&
          -uy(2)+uy(1))&
          +bsjy*(uy(4)-uy(2)&
          -uy(2)+uy(2))&
          +csjy*(uy(5)-uy(2)&
          -uy(2)+uy(3))&
          +dsjy*(uy(6)-uy(2)&
          -uy(2)+uy(4))
     ty(3)=asjy*(uy(4)-uy(3)&
          -uy(3)+uy(2))&
          +bsjy*(uy(5)-uy(3)&
          -uy(3)+uy(1))&
          +csjy*(uy(6)-uy(3)&
          -uy(3)+uy(2))&
          +dsjy*(uy(7)-uy(3)&
          -uy(3)+uy(3))
     ty(4)=asjy*(uy(5)-uy(4)&
          -uy(4)+uy(3))&
          +bsjy*(uy(6)-uy(4)&
          -uy(4)+uy(2))&
          +csjy*(uy(7)-uy(4)&
          -uy(4)+uy(1))&
          +dsjy*(uy(8)-uy(4)&
          -uy(4)+uy(2))
     
     do j=5,ny-4
        ty(j)=asjy*(uy(j+1)-uy(j)&
             -uy(j)+uy(j-1))&
             +bsjy*(uy(j+2)-uy(j)&
             -uy(j)+uy(j-2))&
             +csjy*(uy(j+3)-uy(j)&
             -uy(j)+uy(j-3))&
             +dsjy*(uy(j+4)-uy(j)&
             -uy(j)+uy(j-4))
     enddo
     
     ty(ny-3)=astty*(uy(ny-2)-uy(ny-3)&
          -uy(ny-3)+uy(ny-4))&
          +bstty*(uy(ny-1)-uy(ny-3)&
          -uy(ny-3)+uy(ny-5))&
          +cstty*(uy(ny)-uy(ny-3)&
          -uy(ny-3)+uy(ny-6))
     ty(ny-2)=asty*(uy(ny-1)-uy(ny-2)&
          -uy(ny-2)+uy(ny-3))&
          +bsty*(uy(ny)-uy(ny-2)&
          -uy(ny-2)+uy(ny-4))
     ty(ny-1)=asmy*(uy(ny)-uy(ny-1)&
          -uy(ny-1)+uy(ny-2))
     ty(ny)=asny*uy(ny)+bsny*uy(ny-1)&
                +csny*uy(ny-2)+dsny*uy(ny-3)
     if (iimplicit.ge.1) return
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*ssy(j)
     enddo
     
     ty(ny)=ty(ny)*swy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-sfy(j)*ty(j+1))*swy(j)
     enddo
  endif
  if (npaire==0) then
     
     ty(1)=zero
     ty(2)=asjy*(uy(3)-uy(2)&
          -uy(2)+uy(1))&
          +bsjy*(uy(4)-uy(2)&
          -uy(2)-uy(2))&
          +csjy*(uy(5)-uy(2)&
          -uy(2)-uy(3))&
          +dsjy*(uy(6)-uy(2)&
          -uy(2)-uy(4))
     ty(3)=asjy*(uy(4)-uy(3)&
          -uy(3)+uy(2))&
          +bsjy*(uy(5)-uy(3)&
          -uy(3)+uy(1))&
          +csjy*(uy(6)-uy(3)&
          -uy(3)-uy(2))&
          +dsjy*(uy(7)-uy(3)&
          -uy(3)-uy(3))
     ty(4)=asjy*(uy(5)-uy(4)&
          -uy(4)+uy(3))&
          +bsjy*(uy(6)-uy(4)&
          -uy(4)+uy(2))&
          +csjy*(uy(7)-uy(4)&
          -uy(4)-uy(1))&
          +dsjy*(uy(8)-uy(4)&
          -uy(4)-uy(2))
     
     do j=5,ny-4
        ty(j)=asjy*(uy(j+1)-uy(j)&
             -uy(j)+uy(j-1))&
             +bsjy*(uy(j+2)-uy(j)&
             -uy(j)+uy(j-2))&
             +csjy*(uy(j+3)-uy(j)&
             -uy(j)+uy(j-3))&
             +dsjy*(uy(j+4)-uy(j)&
             -uy(j)+uy(j-4))
     enddo
     
     ty(ny-3)=astty*(uy(ny-2)-uy(ny-3)&
          -uy(ny-3)+uy(ny-4))&
          +bstty*(uy(ny-1)-uy(ny-3)&
          -uy(ny-3)+uy(ny-5))&
          +cstty*(uy(ny)-uy(ny-3)&
          -uy(ny-3)+uy(ny-6))
     ty(ny-2)=asty*(uy(ny-1)-uy(ny-2)&
          -uy(ny-2)+uy(ny-3))&
          +bsty*(uy(ny)-uy(ny-2)&
          -uy(ny-2)+uy(ny-4))
     ty(ny-1)=asmy*(uy(ny)-uy(ny-1)&
          -uy(ny-1)+uy(ny-2))
     ty(ny)=asny*uy(ny)+bsny*uy(ny-1)&
          +csny*uy(ny-2)+dsny*uy(ny-3)
     if (iimplicit.ge.1) return
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*ssy(j)
     enddo
     
     ty(ny)=ty(ny)*swy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-sfy(j)*ty(j+1))*swy(j)
     enddo
  endif

  return
end subroutine deryy1D_12

!---------------------------------------------------------------------------!
subroutine deryy1D_21(ty,uy,ry,sy,sfy,ssy,swy,ny,npaire) 

  use param
  use derivY

  implicit none

  integer :: ny,npaire,j
  real(mytype), dimension(ny) :: ty,uy,ry
  real(mytype) :: sy
  real(mytype), dimension(ny) :: sfy,ssy,swy

  if (npaire==1) then
     
     ty(1)=as1y*uy(1)+bs1y*uy(2)&
          +cs1y*uy(3)+ds1y*uy(4)
     ty(2)=as2y*(uy(3)-uy(2)&
          -uy(2)+uy(1))
     ty(3)=as3y*(uy(4)-uy(3)&
          -uy(3)+uy(2))&
          +bs3y*(uy(5)-uy(3)&
          -uy(3)+uy(1))
     ty(4)=as4y*(uy(5)-uy(4)&
          -uy(4)+uy(3))&
          +bs4y*(uy(6)-uy(4)&
          -uy(4)+uy(2))&
          +cs4y*(uy(7)-uy(4)&
          -uy(4)+uy(1))
     
     do j=5,ny-4
        ty(j)=asjy*(uy(j+1)-uy(j)&
             -uy(j)+uy(j-1))&
             +bsjy*(uy(j+2)-uy(j)&
             -uy(j)+uy(j-2))&
             +csjy*(uy(j+3)-uy(j)&
             -uy(j)+uy(j-3))&
             +dsjy*(uy(j+4)-uy(j)&
             -uy(j)+uy(j-4))
     enddo
        
     ty(ny-3)=asjy*(uy(ny-2)-uy(ny-3)&
          -uy(ny-3)+uy(ny-4))&
          +bsjy*(uy(ny-1)-uy(ny-3)&
          -uy(ny-3)+uy(ny-5))&
          +csjy*(uy(ny)-uy(ny-3)&
          -uy(ny-3)+uy(ny-6))&
          +dsjy*(uy(ny-1)-uy(ny-3)&
          -uy(ny-3)+uy(ny-7))
     ty(ny-2)=asjy*(uy(ny-1)-uy(ny-2)&
          -uy(ny-2)+uy(ny-3))&
          +bsjy*(uy(ny)-uy(ny-2)&
          -uy(ny-2)+uy(ny-4))&
          +csjy*(uy(ny-1)-uy(ny-2)&
          -uy(ny-2)+uy(ny-5))&
          +dsjy*(uy(ny-2)-uy(ny-2)&
          -uy(ny-2)+uy(ny-6))
     ty(ny-1)=asjy*(uy(ny)-uy(ny-1)&
          -uy(ny-1)+uy(ny-2))&
          +bsjy*(uy(ny-1)-uy(ny-1)&
          -uy(ny-1)+uy(ny-3))&
          +csjy*(uy(ny-2)-uy(ny-1)&
          -uy(ny-1)+uy(ny-4))&
          +dsjy*(uy(ny-3)-uy(ny-1)&
          -uy(ny-1)+uy(ny-5))
     ty(ny)=asjy*(uy(ny-1)-uy(ny)&
          -uy(ny)+uy(ny-1))&
          +bsjy*(uy(ny-2)-uy(ny)&
          -uy(ny)+uy(ny-2))&
          +csjy*(uy(ny-3)-uy(ny)&
          -uy(ny)+uy(ny-3))&
          +dsjy*(uy(ny-4)-uy(ny)&
          -uy(ny)+uy(ny-4))
     if (iimplicit.ge.1) return
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*ssy(j)
     enddo
        
     ty(ny)=ty(ny)*swy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-sfy(j)*ty(j+1))*swy(j)
     enddo
  endif
  if (npaire==0) then
     
     ty(1)=as1y*uy(1)+bs1y*uy(2)&
          +cs1y*uy(3)+ds1y*uy(4)
     ty(2)=as2y*(uy(3)-uy(2)&
          -uy(2)+uy(1))
     ty(3)=as3y*(uy(4)-uy(3)&
          -uy(3)+uy(2))&
          +bs3y*(uy(5)-uy(3)&
          -uy(3)+uy(1))
     ty(4)=as4y*(uy(5)-uy(4)&
          -uy(4)+uy(3))&
          +bs4y*(uy(6)-uy(4)&
          -uy(4)+uy(2))&
          +cs4y*(uy(7)-uy(4)&
          -uy(4)+uy(1))
     
     do j=5,ny-4
        ty(j)=asjy*(uy(j+1)-uy(j)&
             -uy(j)+uy(j-1))&
             +bsjy*(uy(j+2)-uy(j)&
             -uy(j)+uy(j-2))&
             +csjy*(uy(j+3)-uy(j)&
             -uy(j)+uy(j-3))&
             +dsjy*(uy(j+4)-uy(j)&
             -uy(j)+uy(j-4))
     enddo
        
     ty(ny-3)=asjy*( uy(ny-2)-uy(ny-3)&
          -uy(ny-3)+uy(ny-4))&
          +bsjy*( uy(ny-1)-uy(ny-3)&
          -uy(ny-3)+uy(ny-5))&
          +csjy*(-uy(ny)-uy(ny-3)&
          -uy(ny-3)+uy(ny-6))&
          +dsjy*(-uy(ny-1)-uy(ny-3)&
          -uy(ny-3)+uy(ny-7))
     ty(ny-2)=asjy*( uy(ny-1)-uy(ny-2)&
          -uy(ny-2)+uy(ny-3))&
          +bsjy*( uy(ny)-uy(ny-2)&
          -uy(ny-2)+uy(ny-4))&
          +csjy*(-uy(ny-1)-uy(ny-2)&
          -uy(ny-2)+uy(ny-5))&
          +dsjy*(-uy(ny-2)-uy(ny-2)&
          -uy(ny-2)+uy(ny-6))
     ty(ny-1)=asjy*( uy(ny)-uy(ny-1)&
          -uy(ny-1)+uy(ny-2))&
          +bsjy*(-uy(ny-1)-uy(ny-1)&
          -uy(ny-1)+uy(ny-3))&
          +csjy*(-uy(ny-2)-uy(ny-1)&
          -uy(ny-1)+uy(ny-4))&
          +dsjy*(-uy(ny-3)-uy(ny-1)&
          -uy(ny-1)+uy(ny-5))
     ty(ny)=zero
     if (iimplicit.ge.1) return
     
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*ssy(j)
     enddo
     
     ty(ny)=ty(ny)*swy(ny)
     
     do j=ny-1,1,-1
        ty(j)=(ty(j)-sfy(j)*ty(j+1))*swy(j)
     enddo
  endif

  return
end subroutine deryy1D_21

!---------------------------------------------------------------------------!
subroutine deryy1D_22(ty,uy,ry,sy,sfy,ssy,swy,ny,npaire) 

  use param
  use derivY

  implicit none

  integer :: ny,npaire,j
  real(mytype), dimension(ny) :: ty,uy,ry
  real(mytype) :: sy
  real(mytype), dimension(ny) :: sfy,ssy,swy
     
  ty(1)=as1y*uy(1)+bs1y*uy(2)&
       +cs1y*uy(3)+ds1y*uy(4)
  ty(2)=as2y*(uy(3)-uy(2)&
       -uy(2)+uy(1))
  ty(3)=as3y*(uy(4)-uy(3)&
       -uy(3)+uy(2))&
       +bs3y*(uy(5)-uy(3)&
       -uy(3)+uy(1))
  ty(4)=as4y*(uy(5)-uy(4)&
       -uy(4)+uy(3))&
       +bs4y*(uy(6)-uy(4)&
       -uy(4)+uy(2))&
       +cs4y*(uy(7)-uy(4)&
       -uy(4)+uy(1))
  
  do j=5,ny-4
     ty(j)=asjy*(uy(j+1)-uy(j)&
         -uy(j)+uy(j-1))&
         +bsjy*(uy(j+2)-uy(j)&
         -uy(j)+uy(j-2))&
         +csjy*(uy(j+3)-uy(j)&
         -uy(j)+uy(j-3))&
         +dsjy*(uy(j+4)-uy(j)&
         -uy(j)+uy(j-4))
  enddo
  
  ty(ny-3)=astty*(uy(ny-2)-uy(ny-3)&
       -uy(ny-3)+uy(ny-4))&
       +bstty*(uy(ny-1)-uy(ny-3)&
       -uy(ny-3)+uy(ny-5))&
       +cstty*(uy(ny)-uy(ny-3)&
       -uy(ny-3)+uy(ny-6))
  ty(ny-2)=asty*(uy(ny-1)-uy(ny-2)&
       -uy(ny-2)+uy(ny-3))&
       +bsty*(uy(ny)-uy(ny-2)&
       -uy(ny-2)+uy(ny-4))
  ty(ny-1)=asmy*(uy(ny)-uy(ny-1)&
       -uy(ny-1)+uy(ny-2))
  ty(ny)=asny*uy(ny)+bsny*uy(ny-1)&
       +csny*uy(ny-2)+dsny*uy(ny-3)
  if (iimplicit.ge.1) return
  
     do j=2,ny
        ty(j)=ty(j)-ty(j-1)*ssy(j)
     enddo
     
     ty(ny)=ty(ny)*swy(ny)
  
     do j=ny-1,1,-1
        ty(j)=(ty(j)-sfy(j)*ty(j+1))*swy(j)
     enddo
  return
end subroutine deryy1D_22
!---------------------------------------------------------------------------!
