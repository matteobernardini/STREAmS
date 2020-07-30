subroutine bcextr(ilat)
!
! Apply extrapolation to fill ghost nodes
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,l,m,ilat
!
 if (ilat==1) then     ! left side
 !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do l=1,ng
     do m=1,nv
      w_gpu(1-l,j,k,m) = w_gpu(1,j,k,m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==2) then ! right side
 !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do l=1,ng
     do m=1,nv
      w_gpu(nx+l,j,k,m) = w_gpu(nx,j,k,m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==3) then  ! lower side
 !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    do l=1,ng
     do m=1,nv
      w_gpu(i,1-l,k,m) = w_gpu(i,1,k,m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==4) then  ! upper side
 !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    do l=1,ng
     do m=1,nv
      w_gpu(i,ny+l,k,m) = w_gpu(i,ny,k,m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==5) then  ! back side
 !$cuf kernel do(2) <<<*,*>>>
  do j=1,ny
   do i=1,nx
    do l=1,ng
     do m=1,nv
      w_gpu(i,j,1-l,m) = w_gpu(i,j,1,m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==6) then  ! fore side
 !$cuf kernel do(2) <<<*,*>>>
  do j=1,ny
   do i=1,nx
    do l=1,ng
     do m=1,nv
      w_gpu(i,j,nz+l,m) = w_gpu(i,j,nz,m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 endif
! 
end subroutine bcextr
