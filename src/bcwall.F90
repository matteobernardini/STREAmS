subroutine bcwall(ilat)
!
! Apply wall boundary conditions
!
 use mod_streams
 implicit none
!
 integer :: i,k,l,ilat
!
 if (ilat==1) then     ! left side
 elseif (ilat==2) then ! right side
 elseif (ilat==3) then ! lower side
 !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    w_gpu(2,i,1,k) = 0._mykind 
    w_gpu(3,i,1,k) = 0._mykind
    w_gpu(4,i,1,k) = 0._mykind 
    w_gpu(5,i,1,k) = w_gpu(1,i,1,k)*gm*taw*trat
    do l=1,ng
     w_gpu(1,i,1-l,k) =    w_gpu(1,i,1+l,k) ! rho   is even
     w_gpu(2,i,1-l,k) = 2._mykind*w_gpu(2,i,1,k) - w_gpu(2,i,1+l,k) ! rho*u is odd	
     w_gpu(3,i,1-l,k) = 2._mykind*w_gpu(3,i,1,k) - w_gpu(3,i,1+l,k) ! rho*v is odd 
     w_gpu(4,i,1-l,k) = 2._mykind*w_gpu(4,i,1,k) - w_gpu(4,i,1+l,k) ! rho*w is odd
     w_gpu(5,i,1-l,k) =    w_gpu(5,i,1+l,k) ! rho*E is even
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==4) then  ! upper side
 !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    w_gpu(2,i,ny,k) = 0._mykind 
    w_gpu(3,i,ny,k) = 0._mykind
    w_gpu(4,i,ny,k) = 0._mykind 
    w_gpu(5,i,ny,k) = w_gpu(1,i,ny,k)*gm*taw*trat
    do l=1,ng
     w_gpu(1,i,ny+l,k) =    w_gpu(1,i,ny-l,k) ! rho   is even
     w_gpu(2,i,ny+l,k) = 2._mykind*w_gpu(2,i,ny,k) - w_gpu(2,i,ny-l,k) ! rho*u is odd	
     w_gpu(3,i,ny+l,k) = 2._mykind*w_gpu(3,i,ny,k) - w_gpu(3,i,ny-l,k) ! rho*v is odd 
     w_gpu(4,i,ny+l,k) = 2._mykind*w_gpu(4,i,ny,k) - w_gpu(4,i,ny-l,k) ! rho*w is odd
     w_gpu(5,i,ny+l,k) =    w_gpu(5,i,ny-l,k) ! rho*E is even
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==5) then  ! back side
 elseif (ilat==6) then  ! fore side
 endif
! 
end subroutine bcwall
