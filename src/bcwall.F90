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
    w_gpu(i,1,k,2) = 0._mykind 
    w_gpu(i,1,k,3) = 0._mykind
    w_gpu(i,1,k,4) = 0._mykind 
    w_gpu(i,1,k,5) = w_gpu(i,1,k,1)*gm*taw*trat
    do l=1,ng
     w_gpu(i,1-l,k,1) =    w_gpu(i,1+l,k,1) ! rho   is even
     w_gpu(i,1-l,k,2) = 2._mykind*w_gpu(i,1,k,2) - w_gpu(i,1+l,k,2) ! rho*u is odd	
     w_gpu(i,1-l,k,3) = 2._mykind*w_gpu(i,1,k,3) - w_gpu(i,1+l,k,3) ! rho*v is odd 
     w_gpu(i,1-l,k,4) = 2._mykind*w_gpu(i,1,k,4) - w_gpu(i,1+l,k,4) ! rho*w is odd
     w_gpu(i,1-l,k,5) =    w_gpu(i,1+l,k,5) ! rho*E is even
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==4) then  ! upper side
 !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    w_gpu(i,ny,k,2) = 0._mykind 
    w_gpu(i,ny,k,3) = 0._mykind
    w_gpu(i,ny,k,4) = 0._mykind 
    w_gpu(i,ny,k,5) = w_gpu(i,ny,k,1)*gm*taw*trat
    do l=1,ng
     w_gpu(i,ny+l,k,1) =    w_gpu(i,ny-l,k,1) ! rho   is even
     w_gpu(i,ny+l,k,2) = 2._mykind*w_gpu(i,ny,k,2) - w_gpu(i,ny-l,k,2) ! rho*u is odd	
     w_gpu(i,ny+l,k,3) = 2._mykind*w_gpu(i,ny,k,3) - w_gpu(i,ny-l,k,3) ! rho*v is odd 
     w_gpu(i,ny+l,k,4) = 2._mykind*w_gpu(i,ny,k,4) - w_gpu(i,ny-l,k,4) ! rho*w is odd
     w_gpu(i,ny+l,k,5) =    w_gpu(i,ny-l,k,5) ! rho*E is even
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==5) then  ! back side
 elseif (ilat==6) then  ! fore side
 endif
! 
end subroutine bcwall
