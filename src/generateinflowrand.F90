subroutine generateinflowrand
!
! Generation of random number at the inflow plane for digital filtering
!
 use mod_streams
 implicit none
!
 integer :: m,j,k
!
! Create a new random field (local from 1 to nz)
!
 do k = 1,nz
  do j = 1-nfmax,ny+nfmax
   do m = 1,3
    call gasdev_s(rf(m,j,k))
   enddo
  enddo
 enddo
!
! MPI rf
!
 if (nblocks(3)==1) then
!
  do k = 1,nfmax
   do j = 1-nfmax,ny+nfmax
    do m = 1,3
     rf(m,j,1-k)  = rf(m,j,nz+1-k)
     rf(m,j,nz+k) = rf(m,j,k)
    enddo
   enddo
  enddo
!
 else
!
  call swapz_extended_rf(rf)
!
 endif
!
end subroutine generateinflowrand
