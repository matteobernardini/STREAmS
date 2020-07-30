subroutine bcfree(ilat)
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
      w_gpu(1-l,j,k,m) = winf_gpu(m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 endif
! 
end subroutine bcfree
