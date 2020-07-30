subroutine bcshk(ilat)
!
! Apply boundary conditions for generation of the oblique shock
!
 use mod_streams
 implicit none
!
 integer :: ilat,i,k,l,m
 real(mykind) :: dwinf,ftanh,xx,yy
!
 if (ilat==4) then     ! left side
 !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    do l=0,ng
     yy  = y_gpu(ny+l)
     xsh = ximp-yy*tan(thetas)
     xx  = x_gpu(i)-xsh
!    xxtan = xx/tanhsfac
!    if (abs(xxtan)>3) then
     if (abs(xx)>.2_mykind) then
      do m=1,nv
       w_gpu(i,ny+l,k,m) = w_gpu(i,ny,k,m)
      enddo
     else
!     ftanh = 0.5_mykind*(1._mykind+tanh(xxtan))
      ftanh = 0.5_mykind*(1._mykind+tanh(xx/tanhsfac))
      do m=1,nv
       dwinf = winf1_gpu(m)-winf_gpu(m)
       w_gpu(i,ny+l,k,m) = winf_gpu(m)+dwinf*ftanh
      enddo
     endif
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 endif
!
end subroutine bcshk
