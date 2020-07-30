subroutine sensor
!
!Evaluation of Ducros shock sensor
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,l
 real(mykind) :: ccl,div,div2,epsi,epsi2,omod,omod2
 real(mykind) :: omx,omy,omz,ri
 real(mykind) :: ux,uy,uz
 real(mykind) :: vx,vy,vz
 real(mykind) :: wx,wy,wz
 real(mykind) :: ducloc
!
 !$cuf kernel do(3) <<<*,*>>> 
 do k=1,nz
  do j=1,ny
   do i=1,nx
!
    ux = 0._mykind
    vx = 0._mykind
    wx = 0._mykind
    uy = 0._mykind
    vy = 0._mykind
    wy = 0._mykind
    uz = 0._mykind
    vz = 0._mykind
    wz = 0._mykind
    do l=1,ivis/2
     ccl = coeff_deriv1_gpu(l)
     ux = ux+ccl*(w_gpu(i+l,j,k,2)/w_gpu(i+l,j,k,1)-w_gpu(i-l,j,k,2)/w_gpu(i-l,j,k,1))
     vx = vx+ccl*(w_gpu(i+l,j,k,3)/w_gpu(i+l,j,k,1)-w_gpu(i-l,j,k,3)/w_gpu(i-l,j,k,1))
     wx = wx+ccl*(w_gpu(i+l,j,k,4)/w_gpu(i+l,j,k,1)-w_gpu(i-l,j,k,4)/w_gpu(i-l,j,k,1))
     uy = uy+ccl*(w_gpu(i,j+l,k,2)/w_gpu(i,j+l,k,1)-w_gpu(i,j-l,k,2)/w_gpu(i,j-l,k,1))
     vy = vy+ccl*(w_gpu(i,j+l,k,3)/w_gpu(i,j+l,k,1)-w_gpu(i,j-l,k,3)/w_gpu(i,j-l,k,1))
     wy = wy+ccl*(w_gpu(i,j+l,k,4)/w_gpu(i,j+l,k,1)-w_gpu(i,j-l,k,4)/w_gpu(i,j-l,k,1))
     uz = uz+ccl*(w_gpu(i,j,k+l,2)/w_gpu(i,j,k+l,1)-w_gpu(i,j,k-l,2)/w_gpu(i,j,k-l,1))
     vz = vz+ccl*(w_gpu(i,j,k+l,3)/w_gpu(i,j,k+l,1)-w_gpu(i,j,k-l,3)/w_gpu(i,j,k-l,1))
     wz = wz+ccl*(w_gpu(i,j,k+l,4)/w_gpu(i,j,k+l,1)-w_gpu(i,j,k-l,4)/w_gpu(i,j,k-l,1))
    enddo
    ux = ux*dcsidx_gpu(i)
    vx = vx*dcsidx_gpu(i)
    wx = wx*dcsidx_gpu(i)
    uy = uy*detady_gpu(j)
    vy = vy*detady_gpu(j)
    wy = wy*detady_gpu(j)
    uz = uz*dzitdz_gpu(k)
    vz = vz*dzitdz_gpu(k)
    wz = wz*dzitdz_gpu(k)
!
    epsi2 = u0**2
    epsi  = u0
!
    omz = vx-uy
    omx = wy-vz
    omy = uz-wx
    omod2 = omx*omx+omy*omy+omz*omz
    omod = sqrt(omod2)
    div = ux+vy+wz
    div2 = div*div
!
    ducloc = max(-div/sqrt(omod2+div2+epsi2),0._mykind)
    ducros_gpu(i,j,k) = .false.
    if (ducloc.gt.tresduc) ducros_gpu(i,j,k) = .true.
!
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
end subroutine sensor
