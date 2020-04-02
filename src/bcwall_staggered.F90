subroutine bcwall_staggered(ilat)
!
! Apply wall boundary conditions (staggered version)
!
 use mod_streams
 implicit none
!
 integer :: i,k,l,ilat
 real(mykind) :: rho,uu,vv,ww,qq,pp,tt,rhoe
!
 if (ilat==1) then     ! left side
 elseif (ilat==2) then ! right side
 elseif (ilat==3) then ! lower side
 !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    do l=1,ng
     rho  = w_gpu(1,i,l,k)
     uu   = w_gpu(2,i,l,k)/w_gpu(1,i,l,k)
     vv   = w_gpu(3,i,l,k)/w_gpu(1,i,l,k)
     ww   = w_gpu(4,i,l,k)/w_gpu(1,i,l,k)
     rhoe = w_gpu(5,i,l,k)
     qq   = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
     pp   = gm1*(rhoe-rho*qq)
     tt   = pp/rho
     tt   = 2._mykind*t0-tt
     rho  = pp/tt
     w_gpu(1,i,1-l,k) =  rho
     w_gpu(2,i,1-l,k) = -rho*uu
     w_gpu(3,i,1-l,k) = -rho*vv
     w_gpu(4,i,1-l,k) = -rho*ww
     w_gpu(5,i,1-l,k) =  pp*gm+qq*rho
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==4) then  ! upper side
 !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    do l=1,ng
     rho  = w_gpu(1,i,ny+1-l,k)
     uu   = w_gpu(2,i,ny+1-l,k)/w_gpu(1,i,ny+1-l,k)
     vv   = w_gpu(3,i,ny+1-l,k)/w_gpu(1,i,ny+1-l,k)
     ww   = w_gpu(4,i,ny+1-l,k)/w_gpu(1,i,ny+1-l,k)
     rhoe = w_gpu(5,i,ny+1-l,k)
     qq   = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
     pp   = gm1*(rhoe-rho*qq)
     tt   = pp/rho
     tt   = 2._mykind*t0-tt
     rho  = pp/tt
     w_gpu(1,i,ny+l,k) =  rho
     w_gpu(2,i,ny+l,k) = -rho*uu
     w_gpu(3,i,ny+l,k) = -rho*vv
     w_gpu(4,i,ny+l,k) = -rho*ww
     w_gpu(5,i,ny+l,k) =  pp*gm+qq*rho
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==5) then  ! back side
 elseif (ilat==6) then  ! fore side
 endif
! 
end subroutine bcwall_staggered
