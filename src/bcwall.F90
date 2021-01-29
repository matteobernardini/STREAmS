subroutine bcwall(ilat)
!
! Apply wall boundary conditions
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
    w_gpu(i,1,k,2) = 0._mykind 
    w_gpu(i,1,k,3) = 0._mykind
    w_gpu(i,1,k,4) = 0._mykind 
    w_gpu(i,1,k,5) = w_gpu(i,1,k,1)*gm*twall
    do l=1,ng
     rho  = w_gpu(i,1+l,k,1)
     uu   = w_gpu(i,1+l,k,2)/w_gpu(i,1+l,k,1)
     vv   = w_gpu(i,1+l,k,3)/w_gpu(i,1+l,k,1)
     ww   = w_gpu(i,1+l,k,4)/w_gpu(i,1+l,k,1)
     rhoe = w_gpu(i,1+l,k,5)
     qq   = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
     pp   = gm1*(rhoe-rho*qq)
     tt   = pp/rho
     tt   = 2._mykind*twall-tt
     rho  = pp/tt
     w_gpu(i,1-l,k,1) =  rho
     w_gpu(i,1-l,k,2) = -rho*uu
     w_gpu(i,1-l,k,3) = -rho*vv
     w_gpu(i,1-l,k,4) = -rho*ww
     w_gpu(i,1-l,k,5) = pp*gm+qq*rho
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
    w_gpu(i,ny,k,5) = w_gpu(i,ny,k,1)*gm*twall
    do l=1,ng
     rho  = w_gpu(i,ny-l,k,1)
     uu   = w_gpu(i,ny-l,k,2)/w_gpu(i,ny-l,k,1)
     vv   = w_gpu(i,ny-l,k,3)/w_gpu(i,ny-l,k,1)
     ww   = w_gpu(i,ny-l,k,4)/w_gpu(i,ny-l,k,1)
     rhoe = w_gpu(i,ny-l,k,5)
     qq   = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
     pp   = gm1*(rhoe-rho*qq)
     tt   = pp/rho
     tt   = 2._mykind*twall-tt
     rho  = pp/tt
     w_gpu(i,ny+l,k,1) =  rho
     w_gpu(i,ny+l,k,2) = -rho*uu
     w_gpu(i,ny+l,k,3) = -rho*vv
     w_gpu(i,ny+l,k,4) = -rho*ww
     w_gpu(i,ny+l,k,5) = pp*gm+qq*rho
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==5) then  ! back side
 elseif (ilat==6) then  ! fore side
 endif
! 
end subroutine bcwall
