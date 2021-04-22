subroutine write_wallpressure
!
! Write wall pressure field and xy slice
!
 use mod_streams
 implicit none
!
 integer :: i,j,k
 real(mykind) :: rho,rhou,rhov,rhow,rhoe,ri,uu,vv,ww,qq,pp
 integer, parameter :: delta_cyc_p=50, delta_cyc_s=50
 integer, parameter :: delta_restart_file=5000
 character(7) :: chicyc
 logical :: async=.true.
 integer, parameter :: i_skip_p=1, k_skip_p=1
 integer, parameter :: i_skip_s=1, j_skip_s=1
!
 if (mod(icyc,delta_cyc_p)==0) then
!
  write(chicyc,"(I7.7)") icyc
!
  !$cuf kernel do(2) <<<*,*>>> 
  do k=1,nz,k_skip_p
   do i=1,nx,i_skip_p
    rho  = w_gpu(i,1,k,1)
    rhou = w_gpu(i,1,k,2)
    rhov = w_gpu(i,1,k,3)
    rhow = w_gpu(i,1,k,4)
    rhoe = w_gpu(i,1,k,5)
    ri   = 1._mykind/rho
    uu   = rhou*ri
    vv   = rhov*ri
    ww   = rhow*ri
    qq   = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
    pp   = gm1*(rhoe-rho*qq)
    wallpfield_gpu(i,k) = pp
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
!
  if(async) then
   if(icyc == ncyc0 + delta_cyc_p) then
    open(122,file='WALLP/wallpressure_'//chx//'_'//chz//'_'//chicyc//'.bin',form='unformatted',asynchronous="yes")
   else
    wait(122)
    if(mod(icyc-ncyc0,delta_restart_file) == 0) then
     close(122)
     open(122,file='WALLP/wallpressure_'//chx//'_'//chz//'_'//chicyc//'.bin',form='unformatted',asynchronous="yes")
    endif
   endif
  endif
!
  wallpfield = wallpfield_gpu
!
  if(async) then
   write(122, asynchronous="no") icyc,telaps,nx,nz !,i_skip_p,k_skip_p
   write(122, asynchronous="yes") wallpfield
  else
   open(122,file='WALLP/wallpressure_'//chx//'_'//chz//'.bin',form='unformatted',position="append")
   write(122) icyc,telaps,nx,nz !,i_skip_p,k_skip_p
   write(122) wallpfield
   close(122)
  endif
!
  if(async) then
   if(icyc == ncyc0+ncyc) then
    wait(122)
    close(122)
   endif
  endif
!
 endif
!
 if (mod(icyc,delta_cyc_s)==0) then
  if (ncoords(3)==0) then
!
   write(chicyc,"(I7.7)") icyc
!
   !$cuf kernel do(2) <<<*,*>>> 
   do j=1,ny,j_skip_s
    do i=1,nx,i_skip_s
     rho  = w_gpu(i,j,1,1)
     rhou = w_gpu(i,j,1,2)
     rhov = w_gpu(i,j,1,3)
     rhow = w_gpu(i,j,1,4)
     rhoe = w_gpu(i,j,1,5)
     ri   = 1._mykind/rho
     uu   = rhou*ri
     vv   = rhov*ri
     ww   = rhow*ri
     qq   = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
     pp   = gm1*(rhoe-rho*qq)
     slicexy_gpu(1,i,j) = rho
     slicexy_gpu(2,i,j) =  uu
     slicexy_gpu(3,i,j) =  vv
     slicexy_gpu(4,i,j) =  ww
     slicexy_gpu(5,i,j) =  pp
    enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
! 
   if(async) then
    if(icyc == ncyc0 + delta_cyc_s) then
     open(133,file='SLICEXY/slicexy_'//chx//'_'//chicyc//'.bin',form='unformatted',asynchronous="yes")
    else
     wait(133)
     if(mod(icyc-ncyc0,delta_restart_file) == 0) then
      close(133)
      open(133,file='SLICEXY/slicexy_'//chx//'_'//chicyc//'.bin',form='unformatted',asynchronous="yes")
     endif
    endif
   endif
! 
   slicexy = slicexy_gpu
!
   if(async) then
     write(133, asynchronous="no") icyc,telaps,nx,ny !,i_skip_s,j_skip_s
     write(133, asynchronous="yes") slicexy
   else
    open(133,file='SLICEXY/slicexy_'//chx//'.bin',form='unformatted',position="append")
     write(133) icyc,telaps,nx,ny !,i_skip_s,j_skip_s
     write(133) slicexy
    close(133)
   endif
!
   if(async) then
    if(icyc == ncyc0+ncyc) then
     wait(133)
     close(133)
    endif
   endif
  endif
 endif
!
end subroutine write_wallpressure
