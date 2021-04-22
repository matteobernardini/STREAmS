subroutine rk
!
! 3-rd order RK solver to advance in time
!
 use mod_streams
 use mod_euler
 implicit none
!
 integer :: i1,i2,j1,j2,k1,k2,lmax,last_call
 logical :: call2,call3
!
 integer :: i,j,k,m,istep,ilat
 real(mykind) :: alp,dt,gam,gamdt,rho,rhodt
 real(mykind) :: elapsed,startTiming,endTiming
!
 lmax = iorder/2
 i1 = 1  ; if (ibc(1)==9.or.ibc(1)==10) i1 = 2    ! Node 1  computed with characteristic decomposition
 i2 = nx ; if (ibc(2)==4.or.ibc(2)==8)  i2 = nx-1 ! Node nx computed with characteristic decomposition
 j1 = 1  ; if (ibc(3)==4.or.ibc(3)==8)  j1 = 2    ! Node 1  computed with characteristic decomposition
 j2 = ny ; if (ibc(4)==4.or.ibc(4)==7)  j2 = ny-1 ! Node ny computed with characteristic decomposition
 k1 = 1  ; if (ibc(5)==4.or.ibc(5)==8)  k1 = 2    ! Node 1  computed with characteristic decomposition
 k2 = nz ; if (ibc(6)==4.or.ibc(6)==7)  k2 = nz-1 ! Node nz computed with characteristic decomposition
 last_call = 3
 call2 = .false.
 call3 = .false.
 if (i1<=lmax) call2 = .true.
 if (i2>=(nx-lmax+1)) call3 = .true.
 if (.not.call3) then
  if (call2) then
   last_call = 2
  else
   last_call = 1
  endif
 endif
!
 if (cfl>0._mykind) then
  dt = dtmin*cfl 
 else
  dt = dtmin
 endif
 dtglobal = dt
!
! Loop on the 'nstage' stages of RK solver
!
 do istep=1,3 ! Alan Wray 3rd order RK method
!
  rho = rhovec(istep) ! Coefficient for nonlinear terms
  gam = gamvec(istep) ! Coefficient for nonlinear terms
  alp = gam+rho       ! Coefficient for linear terms
  rhodt = rho*dt
  gamdt = gam*dt
  alpdt = alp*dt
!
 !$cuf kernel do(3) <<<*,*>>> 
  do k=1,nz
   do j=1,ny
    do i=1,nx
     do m=1,nv
      fln_gpu(i,j,k,m) = -rhodt*fl_gpu(i,j,k,m)
      fl_gpu(i,j,k,m) = 0._mykind
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 if (xrecyc>0._mykind) call recyc
!
#ifdef CUDA_ASYNC
  call bc(0)
  call bcswap_prepare()
  call prims_int()
  call euler_i(1+lmax,nx-lmax,1,last_call,i1,i2)
  call bcswap()
  call prims_ghost()
!
  if (call2) call euler_i(i1,lmax,2,last_call,i1,i2)
  if (call3) call euler_i(nx-lmax+1,i2,3,last_call,i1,i2)
!
#else
  call updateghost()
  call prims()
  call euler_i(i1,i2,0,0,i1,i2)
#endif
!
#ifdef CUDA_ASYNC
  call visflx()
! call visflx_stag()
  call bcswapdiv_prepare()
  call euler_j(j1,j2)
  call bcswapdiv()
  if (ndim==3) call euler_k(k1,k2)
!
  if (istep==3.and.tresduc>0._mykind.and.tresduc<1._mykind) then
   call sensor()
   call bcswapduc_prepare()
   call visflx_div() ! No Cuda Sync here
   call bcswapduc()
  else
   call visflx_div()
   !@cuf iercuda=cudaDeviceSynchronize()
  endif
!
#else
  call euler_j(j1,j2)
  if (ndim==3) call euler_k(k1,k2)
  call visflx()
! call visflx_stag()
  call bcswapdiv_prepare()
  call bcswapdiv()
  if (istep==3.and.tresduc>0._mykind.and.tresduc<1._mykind) then
   call sensor()
   call bcswapduc_prepare()
   call bcswapduc()
  endif
  call visflx_div() ! No Cuda Sync here
  !@cuf iercuda=cudaDeviceSynchronize()
#endif
!
! Call to non-reflecting b.c. (to update f_x, g_y and h_z on the boundaries)
  call bc(1)
!
 !$cuf kernel do(3) <<<*,*>>> 
  do k=1,nz
   do j=1,ny
    do i=1,nx
     do m=1,nv
      fln_gpu(i,j,k,m) = fln_gpu(i,j,k,m)-gamdt*fl_gpu(i,j,k,m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
  if (iflow==0) then
   call pgrad()
   dpdx = -dpdx/alpdt
  endif
!
! Updating solution in inner nodes
!
 !$cuf kernel do(3) <<<*,*>>> 
  do k=1,nz
   do j=1,ny
    do i=1,nx
     do m=1,nv
      w_gpu(i,j,k,m) = w_gpu(i,j,k,m)+fln_gpu(i,j,k,m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
  if (iflow==0) then
   if (trat>0._mykind) call tbforce() ! For channel flow adjust bulk temperature
  endif
!
#ifdef USE_CUDA
 iercuda = cudaGetLastError()
 if (iercuda /= cudaSuccess) then
  call fail_input("CUDA ERROR! Try to reduce the number of Euler threads in cuda_definitions.h: "//cudaGetErrorString(iercuda))
 endif
#endif
!
  alpdtold  = alpdt 
  dfupdated = .false.
!
 enddo
!
 telaps = telaps+dt
!
end subroutine rk
