subroutine rk
!
! 3-rd order RK solver to advance in time
!
 use mod_streams
 use mod_euler
 implicit none
!
 integer :: i,j,k,m,istep,ilat
 real(mykind) :: alp,dt,gam,gamdt,rho,rhodt
 real(mykind) :: elapsed,startTiming,endTiming
!
 real(mykind) :: tt,st,et
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
!st = mpi_wtime()
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
!et = mpi_wtime()
!tt = et-st
!if (masterproc) write(error_unit,*) 'RK-I time =', tt
!
#ifdef CUDA_ASYNC
  call bc(0)
  call bcswap_prepare()
  call prims_int()
  call euler_i(0+iorder/2,nx-iorder/2)
  call bcswap()
  call prims_ghost()
!
! Evaluation of Eulerian fluxes
!
  call euler_i(0,iorder/2-1)
  if (ibc(2)==4.or.ibc(2)==8) then
   call euler_i(nx+1-iorder/2,nx-1)
  else
   call euler_i(nx+1-iorder/2,nx)
  endif
#else
  call updateghost()
  call prims()
  if (ibc(2)==4.or.ibc(2)==8) then
   call euler_i(0,nx-1)
  else
   call euler_i(0,nx)
  endif
#endif
!
#ifdef CUDA_ASYNC
  call visflx()
! call visflx_stag()
  call bcswapdiv_prepare()
  call euler_j()
  call bcswapdiv()
  if (ndim==3) call euler_k()
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
  call euler_j()
  if (ndim==3) call euler_k()
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
