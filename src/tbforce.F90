subroutine tbforce
!
!Adjust bulk temperature
 use mod_streams
!
 implicit none
!
 integer :: i,j,k,l
 real(mykind) :: rho,rhou,rhov,rhow,rhoe
 real(mykind) :: ri,uu,vv,ww,qq,pp,tt
 real(mykind) :: tbtarget,trelax,ttnew
 real(mykind) :: dy
 real(mykind) :: bulkt,bulktg
!
!trelax = 1._mykind
 trelax = min(1._mykind,telaps/10._mykind)
 bulkt  = 0._mykind
! 
 !$cuf kernel do(3) <<<*,*>>> 
 do k=1,nz
  do j=1,ny
   do i=1,nx
    rho  = w_gpu(i,j,k,1)
    rhou = w_gpu(i,j,k,2)
    rhov = w_gpu(i,j,k,3)
    rhow = w_gpu(i,j,k,4)
    rhoe = w_gpu(i,j,k,5)
    ri   = 1._mykind/rho
    uu   = rhou*ri
    vv   = rhov*ri
    ww   = rhow*ri
    qq   = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
    pp   = gm1*(rhoe-rho*qq)
    tt   = pp*ri
    dy = yn_gpu(j+1)-yn_gpu(j)
    bulkt = bulkt + rhou*tt*dy
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
! 
 call mpi_allreduce(bulkt,bulktg,1,mpi_prec,mpi_sum,mp_cart,iermpi)
! 
 bulktg   = bulktg/nxmax/nzmax/(yn(ny+1)-yn(1))
 bulktg   = bulktg/rhobulk/ubulk
 tbtarget = trelax*target_tbulk+(1._mykind-trelax)*bulktg
!
 bulkcooling = 0._mykind
!
 !$cuf kernel do(3) <<<*,*>>> 
 do k=1,nz
  do j=1,ny
   do i=1,nx
    rho   = w_gpu(i,j,k,1)
    rhou  = w_gpu(i,j,k,2)
    rhov  = w_gpu(i,j,k,3)
    rhow  = w_gpu(i,j,k,4)
    rhoe  = w_gpu(i,j,k,5)
    ri    = 1._mykind/rho
    uu    = rhou*ri
    vv    = rhov*ri
    ww    = rhow*ri
    qq    = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
    pp    = gm1*(rhoe-rho*qq)
    tt    = pp*ri
    ttnew = tt+tbtarget-bulktg 
    pp    = rho*ttnew
    w_gpu(i,j,k,5) = pp*gm+rho*qq
    dy = yn_gpu(j+1)-yn_gpu(j)
    bulkcooling = bulkcooling+(w_gpu(i,j,k,5)-rhoe)*dy
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 if (mod(icyc-ncyc0,nprint)==0) then
  call mpi_allreduce(MPI_IN_PLACE,bulkcooling,1,mpi_prec,mpi_sum,mp_cart,iermpi)
  bulkcooling = bulkcooling/nxmax/nzmax/(yn(ny+1)-yn(1))/alpdt
 endif
!
end subroutine tbforce
