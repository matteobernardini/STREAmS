subroutine printres
!
! Printing residual
!
 use mod_streams
 implicit none
!
 integer :: i,j,k
 real(mykind) :: rhomin,rhomax
 real(mykind) :: sumq
 real(mykind) :: xmin,ymin,zmin
!
 rhomin =  100._mykind
 rhomax = -100._mykind
 !$cuf kernel do(3) <<<*,*>>>
 do k=1,nz
  do j=1,ny
   do i=1,nx
    rhomin = min(rhomin,w_gpu(i,j,k,1)*temperature_gpu(i,j,k))
    rhomax = max(rhomax,w_gpu(i,j,k,1)*temperature_gpu(i,j,k))
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 call mpi_allreduce(MPI_IN_PLACE,rhomin,1,mpi_prec,mpi_min,mp_cart,iermpi)
 call mpi_allreduce(MPI_IN_PLACE,rhomax,1,mpi_prec,mpi_max,mp_cart,iermpi)
!
 if (masterproc) then
  if (iflow==0) then
   sumq = heatflux+aeroheat+bulkcooling
   write(* ,100) icyc,telaps,rtrms,dpdx,rhobulk,ubulk,tbulk,heatflux,aeroheat,bulkcooling,sumq
   write(20,100) icyc,telaps,rtrms,dpdx,rhobulk,ubulk,tbulk,heatflux,aeroheat,bulkcooling,sumq
  else
   write(* ,100) icyc,telaps,rtrms,rhomin,rhomax
   write(20,100) icyc,telaps,rtrms,rhomin,rhomax
  endif
 endif
!
 100 format(1I10,40ES20.10)
!
end subroutine printres
