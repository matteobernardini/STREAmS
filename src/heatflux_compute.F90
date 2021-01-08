subroutine heatflux_compute
!
! Evaluation of the viscous fluxes
!
 use mod_streams
 implicit none
!
 integer :: i,j,k
 real(mykind) :: dy
!
 real(mykind) :: ttt,st,et
!
 heatflux = 0._mykind
 aeroheat = 0._mykind
!
 !$cuf kernel do(2) <<<*,*>>>
 do j=1,ny
  do i=1,nx
   do k=1,nz
    dy = yn_gpu(j+1)-yn_gpu(j)
    heatflux = heatflux+fhat_gpu(i,j,k,1)*dy
    aeroheat = aeroheat+fhat_gpu(i,j,k,2)*dy
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 call mpi_allreduce(MPI_IN_PLACE,heatflux,1,mpi_prec,mpi_sum,mp_cart,iermpi)
 call mpi_allreduce(MPI_IN_PLACE,aeroheat,1,mpi_prec,mpi_sum,mp_cart,iermpi)
 heatflux = heatflux/nxmax/nzmax/(yn(ny+1)-yn(1))
 aeroheat = aeroheat/nxmax/nzmax/(yn(ny+1)-yn(1))
!
end subroutine heatflux_compute
