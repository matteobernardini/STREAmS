subroutine computeresidual
!
! Computation of residual
!
 use mod_streams
 implicit none
!
 integer :: i,j,k
 real(mykind) :: rntot,rtrms_ib
 real(mykind), dimension(ny) :: rtrms_ib_1d_cpu
!
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ny
   rtrms_ib_gpu(j,k) = 0._mykind
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ny
   do i=1,nx
    rtrms_ib_gpu(j,k) = rtrms_ib_gpu(j,k)+(fln_gpu(i,j,k,2)/dtglobal)**2
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 !$cuf kernel do(1) <<<*,*>>>
 do j=1,ny
  rtrms_ib_1d_gpu(j) = 0._mykind
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 !$cuf kernel do(1) <<<*,*>>>
 do j=1,ny
  do k=1,nz
   rtrms_ib_1d_gpu(j) = rtrms_ib_1d_gpu(j) + rtrms_ib_gpu(j,k)
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 rtrms_ib_1d_cpu = rtrms_ib_1d_gpu 
!
 rtrms_ib = 0._mykind
 do j=1,ny
  rtrms_ib = rtrms_ib + rtrms_ib_1d_cpu(j)
 enddo
!
 call mpi_allreduce(rtrms_ib,rtrms,1,mpi_prec,mpi_sum,mp_cart,iermpi)
!
 rntot = real(nx,mykind)*real(ny,mykind)*real(nz,mykind)*real(nproc,mykind)
 rtrms = sqrt(rtrms/rntot)
 if (rtrms/=rtrms) then
  if (masterproc) write(*,*) 'BOOM!!!'
  call mpi_barrier(mp_cart,iermpi)
  call mpi_abort(mp_cart,99,iermpi)
 endif
!
end subroutine computeresidual
