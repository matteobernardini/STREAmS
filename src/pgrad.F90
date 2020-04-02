subroutine pgrad
!
! Evaluate pressure gradient for channel flow
 use mod_streams
!
 implicit none
!
 integer :: i,j,k,l
 real(mykind), dimension(5) :: bulk5
 real(mykind), dimension(5) :: bulk5g
 real(mykind), dimension(5) :: bulk5g_gpu
 real(mykind) :: dy,uu
 real(mykind) :: bulk_1,bulk_2,bulk_3,bulk_4,bulk_5
#ifdef USE_CUDA
 attributes(device) :: bulk5g_gpu
#endif
!
 bulk_1 = 0._mykind
 bulk_2 = 0._mykind
 bulk_3 = 0._mykind
 bulk_4 = 0._mykind
 bulk_5 = 0._mykind
!
 !$cuf kernel do(3) <<<*,*>>> 
 do k=1,nz
  do j=1,ny
   do i=1,nx
    dy = yn_gpu(j+1)-yn_gpu(j)
    bulk_1 = bulk_1 + fln_gpu(1,i,j,k)*dy
    bulk_2 = bulk_2 + fln_gpu(2,i,j,k)*dy
    bulk_3 = bulk_3 + w_gpu(1,i,j,k)*dy
    bulk_4 = bulk_4 + w_gpu(2,i,j,k)*dy
    bulk_5 = bulk_5 + w_gpu(1,i,j,k)*temperature_gpu(i,j,k)*dy
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 bulk5(1) = bulk_1
 bulk5(2) = bulk_2
 bulk5(3) = bulk_3
 bulk5(4) = bulk_4
 bulk5(5) = bulk_5
!
 call mpi_allreduce(bulk5,bulk5g,5,mpi_prec,mpi_sum,mp_cart,iermpi)
!
 bulk5g = bulk5g/nxmax/nzmax/(yn_gpu(ny+1)-yn_gpu(1))
!
 dpdx = bulk5g(2)
 rhobulk = bulk5g(3)
 ubulk = bulk5g(4)/rhobulk
 tbulk = bulk5g(5)/rhobulk
 bulk5g_gpu = bulk5g
!
! Add forcing terms in momentum and energy equation
 !$cuf kernel do(3) <<<*,*>>> 
 do k=1,nz
  do j=1,ny
   do i=1,nx
    uu = w_gpu(2,i,j,k)/w_gpu(1,i,j,k)
    fln_gpu(1,i,j,k) = fln_gpu(1,i,j,k) -    bulk5g_gpu(1)
    fln_gpu(2,i,j,k) = fln_gpu(2,i,j,k) -    bulk5g_gpu(2)
    fln_gpu(5,i,j,k) = fln_gpu(5,i,j,k) - uu*bulk5g_gpu(2)
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 end
