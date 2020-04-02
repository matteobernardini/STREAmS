subroutine bcswap_prepare
!
! Preparing arrays for data exchange
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,m
!
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ny
   do i=1,ng
    do m=1,nv
     wbuf1s_gpu(m,i,j,k) = w_gpu(m,i,j,k)
     wbuf2s_gpu(m,i,j,k) = w_gpu(m,nx-ng+i,j,k)
    enddo
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ng
   do i=1,nx
    do m=1,nv
     wbuf3s_gpu(m,i,j,k) = w_gpu(m,i,ny-ng+j,k)
     wbuf4s_gpu(m,i,j,k) = w_gpu(m,i,j,k)
    enddo
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 if (ndim==3) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,ng
   do j=1,ny
    do i=1,nx
     do m=1,nv
      wbuf5s_gpu(m,i,j,k) = w_gpu(m,i,j,k)
      wbuf6s_gpu(m,i,j,k) = w_gpu(m,i,j,nz-ng+k)
     enddo
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
!
end subroutine bcswap_prepare

subroutine bcswap
!
! Apply virtual interface boundary conditions
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,m
 integer :: indx,indy,indz
!
 indx = nv*ng*ny*nz
 indy = nv*nx*ng*nz
 indz = nv*nx*ny*ng
!
#if defined(USE_CUDA) && defined(NOCUDAAWAREMPI)
 iermpi = cudaMemcpyAsync(wbuf1s, wbuf1s_gpu, indx, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(wbuf2s, wbuf2s_gpu, indx, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(wbuf3s, wbuf3s_gpu, indy, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(wbuf4s, wbuf4s_gpu, indy, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(wbuf5s, wbuf5s_gpu, indz, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(wbuf6s, wbuf6s_gpu, indz, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaStreamSynchronize(stream2)
 !wbuf1s = wbuf1s_gpu ; wbuf2s = wbuf2s_gpu ; wbuf3s = wbuf3s_gpu 
 !wbuf4s = wbuf4s_gpu ; wbuf5s = wbuf5s_gpu ; wbuf6s = wbuf6s_gpu 
 call mpi_sendrecv(wbuf1s,indx,mpi_prec,ileftx ,1,wbuf2r,indx,mpi_prec,irightx,1,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(wbuf2s,indx,mpi_prec,irightx,2,wbuf1r,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(wbuf3s,indy,mpi_prec,ilefty ,3,wbuf4r,indy,mpi_prec,irighty,3,mp_carty,istatus,iermpi)       
 call mpi_sendrecv(wbuf4s,indy,mpi_prec,irighty,4,wbuf3r,indy,mpi_prec,ilefty ,4,mp_carty,istatus,iermpi)       
 if (ndim==3) then
   call mpi_sendrecv(wbuf5s,indz,mpi_prec,ileftz ,5,wbuf6r,indz,mpi_prec,irightz,5,mp_cartz,istatus,iermpi)       
   call mpi_sendrecv(wbuf6s,indz,mpi_prec,irightz,6,wbuf5r,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,iermpi)       
 endif
 iermpi = cudaMemcpyAsync(wbuf1r_gpu, wbuf1r, indx, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(wbuf2r_gpu, wbuf2r, indx, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(wbuf3r_gpu, wbuf3r, indy, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(wbuf4r_gpu, wbuf4r, indy, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(wbuf5r_gpu, wbuf5r, indz, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(wbuf6r_gpu, wbuf6r, indz, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaStreamSynchronize(stream2)
 !wbuf1r_gpu = wbuf1r ; wbuf2r_gpu = wbuf2r ; wbuf3r_gpu = wbuf3r
 !wbuf4r_gpu = wbuf4r ; wbuf5r_gpu = wbuf5r ; wbuf6r_gpu = wbuf6r
#else
 call mpi_sendrecv(wbuf1s_gpu,indx,mpi_prec,ileftx ,1,wbuf2r_gpu,indx,mpi_prec,irightx,1,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(wbuf2s_gpu,indx,mpi_prec,irightx,2,wbuf1r_gpu,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(wbuf3s_gpu,indy,mpi_prec,ilefty ,3,wbuf4r_gpu,indy,mpi_prec,irighty,3,mp_carty,istatus,iermpi)       
 call mpi_sendrecv(wbuf4s_gpu,indy,mpi_prec,irighty,4,wbuf3r_gpu,indy,mpi_prec,ilefty ,4,mp_carty,istatus,iermpi)       
 if (ndim==3) then
   call mpi_sendrecv(wbuf5s_gpu,indz,mpi_prec,ileftz ,5,wbuf6r_gpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,iermpi)       
   call mpi_sendrecv(wbuf6s_gpu,indz,mpi_prec,irightz,6,wbuf5r_gpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,iermpi)       
 endif
#endif
!
 if (ileftx/=mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do i=1,ng
     do m=1,nv
      w_gpu(m,i-ng,j,k) = wbuf1r_gpu(m,i,j,k)
     enddo
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (irightx/=mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do i=1,ng
     do m=1,nv
      w_gpu(m,nx+i,j,k) = wbuf2r_gpu(m,i,j,k)
     enddo
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (ilefty/=mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ng
    do i=1,nx
     do m=1,nv
      w_gpu(m,i,j-ng,k) = wbuf3r_gpu(m,i,j,k)
     enddo
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (irighty/=mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ng
    do i=1,nx
     do m=1,nv
      w_gpu(m,i,ny+j,k) = wbuf4r_gpu(m,i,j,k)
     enddo
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (ndim==3) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,ng
   do j=1,ny
    do i=1,nx
     do m=1,nv
      w_gpu(m,i,j,k-ng) = wbuf5r_gpu(m,i,j,k)
      w_gpu(m,i,j,nz+k) = wbuf6r_gpu(m,i,j,k)
     enddo
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
!
end subroutine bcswap
!
subroutine bcswapdiv_prepare
!
! Preparing arrays for data exchange (flow divergence)
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,m
!
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ny
   do i=1,ng
    divbuf1s_gpu(i,j,k) = fhat_gpu(6,i,j,k)
    divbuf2s_gpu(i,j,k) = fhat_gpu(6,nx-ng+i,j,k)
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ng
   do i=1,nx
    divbuf3s_gpu(i,j,k) = fhat_gpu(6,i,ny-ng+j,k)
    divbuf4s_gpu(i,j,k) = fhat_gpu(6,i,j,k)
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 if (ndim==3) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,ng
   do j=1,ny
    do i=1,nx
     divbuf5s_gpu(i,j,k) = fhat_gpu(6,i,j,k)
     divbuf6s_gpu(i,j,k) = fhat_gpu(6,i,j,nz-ng+k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
!
end subroutine bcswapdiv_prepare

subroutine bcswapdiv
!
! Apply virtual interface boundary conditions (flow divergence)
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,m
 integer :: indx,indy,indz
!
 indx = ng*ny*nz
 indy = nx*ng*nz
 indz = nx*ny*ng
!
#if defined(USE_CUDA) && defined(NOCUDAAWAREMPI)
 !divbuf1s = divbuf1s_gpu ; divbuf2s = divbuf2s_gpu ; divbuf3s = divbuf3s_gpu
 !divbuf4s = divbuf4s_gpu ; divbuf5s = divbuf5s_gpu ; divbuf6s = divbuf6s_gpu
 iermpi = cudaMemcpyAsync(divbuf1s, divbuf1s_gpu, indx, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(divbuf2s, divbuf2s_gpu, indx, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(divbuf3s, divbuf3s_gpu, indy, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(divbuf4s, divbuf4s_gpu, indy, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(divbuf5s, divbuf5s_gpu, indz, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(divbuf6s, divbuf6s_gpu, indz, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaStreamSynchronize(stream2)
 call mpi_sendrecv(divbuf1s,indx,mpi_prec,ileftx ,1,divbuf2r,indx,mpi_prec,irightx,1,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(divbuf2s,indx,mpi_prec,irightx,2,divbuf1r,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(divbuf3s,indy,mpi_prec,ilefty ,3,divbuf4r,indy,mpi_prec,irighty,3,mp_carty,istatus,iermpi)       
 call mpi_sendrecv(divbuf4s,indy,mpi_prec,irighty,4,divbuf3r,indy,mpi_prec,ilefty ,4,mp_carty,istatus,iermpi)       
 if (ndim==3) then
  call mpi_sendrecv(divbuf5s,indz,mpi_prec,ileftz ,5,divbuf6r,indz,mpi_prec,irightz,5,mp_cartz,istatus,iermpi)       
  call mpi_sendrecv(divbuf6s,indz,mpi_prec,irightz,6,divbuf5r,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,iermpi)       
 endif
 iermpi = cudaMemcpyAsync(divbuf1r_gpu, divbuf1r, indx, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(divbuf2r_gpu, divbuf2r, indx, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(divbuf3r_gpu, divbuf3r, indy, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(divbuf4r_gpu, divbuf4r, indy, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(divbuf5r_gpu, divbuf5r, indz, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(divbuf6r_gpu, divbuf6r, indz, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaStreamSynchronize(stream2)
 !divbuf1r_gpu = divbuf1r ; divbuf2r_gpu = divbuf2r ; divbuf3r_gpu = divbuf3r
 !divbuf4r_gpu = divbuf4r ; divbuf5r_gpu = divbuf5r ; divbuf6r_gpu = divbuf6r
#else
 call mpi_sendrecv(divbuf1s_gpu,indx,mpi_prec,ileftx ,1,divbuf2r_gpu,indx,mpi_prec,irightx,1,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(divbuf2s_gpu,indx,mpi_prec,irightx,2,divbuf1r_gpu,indx,mpi_prec,ileftx ,2,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(divbuf3s_gpu,indy,mpi_prec,ilefty ,3,divbuf4r_gpu,indy,mpi_prec,irighty,3,mp_carty,istatus,iermpi)       
 call mpi_sendrecv(divbuf4s_gpu,indy,mpi_prec,irighty,4,divbuf3r_gpu,indy,mpi_prec,ilefty ,4,mp_carty,istatus,iermpi)       
 if (ndim==3) then
  call mpi_sendrecv(divbuf5s_gpu,indz,mpi_prec,ileftz ,5,divbuf6r_gpu,indz,mpi_prec,irightz,5,mp_cartz,istatus,iermpi)       
  call mpi_sendrecv(divbuf6s_gpu,indz,mpi_prec,irightz,6,divbuf5r_gpu,indz,mpi_prec,ileftz ,6,mp_cartz,istatus,iermpi)       
 endif
#endif
!
 if (ileftx==mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do i=1,ng
     fhat_gpu(6,1-i,j,k) = 2._mykind*fhat_gpu(6,2-i,j,k)-fhat_gpu(6,3-i,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 else
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do i=1,ng
     fhat_gpu(6,i-ng,j,k) = divbuf1r_gpu(i,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (irightx==mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do i=1,ng
     fhat_gpu(6,nx+i,j,k) = 2._mykind*fhat_gpu(6,nx+i-1,j,k)-fhat_gpu(6,nx+i-2,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 else
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do i=1,ng
     fhat_gpu(6,nx+i,j,k) = divbuf2r_gpu(i,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (ilefty==mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    do j=1,ng
     fhat_gpu(6,i,1-j,k) = 2._mykind*fhat_gpu(6,i,2-j,k)-fhat_gpu(6,i,3-j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 else
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ng
    do i=1,nx
     fhat_gpu(6,i,j-ng,k) = divbuf3r_gpu(i,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (irighty==mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    do j=1,ng
     fhat_gpu(6,i,ny+j,k) = 2._mykind*fhat_gpu(6,i,ny+j-1,k)-fhat_gpu(6,i,ny+j-2,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 else
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ng
    do i=1,nx
     fhat_gpu(6,i,ny+j,k) = divbuf4r_gpu(i,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (ndim==3) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,ng
   do j=1,ny
    do i=1,nx
     fhat_gpu(6,i,j,k-ng) = divbuf5r_gpu(i,j,k)
     fhat_gpu(6,i,j,nz+k) = divbuf6r_gpu(i,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
!
end subroutine bcswapdiv
!
subroutine bcswapduc_prepare
!
! Preparing arrays for data exchange (flow ducros)
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,m
!
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ny
   do i=1,ng
    ducbuf1s_gpu(i,j,k) = ducros_gpu(i,j,k)
    ducbuf2s_gpu(i,j,k) = ducros_gpu(nx-ng+i,j,k)
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ng
   do i=1,nx
    ducbuf3s_gpu(i,j,k) = ducros_gpu(i,ny-ng+j,k)
    ducbuf4s_gpu(i,j,k) = ducros_gpu(i,j,k)
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 if (ndim==3) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,ng
   do j=1,ny
    do i=1,nx
     ducbuf5s_gpu(i,j,k) = ducros_gpu(i,j,k)
     ducbuf6s_gpu(i,j,k) = ducros_gpu(i,j,nz-ng+k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
!
end subroutine bcswapduc_prepare

subroutine bcswapduc
!
! Apply virtual interface boundary conditions (flow ducros)
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,m
 integer :: indx,indy,indz
!
 indx = ng*ny*nz
 indy = nx*ng*nz
 indz = nx*ny*ng
!
#if defined(USE_CUDA) && defined(NOCUDAAWAREMPI)
 iermpi = cudaMemcpyAsync(ducbuf1s, ducbuf1s_gpu, indx, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(ducbuf2s, ducbuf2s_gpu, indx, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(ducbuf3s, ducbuf3s_gpu, indy, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(ducbuf4s, ducbuf4s_gpu, indy, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(ducbuf5s, ducbuf5s_gpu, indz, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaMemcpyAsync(ducbuf6s, ducbuf6s_gpu, indz, cudaMemcpyDeviceToHost, stream2)
 iermpi = cudaStreamSynchronize(stream2)
 call mpi_sendrecv(ducbuf1s,indx,mpi_logical,ileftx ,1,ducbuf2r,indx,mpi_logical,irightx,1,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(ducbuf2s,indx,mpi_logical,irightx,2,ducbuf1r,indx,mpi_logical,ileftx ,2,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(ducbuf3s,indy,mpi_logical,ilefty ,3,ducbuf4r,indy,mpi_logical,irighty,3,mp_carty,istatus,iermpi)       
 call mpi_sendrecv(ducbuf4s,indy,mpi_logical,irighty,4,ducbuf3r,indy,mpi_logical,ilefty ,4,mp_carty,istatus,iermpi)       
 if (ndim==3) then
  call mpi_sendrecv(ducbuf5s,indz,mpi_logical,ileftz ,5,ducbuf6r,indz,mpi_logical,irightz,5,mp_cartz,istatus,iermpi)       
  call mpi_sendrecv(ducbuf6s,indz,mpi_logical,irightz,6,ducbuf5r,indz,mpi_logical,ileftz ,6,mp_cartz,istatus,iermpi)       
 endif
 iermpi = cudaMemcpyAsync(ducbuf1r_gpu, ducbuf1r, indx, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(ducbuf2r_gpu, ducbuf2r, indx, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(ducbuf3r_gpu, ducbuf3r, indy, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(ducbuf4r_gpu, ducbuf4r, indy, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(ducbuf5r_gpu, ducbuf5r, indz, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaMemcpyAsync(ducbuf6r_gpu, ducbuf6r, indz, cudaMemcpyHostToDevice, stream2)
 iermpi = cudaStreamSynchronize(stream2)
#else
 call mpi_sendrecv(ducbuf1s_gpu,indx,mpi_logical,ileftx ,1,ducbuf2r_gpu,indx,mpi_logical,irightx,1,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(ducbuf2s_gpu,indx,mpi_logical,irightx,2,ducbuf1r_gpu,indx,mpi_logical,ileftx ,2,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(ducbuf3s_gpu,indy,mpi_logical,ilefty ,3,ducbuf4r_gpu,indy,mpi_logical,irighty,3,mp_carty,istatus,iermpi)       
 call mpi_sendrecv(ducbuf4s_gpu,indy,mpi_logical,irighty,4,ducbuf3r_gpu,indy,mpi_logical,ilefty ,4,mp_carty,istatus,iermpi)       
 if (ndim==3) then
  call mpi_sendrecv(ducbuf5s_gpu,indz,mpi_logical,ileftz ,5,ducbuf6r_gpu,indz,mpi_logical,irightz,5,mp_cartz,istatus,iermpi)       
  call mpi_sendrecv(ducbuf6s_gpu,indz,mpi_logical,irightz,6,ducbuf5r_gpu,indz,mpi_logical,ileftz ,6,mp_cartz,istatus,iermpi)       
 endif
#endif
!
 if (ileftx/=mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do i=1,ng
     ducros_gpu(i-ng,j,k) = ducbuf1r_gpu(i,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (irightx/=mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do i=1,ng
     ducros_gpu(nx+i,j,k) = ducbuf2r_gpu(i,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (ilefty/=mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ng
    do i=1,nx
     ducros_gpu(i,j-ng,k) = ducbuf3r_gpu(i,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (irighty/=mpi_proc_null) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ng
    do i=1,nx
     ducros_gpu(i,ny+j,k) = ducbuf4r_gpu(i,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
 if (ndim==3) then
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,ng
   do j=1,ny
    do i=1,nx
     ducros_gpu(i,j,k-ng) = ducbuf5r_gpu(i,j,k)
     ducros_gpu(i,j,nz+k) = ducbuf6r_gpu(i,j,k)
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
!
end subroutine bcswapduc
