subroutine recyc
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,m,indx
 integer, dimension(4) :: req
 integer :: kcoordsendto1,kcoordsendto2
 integer :: kcoordrecvfrom1,kcoordrecvfrom2
 integer :: sendto1,sendto2
 integer :: recvfrom1,recvfrom2
 integer :: kshiftglob
 integer :: n1_start_send,n1_end_send,n2_start_send,n2_end_send
 integer :: n1_start_recv,n1_end_recv,n2_start_recv,n2_end_recv
!
 kshiftglob    = nzmax/2 ! global shift in the spanwise direction (between 0 and nzmax-1)
 n1_start_send = 1 
 n1_end_send   = nz-mod(kshiftglob,nz)
 n2_start_send = n1_end_send+1
 n2_end_send   = nz
 n1_start_recv = 1+mod(kshiftglob,nz)
 n1_end_recv   = nz
 n2_start_recv = 1
 n2_end_recv   = mod(kshiftglob,nz)
!
 req = mpi_request_null
!
 if (ncoords(1)==ibrecyc) then ! Send data
  kcoordsendto1 = ncoords(3)+kshiftglob/nz
  kcoordsendto2 = kcoordsendto1+1
  kcoordsendto1 = mod(kcoordsendto1,nblocks(3))
  kcoordsendto2 = mod(kcoordsendto2,nblocks(3))
  call mpi_cart_rank(mp_cart,[0,0,kcoordsendto1],sendto1,iermpi)
  call mpi_cart_rank(mp_cart,[0,0,kcoordsendto2],sendto2,iermpi)
! print *, 'I am', ncoords(1),ncoords(3),kcoordsendto1,kcoordsendto2,n1_start_send,n1_end_send,n2_start_send,n2_end_send
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do i=1,ng
     do m=1,nv
      wbuf1s_gpu(i,j,k,m) = w_gpu(irecyc+1-i,j,k,m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
  indx = nv*ng*ny*nz
  call mpi_isend(wbuf1s_gpu,indx,mpi_prec,sendto1,2000,mp_cart,req(1),iermpi)
  call mpi_isend(wbuf1s_gpu,indx,mpi_prec,sendto2,3000,mp_cart,req(2),iermpi)
!! call mpi_ssend(wbuf1s_gpu,indx,mpi_prec,0,2000,mp_cartx,iermpi)
 endif
 if (ncoords(1)==0) then ! Receive data
  kcoordrecvfrom1 = ncoords(3)-kshiftglob/nz+nblocks(3)
  kcoordrecvfrom2 = kcoordrecvfrom1-1
  kcoordrecvfrom1 = mod(kcoordrecvfrom1,nblocks(3))
  kcoordrecvfrom2 = mod(kcoordrecvfrom2,nblocks(3))
  call mpi_cart_rank(mp_cart,[ibrecyc,0,kcoordrecvfrom1],recvfrom1,iermpi)
  call mpi_cart_rank(mp_cart,[ibrecyc,0,kcoordrecvfrom2],recvfrom2,iermpi)
! print *, 'I am recv', ncoords(1),ncoords(3),kcoordrecvfrom1,kcoordrecvfrom2,n1_start_recv,n1_end_recv,n2_start_recv,n2_end_recv
  indx = nv*ng*ny*nz
  call mpi_irecv(wbuf1r_gpu,indx,mpi_prec,recvfrom1,2000,mp_cart,req(3),iermpi)
  call mpi_irecv(wbuf2r_gpu,indx,mpi_prec,recvfrom2,3000,mp_cart,req(4),iermpi)
 endif
 call mpi_waitall(4,req,mpi_statuses_ignore,iermpi)
 if (ncoords(1)==0) then
 !$cuf kernel do(2) <<<*,*>>>
  do k=n1_start_recv,n1_end_recv
   do j=1,ny
    do i=1,ng
     do m=1,nv
      wrecyc_gpu(i,j,k,m) = wbuf1r_gpu(i,j,k-n1_start_recv+n1_start_send,m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 !$cuf kernel do(2) <<<*,*>>>
  do k=n2_start_recv,n2_end_recv
   do j=1,ny
    do i=1,ng
     do m=1,nv
      wrecyc_gpu(i,j,k,m) = wbuf2r_gpu(i,j,k-n2_start_recv+n2_start_send,m)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 endif
!
end subroutine recyc
!
subroutine recyc_prepare
!
 use mod_streams
 implicit none
!
 integer :: igrecyc
!
 call locateval(xg(1:nxmax),nxmax,xrecyc,igrecyc) ! xrecyc is between xg(ii) and xg(ii+1), ii is between 0 and nxmax
 ibrecyc = (igrecyc-1)/nx
 irecyc  = igrecyc-nx*ibrecyc
 if (irecyc<ng) then
  irecyc  = ng
  igrecyc = irecyc+nx*ibrecyc
 endif
 if (masterproc) write(*,*) 'Recycling station exactly at = ', xg(igrecyc)
!
end subroutine recyc_prepare
