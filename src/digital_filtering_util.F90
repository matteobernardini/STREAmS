subroutine trasp_xz_to_yz(arr_xz,arr_yz,ndir,nx_1,ny_1,nz_3,nb1)
!
! Trasposition of data from a xz split to a yz split or viceversa
!
 use mod_streams
 implicit none
!
 integer, intent(in) :: ndir,nx_1,ny_1,nz_3,nb1
 real(mykind), dimension(3,nx_1,ny_1*nb1,nz_3), intent(inout) :: arr_xz
 real(mykind), dimension(3,nx_1*nb1,ny_1,nz_3), intent(inout) :: arr_yz
!
 real(mykind), dimension(3*nb1*nx_1*ny_1*nz_3) :: sendbuf,recvbuf
!
 integer :: i,j,k,l,m,im,jm,ii,jj,ind,nportion
!
 if (ndir==1) then ! Direct transposition from xz to yz
!
  ind = 0
  do m=1,nb1
   jm = (m-1)*ny_1
   do k=1,nz_3
    do j=1,ny_1
     jj = j+jm
     do i=1,nx_1
      do l=1,3
       ind = ind+1
       sendbuf(ind) = arr_xz(l,i,jj,k)
      enddo
     enddo
    enddo
   enddo
  enddo
  nportion = ind/nb1
!
  call mpi_alltoall(sendbuf,nportion,mpi_prec,recvbuf,nportion,mpi_prec,mp_cartx,iermpi)
!
  ind = 0
  do m=1,nb1
   im = (m-1)*nx_1
   do k=1,nz_3
    do j=1,ny_1
     do i=1,nx_1
      ii = i+im
      do l=1,3
       ind = ind + 1
       arr_yz(l,ii,j,k) = recvbuf(ind)
      enddo
     enddo
    enddo
   enddo
  enddo
!
 elseif (ndir==-1) then ! Inverse transposition from yz to xz
!
  ind = 0
  do m=1,nb1
   im = (m-1)*nx_1
   do k=1,nz_3
    do j=1,ny_1
     do i=1,nx_1
      ii = i+im
      do l=1,3
       ind = ind+1
       sendbuf(ind) = arr_yz(l,ii,j,k)
      enddo
     enddo
    enddo
   enddo
  enddo
  nportion = ind/nb1
!
  call mpi_alltoall(sendbuf,nportion,mpi_prec,recvbuf,nportion,mpi_prec,mp_cartx,iermpi)
!
  ind = 0
  do m=1,nb1
   jm = (m-1)*ny_1
   do k=1,nz_3
    do j=1,ny_1
     jj = j+jm
     do i=1,nx_1
      do l=1,3
       ind = ind + 1
       arr_xz(l,i,jj,k) = recvbuf(ind)
      enddo
     enddo
    enddo
   enddo
  enddo
!
 endif
!
end subroutine trasp_xz_to_yz
!
subroutine trasp_yz_to_xy(arr_yz,arr_xy,ndir,nx_3,ny_1,nz_3,nb3)
!
! Trasposition of data from a yz split to a xy split or viceversa
!
 use mod_streams
 implicit none
!
 integer, intent(in) :: ndir,nx_3,ny_1,nz_3,nb3
 real(mykind), dimension(3,nx_3*nb3,ny_1,nz_3), intent(inout) :: arr_yz
 real(mykind), dimension(3,nx_3,ny_1,nz_3*nb3), intent(inout) :: arr_xy
!
 real(mykind), dimension(3*nb3*nx_3*ny_1*nz_3) :: sendbuf,recvbuf
!
 integer :: i,j,k,l,m,im,km,ii,kk,ind,nportion
!
 if (ndir==1) then ! Direct transposition from xz to yz
!
  ind = 0
  do m=1,nb3
   im = (m-1)*nx_3
   do k=1,nz_3
    do j=1,ny_1
     do i=1,nx_3
      ii = i+im
      do l=1,3
       ind = ind+1
       sendbuf(ind) = arr_yz(l,ii,j,k)
      enddo
     enddo
    enddo
   enddo
  enddo
  nportion = ind/nb3
!
  call mpi_alltoall(sendbuf,nportion,mpi_prec,recvbuf,nportion,mpi_prec,mp_cartz,iermpi)
!
  ind = 0
  do m=1,nb3
   km = (m-1)*nz_3
   do k=1,nz_3
    kk = k+km
    do j=1,ny_1
     do i=1,nx_3
      do l=1,3
       ind = ind + 1
       arr_xy(l,i,j,kk) = recvbuf(ind)
      enddo
     enddo
    enddo
   enddo
  enddo
!
 elseif (ndir==-1) then ! Inverse transposition from yz to xz
!
  ind = 0
  do m=1,nb3
   km = (m-1)*nz_3
   do k=1,nz_3
    kk = k+km
    do j=1,ny_1
     do i=1,nx_3
      do l=1,3
       ind = ind+1
       sendbuf(ind) = arr_xy(l,i,j,kk)
      enddo
     enddo
    enddo
   enddo
  enddo
  nportion = ind/nb3
!
  call mpi_alltoall(sendbuf,nportion,mpi_prec,recvbuf,nportion,mpi_prec,mp_cartz,iermpi)
!
  ind = 0
  do m=1,nb3
   im = (m-1)*nx_3
   do k=1,nz_3
    do j=1,ny_1
     do i=1,nx_3
      ii = i+im
      do l=1,3
       ind = ind + 1
       arr_yz(l,ii,j,k) = recvbuf(ind)
      enddo
     enddo
    enddo
   enddo
  enddo
!
 endif
!
end subroutine trasp_yz_to_xy
!
subroutine swap_ranfxz(wxz)
!
 use mod_streams
 implicit none
!
 real(mykind), dimension(3,1-ngdf:nx+ngdf,ny,nz), intent(inout) :: wxz
!
 integer :: i,j,k,m,ind,procleft,procright
 real(mykind), dimension(3,ngdf,ny,nz) :: buf1s,buf2s,buf1r,buf2r
!
 do k=1,nz
  do j=1,ny
   do i=1,ngdf
    do m=1,3
     buf1s(m,i,j,k) = wxz(m,i,j,k)
     buf2s(m,i,j,k) = wxz(m,nx-ngdf+i,j,k)
    enddo
   enddo
  enddo
 enddo
!
 ind = 3*ngdf*ny*nz
!
 procleft  = ncoords(1)-1
 procright = ncoords(1)+1
 if (ncoords(1)==0) procleft = nblocks(1)-1
 if (ncoords(1)==(nblocks(1)-1)) procright = 0
 call mpi_sendrecv(buf1s,ind,mpi_prec,procleft ,1,buf2r,ind,mpi_prec,procright,1,mp_cartx,istatus,iermpi)       
 call mpi_sendrecv(buf2s,ind,mpi_prec,procright,2,buf1r,ind,mpi_prec,procleft ,2,mp_cartx,istatus,iermpi)       
!
 do k=1,nz
  do j=1,ny
   do i=1,ngdf
    do m=1,3
     wxz(m,i-ngdf,j,k) = buf1r(m,i,j,k)
    enddo
   enddo
  enddo
 enddo
 do k=1,nz
  do j=1,ny
   do i=1,ngdf
    do m=1,3
     wxz(m,nx+i,j,k) = buf2r(m,i,j,k)
    enddo
   enddo
  enddo
 enddo
!
end subroutine swap_ranfxz
!
subroutine swapz_rf(rfloc)
!
 use mod_streams
 implicit none
!
 real(mykind), dimension(3,1-nfmax:ny+nfmax,1-nfmax:nz+nfmax), intent(inout) :: rfloc
!
 integer :: i,j,k,m,ind
 real(mykind), dimension(3,1-nfmax:ny+nfmax,nfmax) :: buf5s,buf6s,buf5r,buf6r
!
 do k=1,nfmax
  do j=1-nfmax,ny+nfmax
   do m=1,3
    buf5s(m,j,k) = rfloc(m,j,k)
    buf6s(m,j,k) = rfloc(m,j,nz-nfmax+k)
   enddo
  enddo
 enddo
!
 ind = 3*nfmax*(ny+2*nfmax)
!
 call mpi_sendrecv(buf5s,ind,mpi_prec,ileftz ,5,buf6r,ind,mpi_prec,irightz,5,mp_cartz,istatus,iermpi)
 call mpi_sendrecv(buf6s,ind,mpi_prec,irightz,6,buf5r,ind,mpi_prec,ileftz ,6,mp_cartz,istatus,iermpi)
!
 do k=1,nfmax
  do j=1-nfmax,ny+nfmax
   do m=1,3
    rfloc(m,j,k-nfmax) = buf5r(m,j,k)
    rfloc(m,j,nz+k)    = buf6r(m,j,k)
   enddo
  enddo
 enddo
!
end subroutine swapz_rf
!
subroutine swapz_extended_rf(rfloc)
!
 use mod_streams
 implicit none
!
 real(mykind), dimension(3,1-nfmax:ny+nfmax,1-nfmax:nz+nfmax), intent(inout) :: rfloc
!
 integer :: i,j,k,m,ind,kend,kswap,nplanes,procleft,procright
 real(mykind), dimension(:,:,:), allocatable :: buf5s,buf6s,buf5r,buf6r
!
 kend = 1+int((nfmax-1)/nz)
 do kswap = 1,kend
  nplanes = nz
  if (kswap==kend) nplanes = mod(nfmax-1,nz)+1
  allocate(buf5s(3,1-nfmax:ny+nfmax,nplanes))
  allocate(buf5r(3,1-nfmax:ny+nfmax,nplanes))
  allocate(buf6s(3,1-nfmax:ny+nfmax,nplanes))
  allocate(buf6r(3,1-nfmax:ny+nfmax,nplanes))
!
  do k=1,nplanes
   do j=1-nfmax,ny+nfmax
    do m=1,3
     buf5s(m,j,k) = rfloc(m,j,k)
     buf6s(m,j,k) = rfloc(m,j,nz-nplanes+k)
    enddo
   enddo
  enddo
!
  ind = 3*nplanes*(ny+2*nfmax)
!
  procleft  = ncoords(3)-kswap
  procright = ncoords(3)+kswap
  if (procleft<0) procleft = nblocks(3)+procleft
  if (procright>(nblocks(3)-1)) procright = procright-nblocks(3)
  call mpi_sendrecv(buf5s,ind,mpi_prec,procleft ,5,buf6r,ind,mpi_prec,procright,5,mp_cartz,istatus,iermpi)
  call mpi_sendrecv(buf6s,ind,mpi_prec,procright,6,buf5r,ind,mpi_prec,procleft ,6,mp_cartz,istatus,iermpi)
!
  do k=1,nplanes
   do j=1-nfmax,ny+nfmax
    do m=1,3
     rfloc(m,j,k-nz*(kswap-1)-nplanes) = buf5r(m,j,k)
     rfloc(m,j,nz+nz*(kswap-1)+k) = buf6r(m,j,k)
    enddo
   enddo
  enddo
!
  deallocate(buf5s,buf6s,buf5r,buf6r)
!
 enddo
!
end subroutine swapz_extended_rf
