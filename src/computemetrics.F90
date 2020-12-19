subroutine computemetrics
!
! Computing mesh metrics
!
 use mod_streams
 implicit none
!
 integer, dimension(6) :: itag
!
 real(mykind), dimension(nxmax) :: d2xg
 real(mykind), dimension(nymax) :: d2yg
 real(mykind), dimension(nzmax) :: d2zg
 real(mykind), dimension(0:nxmax) :: dxhg
 real(mykind), dimension(0:nymax) :: dyhg
 real(mykind), dimension(0:nzmax) :: dzhg
!
 real(mykind), dimension(4)   :: c
 real(mykind), dimension(0:4) :: cc
 real(mykind), dimension(4)   :: cs
!
 integer :: i,ii,j,jj,k,kk,l,mm,iend,jend,kend
 real(mykind) :: xloc
!
! z
 if (ndim==3) then
 else
  zg   = 0._mykind
  dzg  = 0._mykind
  d2zg = 0._mykind
 endif
!
 rlx = xg(nxmax)-xg(1)
 if (iflow==0) rlx = xg(nxmax+1)-xg(1)
 rly = yg(nymax)-yg(1)
 if (ndim==3) then
  rlz = zg(nzmax+1)-zg(1)
 else
  rlz = 1._mykind
 endif
!
! local coordinates (nodes)
 ii = nx*ncoords(1)
 do i=1-ng,nx+ng
  x(i) = xg(ii+i)
 enddo
 jj = ny*ncoords(2)
 do j=1-ng,ny+ng
  y(j) = yg(jj+j)
 enddo
 kk = nz*ncoords(3)
 do k=1-ng,nz+ng
  z(k) = zg(kk+k)
 enddo
!
! Defining tags for boundary conditions
 itag   = ibc
 ibc    = 0
 iend   = nblocks(1)-1
 jend   = nblocks(2)-1
 kend   = nblocks(3)-1
 if (ncoords(1)==   0) ibc(1) = itag(1)
 if (ncoords(1)==iend) ibc(2) = itag(2)
 if (ncoords(2)==   0) ibc(3) = itag(3)
 if (ncoords(2)==jend) ibc(4) = itag(4)
 if (ncoords(3)==   0) ibc(5) = itag(5)
 if (ncoords(3)==kend) ibc(6) = itag(6)
!
! Evaluation of metric terms
 mm = ivis/2
 select case (mm) ! coefficient for first derivatives
 case (1)
  c(1) = 0.5_mykind
 case (2)
  c(1) = 2._mykind/3._mykind
  c(2) = -1._mykind/12._mykind
 case (3)
  c(1) = 0.75_mykind
  c(2) = -0.15_mykind
  c(3) = 1._mykind/60._mykind
 case (4)
  c(1) = 0.8_mykind
  c(2) = -0.2_mykind
  c(3) = 4._mykind/105._mykind
  c(4) = -1._mykind/280._mykind
 endselect

 select case (mm) ! coefficient for second derivatives
 case (1)
  cc(0) = -2._mykind
  cc(1) =  1._mykind
 case (2)
  cc(0) = -2.5_mykind
  cc(1) = 4._mykind/3._mykind
  cc(2) = -1._mykind/12._mykind
 case (3)
  cc(0) = -245._mykind/90._mykind
  cc(1) = 1.5_mykind
  cc(2) = -0.15_mykind
  cc(3) = 1._mykind/90._mykind
 case (4)
  cc(0) = -25625._mykind/9000._mykind
  cc(1) = 1.6_mykind
  cc(2) = -0.2_mykind
  cc(3) = 253968._mykind/9999990._mykind
  cc(4) = -17857125._mykind/9999990000._mykind
 endselect

 dxg = 0._mykind
 do i=1,nxmax
  do l=1,mm
   dxg(i) = dxg(i)+c(l)*(xg(i+l)-xg(i-l))
  enddo
  d2xg(i) = cc(0)*xg(i)
  do l=1,mm
   d2xg(i) = d2xg(i)+cc(l)*(xg(i+l)+xg(i-l))
  enddo
 enddo

 dyg = 0._mykind
 do j=1,nymax
  do l=1,mm
   dyg(j) = dyg(j)+c(l)*(yg(j+l)-yg(j-l))
  enddo
  d2yg(j) = cc(0)*yg(j)
  do l=1,mm
   d2yg(j) = d2yg(j)+cc(l)*(yg(j+l)+yg(j-l))
  enddo
 enddo

 if (ndim==3) then
  dzg = 0._mykind
  do k=1,nzmax
   do l=1,mm
    dzg(k) = dzg(k)+c(l)*(zg(k+l)-zg(k-l))
   enddo
   d2zg(k) = cc(0)*zg(k)
   do l=1,mm
    d2zg(k) = d2zg(k)+cc(l)*(zg(k+l)+zg(k-l))
   enddo
  enddo
 endif
!
 ii = nx*ncoords(1)
 do i=1,nx
  dcsidx (i) = 1._mykind/(dxg(ii+i))
  dcsidxs(i) = dcsidx(i)*dcsidx(i)
  dcsidx2(i) = -d2xg(ii+i)*dcsidxs(i)
 enddo
 jj = ny*ncoords(2)
 do j=1,ny
  detady (j) = 1._mykind/dyg(jj+j)
  detadys(j) = detady(j)*detady(j)
  detady2(j) = -d2yg(jj+j)*detadys(j)
 enddo
 if (ndim==3) then
  kk = nz*ncoords(3)
  do k=1,nz
   dzitdz (k) = 1._mykind/dzg(kk+k)
   dzitdzs(k) = dzitdz(k)*dzitdz(k)
   dzitdz2(k) = -d2zg(kk+k)*dzitdzs(k)
  enddo
 else
  do k=1,nz
   dzitdz (k) = 0._mykind
   dzitdzs(k) = 0._mykind
   dzitdz2(k) = 0._mykind
  enddo
 endif
!
 if (iflow>0) then
  nstatloc = 0
  do i=1,nstat
   xloc = xstat(i)
   call locateval(x(1:nx+1),nx+1,xloc,ii)
   if (ii>=1.and.ii<=nx) then
    nstatloc = nstatloc+1
    ixstat(nstatloc) = ii
    igxstat(nstatloc) = i
   endif
  enddo
 endif
!
 if (masterproc) then
  open(18,file='dxg.dat')
  do i=1,nxmax
   write(18,*) xg(i),dxg(i)
  enddo
  close(18)
  open(18,file='dyg.dat')
  do j=1,nymax
   write(18,*) yg(j),dyg(j)
  enddo
  close(18)
  if (ndim==3) then
   open(18,file='dzg.dat')
   do k=1,nzmax
    write(18,*) zg(k),dzg(k)
   enddo
   close(18)
  endif
 endif
!
! New coordinates and metrics for improved viscous fluxes
! Coordinates of faces
!
 do i=0,nxmax
  xgh(i) = 0.5_mykind*(xg(i)+xg(i+1))
 enddo
 do j=0,nymax
  ygh(j) = 0.5_mykind*(yg(j)+yg(j+1))
 enddo
 do k=0,nzmax
  zgh(k) = 0.5_mykind*(zg(k)+zg(k+1))
 enddo
!
! local coordinates (faces)
 ii = nx*ncoords(1)
 do i=0,nx
  xh(i) = xgh(ii+i)
 enddo
 jj = ny*ncoords(2)
 do j=0,ny
  yh(j) = ygh(jj+j)
 enddo
 kk = nz*ncoords(3)
 do k=0,nz
  zh(k) = zgh(kk+k)
 enddo
!
! Evaluation of metric terms
 mm = ivis/2
 select case (mm) ! coefficient for first derivatives
 case (1)
  cs(1) = 1._mykind
 case (2)
  cs(1) =   5._mykind/4._mykind
  cs(2) = -1._mykind/12._mykind
 case (3)
  cs(1) = 49._mykind/36._mykind
  cs(2) = -5._mykind/36._mykind
  cs(3) =  1._mykind/90._mykind
 case (4)
  call fail_input("new viscous fluxes not implemented for this order")
 endselect
!
 dxhg = 0._mykind
 do i=0,nxmax
  do l=1,mm
   dxhg(i) = dxhg(i)+cs(l)*(xg(i+l)-xg(i-l+1))
  enddo
 enddo
 dyhg = 0._mykind
 do j=0,nymax
  do l=1,mm
   dyhg(j) = dyhg(j)+cs(l)*(yg(j+l)-yg(j-l+1))
  enddo
 enddo
 dzhg = 0._mykind
 do k=0,nzmax
  do l=1,mm
   dzhg(k) = dzhg(k)+cs(l)*(zg(k+l)-zg(k-l+1))
  enddo
 enddo
 ii = nx*ncoords(1)
 do i=0,nx
  dcsidxh(i) = 1._mykind/(dxhg(ii+i))
 enddo
 jj = ny*ncoords(2)
 do j=0,ny
  detadyh(j) = 1._mykind/(dyhg(jj+j))
 enddo
 kk = nz*ncoords(3)
 do k=0,nz
  dzitdzh(k) = 1._mykind/(dzhg(kk+k))
 enddo
!
end subroutine computemetrics
