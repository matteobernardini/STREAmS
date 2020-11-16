subroutine readgrid
!
! Reading the mesh
!
 use mod_streams
 implicit none
!
 integer :: i,j,k
 integer :: i1,i2,j1,j2,k1,k2
!
! x
 i1 = 1-ng
 i2 = nxmax+ng+1
 if (iflow==-1) then
  i1 = 1
  i2 = nxmax
 endif
 open(10,file='x.dat')
 do i=i1,i2
  read(10,*) xg(i)
 enddo
 close(10)
 if (iflow==-1) then
  do i=1,ng
   xg(1-i)     = 2._mykind*xg(2-i)-xg(3-i)
   xg(nxmax+i) = 2._mykind*xg(nxmax+i-1)-xg(nxmax+i-2)
  enddo
  xg(nxmax+ng+1) = 2._mykind*xg(nxmax+ng)-xg(nxmax+ng-1)
 endif
! y
 open(10,file='y.dat')
 j1 = 1-ng
 j2 = nymax+ng
 if (iflow==-1) then
  j1 = 1
  j2 = nymax
 endif
 do j=j1,j2
  read(10,*) yg(j)
 enddo
 close(10)
 if (iflow==-1) then
  do j=1,ng
   yg(1-j)     = 2._mykind*yg(2-j)-yg(3-j)
   yg(nymax+j) = 2._mykind*yg(nymax+j-1)-yg(nymax+j-2)
  enddo
 endif
! z
 if (ndim==3) then
  k1 = 1-ng
  k2 = nzmax+ng
  if (iflow==-1) then
   k1 = 1
   k2 = nzmax
  endif
  open(10,file='z.dat')
  do k=k1,k2
   read(10,*) zg(k)
  enddo
  close(10)
  if (iflow==-1) then
   do k=1,ng
    zg(1-k)     = 2._mykind*zg(2-k)-zg(3-k)
    zg(nzmax+k) = 2._mykind*zg(nzmax+k-1)-zg(nzmax+k-2)
   enddo
  endif
 endif
!
 if (iflow==0) then
  open(10,file='yn.dat')
  do j=1,nymax+1
   read(10,*) yn(j)
  enddo
  close(10)
 endif
!
end subroutine readgrid
