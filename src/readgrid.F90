subroutine readgrid
!
! Reading the mesh
!
 use mod_streams
 implicit none
!
 integer :: i,j,k
!
! x
 open(10,file='x.dat')
 do i=1-ng,nxmax+ng+1
  read(10,*) xg(i)
 enddo
!xg(nxmax+ng+1) = 2._mykind*xg(nxmax+ng)-xg(nxmax+ng-1)
 close(10)
! y
 open(10,file='y.dat')
 do j=1-ng,nymax+ng
  read(10,*) yg(j)
 enddo
 close(10)
! z
 if (ndim==3) then
  open(10,file='z.dat')
  do k=1-ng,nzmax+ng
   read(10,*) zg(k)
  enddo
  close(10)
 endif
!
end subroutine readgrid
