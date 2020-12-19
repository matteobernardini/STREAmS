subroutine generategrid
!
!Generating mesh (x.dat, y.dat and z.dat)
!
!
! Mesh is generated in this subroutine for standard cases (iflow = 0,1,2).
! When iflow = -1, an external program is used for mesh generation
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,l
 integer :: nygp
 real(mykind) :: dx,dy,dz,dcsi,csi,dy1
 real(mykind) :: b,bold
 real(mykind) :: rlygp,rgp,rgpold
 real(mykind) :: rlywrold
!
 real(mykind), dimension(nymax+ng) :: tmpyg
!
! Uniform spacing in x and z
!
 dx = rlx/(nxmax-1)
 if (iflow==0) dx = rlx/nxmax
 do i=1-ng,nxmax+ng+1 ! One extra ghost node needed for bl cases for estimation of initial wall-normal velocity by means of continuity equation
  xg(i) = (i-1)*dx
 enddo
!
 if (ndim==3) then
  dz = rlz/nzmax
  do k=1-ng,nzmax+ng
   zg(k) = (k-1)*dz
  enddo
 endif
!
! Wall-normal distribution is case dependent
!
 select case (iflow)
!
 case(-1) ! wind tunnel
  call readgrid()
 case(0)  ! Channel Flow
!
  dcsi = 1._mykind/nymax
  b    = 1._mykind
  do
   bold = b
   b = atanh(tanh(0.5_mykind*b)*(dyp_target/retauinflow-1._mykind))/(dcsi-0.5_mykind)
   if (abs(b-bold) < tol_iter) exit
  enddo
  if (masterproc) write(*,*) 'Stretching parameter b =', b
  do j=1,ny+1
   csi   = (real(j,mykind)-1._mykind)*dcsi
   yn(j) = tanh(b*(csi-0.5_mykind))/tanh(b*0.5_mykind)
  enddo
  do j=1,ny
   yg(j) = 0.5_mykind*(yn(j)+yn(j+1))
  enddo
  do j=1,ng
   yg(nymax+j) =  2._mykind-yg(nymax+1-j)
   yg(1-j)     = -2._mykind-yg(j)
  enddo
!
 case(1) ! BL
!
! Create sinh mesh from 0 up to rlywr, then geometric progression
!
  dcsi = 1._mykind/(nymaxwr-1)
!
  b = 1._mykind
  do
   bold = b
   b = asinh(sinh(b*dcsi)*rlywr*retauinflow/dyp_target)
   if (abs(b-bold) < tol_iter) exit
!  write(*,*) 'Stretching parameter b =', b 
  enddo
  if (masterproc) write(*,*) 'Stretching parameter b =', b 
!
  do j=1,nymaxwr+1
   csi = (j-1)*dcsi
   yg(j) = rlywr*sinh(b*csi)/sinh(b)
  enddo
  dy1 = yg(nymaxwr+1)-yg(nymaxwr)
  if (masterproc) write(*,*) 'Delta y+_w = ', yg(2)*retauinflow
!
  rlygp = rly-rlywr
  nygp  = nymax-nymaxwr
  rgp = 1.01_mykind
  do
   rgpold = rgp
   rgp = 1._mykind-(1._mykind-rgp)/dy1*rlygp
   rgp = rgp**(1._mykind/nygp)
   if (abs(rgp-rgpold) < tol_iter) exit
!  write(*,*) 'Geometric progression factor rgp =', rgp
  enddo
  if (masterproc) write(*,*) 'Geometric progression factor rgp =', rgp
!
  do j=1,nygp+ng
   yg(nymaxwr+j) = yg(nymaxwr+j-1) + dy1*rgp**(j-1._mykind)
  enddo
!
  do l=1,nsmoo
   tmpyg = yg(1:nymax+ng)
   do j=2,nymax+ng-1
    yg(j) = tmpyg(j-1)+4._mykind*tmpyg(j)+tmpyg(j+1)
    yg(j) = yg(j)/6._mykind
   enddo
  enddo
!
  do j=1,ng
   yg(1-j) = -yg(1+j)
  enddo
!
 case(2) ! SBLI
!
! Create sinh mesh from 0 up to rlywr (to be identified in order to have yg(nymax) = rly), then constant spacing
!
  do
   rlywrold = rlywr
  
   dcsi = 1._mykind/(nymaxwr-1)
!
   b = 1._mykind
   do
    bold = b
    b = asinh(sinh(b*dcsi)*rlywr*retauinflow/dyp_target)
    if (abs(b-bold) < tol_iter) exit
!   write(*,*) 'Stretching parameter b =', b 
   enddo
!
   do j=1,nymaxwr+1
    csi = (j-1)*dcsi
    yg(j) = rlywr*sinh(b*csi)/sinh(b)
   enddo
   dy1 = yg(nymaxwr+1)-yg(nymaxwr)
!
   nygp  = nymax-nymaxwr
!
   do j=1,nygp+ng
    yg(nymaxwr+j) = yg(nymaxwr+j-1) + dy1
   enddo
!
   do l=1,nsmoo
    tmpyg = yg(1:nymax+ng)
    do j=2,nymax+ng-1
     yg(j) = tmpyg(j-1)+4._mykind*tmpyg(j)+tmpyg(j+1)
     yg(j) = yg(j)/6._mykind
    enddo
   enddo
!
   do j=1,ng
    yg(1-j) = -yg(1+j)
   enddo
!
   rlywr = rlywr*rly/yg(nymax)
   if (abs(rlywr-rlywrold) < tol_iter) exit
!
  enddo
  if (masterproc) write(*,*) 'Stretching parameter b =', b 
  if (masterproc) write(*,*) 'Delta y+_w = ', yg(2)*retauinflow
!
 end select
!
 if (masterproc) then
  if (iflow>=0) then
   write(*,*) 'Delta x+ = ', dx*retauinflow
   open(18,file='x.dat')
   do i=1-ng,nxmax+ng+1
    write(18,*) xg(i)
   enddo
   close(18)
   open(18,file='y.dat')
   do j=1-ng,nymax+ng
    write(18,*) yg(j)
   enddo
   close(18)
   write(*,*) 'Delta z+ = ', dz*retauinflow
   open(18,file='z.dat')
   do k=1-ng,nzmax+ng
    write(18,*) zg(k)
   enddo
   close(18)
   if (iflow==0) then
    open(18,file='yn.dat')
    do j=1,nymax+1
     write(18,*) yn(j)
    enddo
    close(18)
   endif
  endif
 endif
!
end subroutine generategrid
