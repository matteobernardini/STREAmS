subroutine generatewmean
!
! Generating field for initialization of BL and SBLI simulations
!
 use mod_streams
 implicit none
!
 real(mykind), dimension(1-ng:nxmax+ng+1) :: deltavec,deltavvec,cfvec,thvec
 real(mykind), dimension(ny) :: uvec,dvec,tvec,yvec
 real(mykind) :: cf,delta,deltav,retau,retauold,rhoi,th,thrat,ui,vi,vi_j,vi_jm,yl
 integer :: i,j,k,m,ii,jj,jjj
 real(mykind) :: ddrho,ddt,ddu,ti
!
 do j=1,ny
  yvec(j) = y(j)
 enddo
!
! Compute globally deltavec and deltavvec
!
 delta        = 1._mykind
 deltavec(1)  = delta
 deltav       = delta/retauinflow
 deltavvec(1) = deltav
 call meanvelocity(ny,yvec,uvec,dvec,tvec,retauinflow,rm,trat,thrat,cf)
 cfvec(1)     = cf
 thvec(1)     = thrat*deltavec(1)
 retauold     = retauinflow
!
 do i=1,ng
  retau = retauold
  do
   call meanvelocity(ny,yvec,uvec,dvec,tvec,retau,rm,trat,thrat,cf)
   th = thvec(2-i)-0.25_mykind*abs((xg(1-i)-xg(2-i)))*(cf+cfvec(2-i)) ! dthdx evaluated with second order accuracy
   delta  = th/thrat
   deltav = deltavvec(2-i)*sqrt(cfvec(2-i)/cf)
   retau  = delta/deltav
   if (abs(retau-retauold) < 0.01_mykind) exit
   retauold = retau
  enddo
  thvec    (1-i) = th
  cfvec    (1-i) = cf
  deltavec (1-i) = delta
  deltavvec(1-i) = deltav
 enddo
!
 retauold = deltavec(1)/deltavvec(1)
!
 if (masterproc) open(182,file='cfstart.dat')
 do i=2,nxmax+ng+1
  retau = retauold
  do
   call meanvelocity(ny,yvec,uvec,dvec,tvec,retau,rm,trat,thrat,cf)
   th     = thvec(i-1)+0.25_mykind*(xg(i)-xg(i-1))*(cf+cfvec(i-1))
   delta  = th/thrat
   deltav = deltavvec(i-1)*sqrt(cfvec(i-1)/cf)
   retau  = delta/deltav
   if (abs(retau-retauold)<0.01_mykind) exit
   retauold = retau
  enddo
   thvec    (i) = th
   cfvec    (i) = cf
   deltavec (i) = delta
   deltavvec(i) = deltav
   if (masterproc) write(182,100) xg(i),delta,deltav,cf,th
 100  format(20ES20.10)
 enddo
 if (masterproc) close(182)
!
!Compute locally wmean from 1-ng to nx+ng+1
!
 wmean = 0._mykind
 do i=1-ng,nx+ng+1
  ii     = ncoords(1)*nx+i
  delta  = deltavec(ii)
  deltav = deltavvec(ii)
  retau  = delta/deltav
  call meanvelocity(ny,yvec,uvec,dvec,tvec,retau,rm,trat,thrat,cf)
  do j=1,ny
   yl = yvec(j)/delta
   call locateval(yvec,ny,yl,jj) ! yl is between yvec(jj) and yvec(jj+1)
   m = 4
   jjj = min(max(jj-(m-1)/2,1),ny+1-m)
!  call polint(yvec(jjj),dvec(jjj),m,yl,rhoi,ddrho)
!  call polint(yvec(jjj),uvec(jjj),m,yl,ui  ,ddu  )
!  call polint(yvec(jjj),tvec(jjj),m,yl,ti  ,ddt  )
   call pol_int(yvec(jjj),dvec(jjj),m,yl,rhoi)
   call pol_int(yvec(jjj),uvec(jjj),m,yl,ui  )
   call pol_int(yvec(jjj),tvec(jjj),m,yl,ti  )
   wmean(1,i,j) = rhoi
   wmean(2,i,j) = rhoi*ui*u0
  enddo
 enddo
!
 do i=1-ng,nx+ng
  ii = ncoords(1)*nx+i
  do j=2,ny
   vi_j  = -(wmean(2,i+1,j)-wmean(2,i,j))/(xg(ii+1)-xg(ii))
   vi_jm = -(wmean(2,i+1,j-1)-wmean(2,i,j-1))/(xg(ii+1)-xg(ii))
   vi    = 0.5_mykind*(vi_j+vi_jm)
   wmean(3,i,j) = wmean(3,i,j-1)+vi*(yvec(j)-yvec(j-1))
  enddo
 enddo
!
 do j=1,ny
  do i=1-ng,nx+ng
   wmean(5,i,j) = p0*gm+0.5_mykind*(wmean(2,i,j)**2+wmean(3,i,j)**2)/wmean(1,i,j)
  enddo
 enddo
!
end subroutine generatewmean
