subroutine target_reystress
!
! Interpolate DNS database to obtain target Reynolds stress for digital filtering 
!
 use mod_streams
 implicit none
!
 integer :: nstatbl,nydatamax
 integer :: l,j,jj,jjj,m,mm,ltau,lltau
 real(mykind) :: yy,fyy
 real(mykind) :: uui,vvi,wwi,uvi
 real(mykind) :: duui,dvvi,dwwi,duvi
 real(mykind) :: rhofac
 integer, dimension(:), allocatable :: nyvec
 real(mykind), dimension(:), allocatable :: retauvec
 real(mykind), dimension(:,:,:), allocatable :: reydata
 real(mykind), dimension(:,:,:), allocatable :: reydataycode
 real(mykind), dimension(3,3,ny) :: reytarget
!
 open(12,file='database_bl.dat')
 read(12,*) nstatbl,nydatamax
 allocate(nyvec(nstatbl))
 allocate(retauvec(nstatbl))
 allocate(reydata(nstatbl,nydatamax,5))
 allocate(reydataycode(nstatbl,ny,5))
 do l=1,nstatbl
  read(12,*) nyvec(l),retauvec(l)
  do j=1,nyvec(l)
   read(12,*) (reydata(l,j,m),m=1,5)
  enddo
 enddo
 close(12)
!
 reydataycode = 0._mykind
 mm = 2
 do l=1,nstatbl
  do j=2,ny
   yy = y(j)
   call locateval(reydata(l,1:nyvec(l),1),nyvec(l),yy,jj)
   jjj = min(max(jj-(mm-1)/2,1),nyvec(l)+1-mm)
   call pol_int(reydata(l,jjj:jjj+mm-1,1),reydata(l,jjj:jjj+mm-1,2),mm,yy,uui)
   call pol_int(reydata(l,jjj:jjj+mm-1,1),reydata(l,jjj:jjj+mm-1,3),mm,yy,vvi)
   call pol_int(reydata(l,jjj:jjj+mm-1,1),reydata(l,jjj:jjj+mm-1,4),mm,yy,wwi)
   call pol_int(reydata(l,jjj:jjj+mm-1,1),reydata(l,jjj:jjj+mm-1,5),mm,yy,uvi)
   reydataycode(l,j,1) = yy
   fyy = 0.5*(1-tanh((yy-2.)/.1))
   reydataycode(l,j,2) = uui*fyy
   reydataycode(l,j,3) = vvi*fyy
   reydataycode(l,j,4) = wwi*fyy
   reydataycode(l,j,5) = uvi*fyy
  enddo
 enddo
!
 reytarget = 0._mykind
 open(12,file='rey_imposed.dat')
 call locateval(retauvec,nstatbl,retauinflow,ltau) ! retau is between retauvec(ltau) and retauvec(ltau+1)
 mm = 2
 lltau = min(max(ltau-(mm-1)/2,1),nstatbl+1-mm)
 do j=2,ny
  rhofac = wmean(1,1,j)/wmean(1,1,1)
  call pol_int(retauvec(lltau),reydataycode(lltau:lltau+mm-1,j,2),mm,retauinflow,uui)
  reytarget(1,1,j) = uui/rhofac
  call pol_int(retauvec(lltau),reydataycode(lltau:lltau+mm-1,j,3),mm,retauinflow,vvi)
  reytarget(2,2,j) = vvi/rhofac
  call pol_int(retauvec(lltau),reydataycode(lltau:lltau+mm-1,j,4),mm,retauinflow,wwi)
  reytarget(3,3,j) = wwi/rhofac
  call pol_int(retauvec(lltau),reydataycode(lltau:lltau+mm-1,j,5),mm,retauinflow,uvi)
  reytarget(1,2,j) = uvi/rhofac
  reytarget(2,1,j) = uvi/rhofac
  reytarget(1,1,j) = max(reytarget(1,1,j),0.000000001_mykind)
  reytarget(2,2,j) = max(reytarget(2,2,j),0.000000001_mykind)
  reytarget(3,3,j) = max(reytarget(3,3,j),0.000000001_mykind)
  write(12,100) y(j),(reytarget(l,l,j),l=1,3),reytarget(1,2,j)
 100 format(20ES20.10)
 enddo
 close(12)
!
 amat_df = 0._mykind
 do j=2,ny
  amat_df(1,1,j) = sqrt(reytarget(1,1,j))
  amat_df(2,1,j) = reytarget(2,1,j)/amat_df(1,1,j)
  amat_df(2,2,j) = sqrt(reytarget(2,2,j)-amat_df(2,1,j)**2)
  amat_df(3,1,j) = reytarget(3,1,j)/amat_df(1,1,j)
  amat_df(3,2,j) = (reytarget(3,2,j)-amat_df(2,1,j)*amat_df(3,1,j))/amat_df(2,2,j)
  amat_df(3,3,j) = sqrt(reytarget(3,3,j)-amat_df(3,1,j)**2-amat_df(3,2,j)**2)
 enddo
!
end subroutine target_reystress
