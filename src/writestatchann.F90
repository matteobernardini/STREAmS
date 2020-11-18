subroutine writestatchann
!
! Writing channel flow statistics
!
 use mod_streams
 implicit none
!
 integer :: j,m
 real(mykind), dimension(nvmean,ny/2) :: w1dh
 real(mykind), dimension(ny/2) :: ufav,vfav,wfav
 real(mykind), dimension(ny/2) :: uvd,ut,yt,dft,ft,gt
 real(mykind), dimension(3)    :: cc
 real(mykind) :: rhow,rmuw,rnuw,deltav
 real(mykind) :: d1,d2
 real(mykind) :: dudyw,tauw,utau,retau
 real(mykind) :: rr,rn
 real(mykind) :: yy,uu2,vv2,ww2,uv,rho2,tt2
 logical, dimension(nvmean) :: symm
!
 cc(1) =  1.83333333333333_mykind
 cc(2) = -1.16666666666667_mykind
 cc(3) =  0.33333333333333_mykind
!
 if (masterproc) then
!
  symm     = .true.
  symm(3 ) = .false.
  symm(14) = .false.
  symm(19) = .false.
!
  do j=1,ny/2
   do m=1,nvmean
   if (symm(m)) then
    w1dh(m,j) = 0.5_mykind*(w_av_1d(m,j)+w_av_1d(m,ny-j+1))
   else
    w1dh(m,j) = 0.5_mykind*(w_av_1d(m,j)-w_av_1d(m,ny-j+1))
   endif
  enddo
 enddo
!
! Wall density and viscosity
  rhow = cc(1)*w1dh( 1,1)+cc(2)*w1dh( 1,2)+cc(3)*w1dh( 1,3)
  rmuw = cc(1)*w1dh(20,1)+cc(2)*w1dh(20,2)+cc(3)*w1dh(20,3)
  rnuw = rmuw/rhow
!
  do j=1,ny/2
   ufav(j) = w1dh(13,j)/w1dh(1,j)
   vfav(j) = w1dh(14,j)/w1dh(1,j)
   wfav(j) = w1dh(15,j)/w1dh(1,j)
  enddo
  d1 = y(1)+1._mykind
  d2 = y(2)+1._mykind
  dudyw = (ufav(1)*d2**2-ufav(2)*d1**2)/(d1*d2*(d2-d1))
  tauw  = rmuw*dudyw
  utau  = sqrt(tauw/rhow)
  retau = rhow*utau/rmuw
!
! Van Driest velocity
  uvd(1) = ufav(1)*sqrt(w1dh(1,1)/rhow) 
  do j=2,ny/2
   uvd(j) = uvd(j-1)+(ufav(j)-ufav(j-1))*0.5_mykind*(sqrt(w1dh(1,j)/rhow)+sqrt(w1dh(1,j-1)/rhow))
  enddo
!
! Trettel and Larsson velocity
  do j=1,ny/2
   rr    = sqrt(w1dh(1,j)/rhow)
   rn    = (w1dh(20,j)/w1dh(1,j)/rnuw)
   ft(j) = (y(j)+1._mykind)/(rr*rn)
  enddo
!
  do j=2,ny/2-1
!  dft(j)   = 0.5_mykind*(ft(j+1)-ft(j-1))*detady(j)
   dft(j)   = (ft(j+1)-ft(j-1))/(y(j+1)-y(j-1))
  enddo
! dft(1)    = (-1.5_mykind*ft(1)+2._mykind*ft(2)-0.5_mykind*ft(3))*detady(1)
  dft(1)    = (-1.5_mykind*ft(1)+2._mykind*ft(2)-0.5_mykind*ft(3))/&
              (-1.5_mykind* y(1)+2._mykind* y(2)-0.5_mykind* y(3))
! dft(ny/2) = ( 0.5_mykind*ft(ny/2-2)-2._mykind*ft(ny/2-1)+1.5_mykind*ft(ny/2))*detady(ny/2)
  dft(ny/2) = ( 0.5_mykind*ft(ny/2-2)-2._mykind*ft(ny/2-1)+1.5_mykind*ft(ny/2))/&
              ( 0.5_mykind* y(ny/2-2)-2._mykind* y(ny/2-1)+1.5_mykind* y(ny/2))
!
  do j=1,ny/2
   gt(j) = dft(j)*w1dh(1,j)/rhow*w1dh(20,j)/w1dh(1,j)/rnuw
  enddo
!
  ut(1) = ufav(1)*gt(1)
  yt(1) = (1._mykind+y(1))*dft(1)
  do j=2,ny/2
   ut(j) = ut(j-1)+(ufav(j)-ufav(j-1))*0.5_mykind*(gt (j)+gt (j-1))
   yt(j) = yt(j-1)+(y   (j)-y   (j-1))*0.5_mykind*(dft(j)+dft(j-1))
  enddo
!
  open(unit=10,file='channstat.prof',form='formatted')
  deltav = 1._mykind/retau
  do j=1,ny/2
   yy   = 1.+y(j)
   uu2  = w1dh(16,j)/w1dh(1,j) - ufav(j)*ufav(j)
   vv2  = w1dh(17,j)/w1dh(1,j) - vfav(j)*vfav(j)
   ww2  = w1dh(18,j)/w1dh(1,j) - wfav(j)*wfav(j)
   uv   = w1dh(19,j)/w1dh(1,j) - ufav(j)*vfav(j)
   rho2 = w1dh(7 ,j) - w1dh(1,j)**2
   tt2  = w1dh(12,j) - w1dh(6,j)**2
   write(10,100) yy, &                !y/h
                 yy/deltav,&          !y^+
                 yt(j)/deltav,&       !y_T^+
                 ufav(j)/u0,&         !u/u_b
                 ufav(j)/utau,&       !u^+
                 uvd(j)/utau,&        !u_VD^+
                 ut(j)/utau,&         !u_T^+
                 w1dh(1,j)/rhow,&     !rho/rho_w
                 uu2/utau**2,&        !u2^+
                 vv2/utau**2,&        !v2^+
                 ww2/utau**2,&        !w2^+
                 uv /utau**2,&        !uv^+
                 w1dh(6,j),  &        !T/T_w
                 rho2/rhow,  &        !rho2/rhow
                 tt2                  !T2/T_w

  enddo
  close(10)
 endif
!
 100  format(200ES20.10)
!
end subroutine writestatchann
