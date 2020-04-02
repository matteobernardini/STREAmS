subroutine df_par
!
! Definition of parameters for digital filtering
!
 use mod_streams
 implicit none
!
 integer :: i,ii,j,jj,kk,m
 real(mykind), dimension(3,ny) :: ylen ! Integral lengthscale in y
 real(mykind), dimension(3,ny) :: zlen ! Integral lengthscale in z
 real(mykind) :: sumbx,sumby,sumbz
!      
 do j=1,ny
  zlen(1,j) = min(0.05_mykind+(y(j)/0.3_mykind)*0.20_mykind,0.25_mykind)
  zlen(2,j) = min(0.05_mykind+(y(j)/0.3_mykind)*0.15_mykind,0.20_mykind)
  zlen(3,j) = min(0.05_mykind+(y(j)/0.3_mykind)*0.20_mykind,0.25_mykind)
 enddo
 ylen = 0.7_mykind*zlen
 xlen_df(1) = .50_mykind
 xlen_df(2) = .15_mykind
 xlen_df(3) = .15_mykind
!
 bx_df = 0.
 by_df = 0.
 bz_df = 0.
!
 do m=1,3
!
  do i=1,nxmax
   sumbx = 0._mykind
   do ii=-nfmax,nfmax
    bx_df(m,i,ii) = exp(-pi*abs(ii)/(xlen_df(m)/dxg(i)))
    sumbx = sumbx+bx_df(m,i,ii)**2
   enddo
   sumbx=sqrt(sumbx)
   do ii=-nfmax,nfmax
    bx_df(m,i,ii)=bx_df(m,i,ii)/sumbx
   enddo
  enddo
!
  do j=1,ny
   sumby = 0._mykind
   do jj=-nfmax,nfmax
!   by_df(m,j,jj) =
!    exp(-pi*jj**2/(2._mykind*(ylen(m,j)/(y(j)-y(j-1)))**2))
    by_df(m,j,jj) = exp(-pi*abs(jj)/(ylen(m,j)/(y(j)-y(j-1))))
    sumby = sumby+by_df(m,j,jj)**2
   enddo
   sumby=sqrt(sumby)
   do jj=-nfmax,nfmax
    by_df(m,j,jj)=by_df(m,j,jj)/sumby
   enddo
  enddo
!
  do j=1,ny
   sumbz = 0._mykind
   do kk=-nfmax,nfmax
!   bz_df(m,j,kk) = exp(-pi*kk**2/(2._mykind*(zlen(m,j)/dzg(1))**2))
    bz_df(m,j,kk) = exp(-pi*abs(kk)/(zlen(m,j)/dzg(1)))
    sumbz = sumbz+bz_df(m,j,kk)**2
   enddo
   sumbz=sqrt(sumbz)
   do kk=-nfmax,nfmax
    bz_df(m,j,kk)=bz_df(m,j,kk)/sumbz
   enddo
  enddo
!
 enddo
!
end subroutine df_par
