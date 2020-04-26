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
 real(mykind) :: zlenin_u,zlenin_v,zlenin_w
 real(mykind) :: zlenou_u,zlenou_v,zlenou_w
 real(mykind) :: ftany
!      
!Outer z scale for u,v,w
 zlenou_u =  0.40_mykind
 zlenou_v =  0.30_mykind
 zlenou_w =  0.40_mykind
!Inner z scale for u,v,w
 zlenin_u = min(150._mykind/retauinflow,zlenou_u)
 zlenin_v = min( 75._mykind/retauinflow,zlenou_v)
 zlenin_w = min(150._mykind/retauinflow,zlenou_w)
 do j=1,ny
  ftany      = 0.5_mykind*(1._mykind+tanh((y(j)-.2_mykind)/.03_mykind)) ! blending function
  zlen(1,j)  = zlenin_u+ftany*(zlenou_u-zlenin_u)
  zlen(2,j)  = zlenin_v+ftany*(zlenou_v-zlenin_v)
  zlen(3,j)  = zlenin_w+ftany*(zlenou_w-zlenin_w)
 enddo
 ylen = 0.7_mykind*zlen ! Xie & Castro, FTC 2008
 xlen_df(1) = .80_mykind
 xlen_df(2) = .30_mykind
 xlen_df(3) = .30_mykind
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
