subroutine bcdf(ilat)
!
! Digital filtering for turbulent inflow
!
 use mod_streams
 implicit none
!
 integer :: ind,i,j,k,l,m,jj,kk,ll,kb,kg,ilat
 real(mykind) :: avrf,avrf2,dwdxfac,rho,rhofac,rhold,rhomean
 real(mykind) :: rhouu,rhovv,rhoww,rhowall,rmsrf,tfluc,tlen,tmean,tt,vfrms
 real(mykind) :: uu,ufluc,uumean
 real(mykind) :: vv,vfluc,vvmean
 real(mykind) :: ww,wfluc,wwmean
 real(mykind) :: exparg
 real(mykind) :: exp1_1,sqrtexp_1
 real(mykind) :: exp1_2,sqrtexp_2
 real(mykind) :: exp1_3,sqrtexp_3
!
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ny
   do m=1,3
    vf_df_old(m,j,k) = vf_df_gpu(m,j,k)
    vf_df_gpu(m,j,k) = 0._mykind
    uf(m,j,k) = 0._mykind
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 !$cuf kernel do(2) <<<*,*>>>
 do k=1-nfmax,nz+nfmax
  do j=1,ny
   do m=1,3
    rfy_gpu(m,j,k) = 0._mykind
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
! Convolution in y (from rf_gpu to rfy_gpu)
!
 !$cuf kernel do(2) <<<*,*>>>
 do k=1-nfmax,nz+nfmax
  do j=1,ny
   do m=1,3
    do jj=-nfmax,nfmax
     rfy_gpu(m,j,k) = rfy_gpu(m,j,k)+by_df_gpu(m,j,jj)*rf_gpu(m,j+jj,k)
    enddo
   enddo
  enddo
 enddo      
 !@cuf iercuda=cudaDeviceSynchronize()
!
! Convolution in z (from rfy_gpu to vf_df_gpu)
!
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ny
   do m=1,3
    do kk=-nfmax,nfmax
     vf_df_gpu(m,j,k) = vf_df_gpu(m,j,k)+bz_df_gpu(m,j,kk)*rfy_gpu(m,j,k+kk)
    enddo
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
!Provide correlation in time
!
 tlen      = xlen_df(1)/u0
 exparg    = -pi*alpdtold/tlen
 sqrtexp_1 = sqrt(1._mykind-exp(exparg))
 exp1_1    = exp(0.5_mykind*exparg) 
 tlen      = xlen_df(2)/u0
 exparg    = -pi*alpdtold/tlen
 sqrtexp_2 = sqrt(1._mykind-exp(exparg))
 exp1_2    = exp(0.5_mykind*exparg) 
 tlen      = xlen_df(3)/u0
 exparg    = -pi*alpdtold/tlen
 sqrtexp_3 = sqrt(1._mykind-exp(exparg))
 exp1_3    = exp(0.5_mykind*exparg) 
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ny
   vf_df_gpu(1,j,k) = vf_df_old(1,j,k)*exp1_1+vf_df_gpu(1,j,k)*sqrtexp_1
   vf_df_old(1,j,k) = vf_df_gpu(1,j,k)
   vf_df_gpu(2,j,k) = vf_df_old(2,j,k)*exp1_2+vf_df_gpu(2,j,k)*sqrtexp_2
   vf_df_old(2,j,k) = vf_df_gpu(2,j,k)
   vf_df_gpu(3,j,k) = vf_df_old(3,j,k)*exp1_3+vf_df_gpu(3,j,k)*sqrtexp_3
   vf_df_old(3,j,k) = vf_df_gpu(3,j,k)
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
! Remove spanwise mean and normalize
!
!vfav  = 0._mykind
!do k=1,nz
! do j=2,ny
!  do m=1,3
!   vfav(1,m,j) = vfav(1,m,j)+vf_df(m,j,k)
!   vfav(2,m,j) = vfav(2,m,j)+vf_df(m,j,k)**2
!  enddo
! enddo
!enddo
!call mpi_allreduce(vfav,vfavg,6*ny,mpi_prec,mpi_sum,mp_cartz,iermpi)
!vfav = vfavg/nzmax
!do k=1,nz
! do j=2,ny
!  do m=1,3
!   vfrms = sqrt(vfav(2,m,j)-vfav(1,m,j)**2)
!   vf_df_tmp(m,j,k) = (vf_df(m,j,k)-vfav(1,m,j))/vfrms
!  enddo
! enddo
!enddo
!
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=2,ny
   do m=1,3
    do l=1,3
     uf(m,j,k) = uf(m,j,k)+amat_df_gpu(m,l,j)*vf_df_gpu(l,j,k)
    enddo
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 if (ilat==1) then     ! left side
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
    do i=1,ngdf-1
!
     dwdxfac = alpdtold*u0/(xg_gpu(1-i)-xg_gpu(-i))
!
! Interpolate to apply Taylor hypothesis
!
     do m=1,nv
      w_gpu(1-i,j,k,m) = w_gpu(1-i,j,k,m)-dwdxfac*(w_gpu(1-i,j,k,m)-w_gpu(-i,j,k,m))
     enddo        
!
    enddo
    do i=ngdf,ngdf
     rhomean = wmean_gpu(1,1-i,j)
     uumean  = wmean_gpu(2,1-i,j)/rhomean
     vvmean  = wmean_gpu(3,1-i,j)/rhomean
     wwmean  = wmean_gpu(4,1-i,j)/rhomean
     tmean   = p0/rhomean
     ufluc   = uf(1,j,k)*u0
     vfluc   = uf(2,j,k)*u0
     wfluc   = uf(3,j,k)*u0
     tfluc   = -gm1/gamma*ufluc*uumean/tmean
     tfluc   = tfluc*dftscaling
     tt      = tmean*(1._mykind+tfluc)
     rho     = p0/tt
     uu      = uumean +ufluc
     vv      = vvmean +vfluc
     ww      = wwmean +wfluc
     rhouu   = rho*uu
     rhovv   = rho*vv
     rhoww   = rho*ww
     w_gpu(1-i,j,k,1) = rho
     w_gpu(1-i,j,k,2) = rhouu
     w_gpu(1-i,j,k,3) = rhovv
     w_gpu(1-i,j,k,4) = rhoww
     w_gpu(1-i,j,k,5) = p0*gm + 0.5_mykind*(rhouu**2+ rhovv**2+ rhoww**2)/rho
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==2) then  ! right
 elseif (ilat==3) then  ! bottom side
 elseif (ilat==4) then  ! top side
 elseif (ilat==5) then  ! back side
 elseif (ilat==6) then  ! fore side
 endif
! 
end subroutine bcdf
