subroutine initurb
!
! BL and SBLI initialization (wmean + perturbations from digital filtering)
!
 use mod_streams
 implicit none
!
 integer :: nx_1,nx_3,ny_1,ny_3,nz_1,nz_3
 real(mykind), dimension(:,:,:,:), allocatable :: ranf,ranf_x,ranf_y,ranf_z
 real(mykind), dimension(:,:,:,:), allocatable :: ranf_y_trasp,ranf_x_trasp,ranf_yz,ranf_xz,ranf_xz_extended
 real(mykind), dimension(:,:,:,:), allocatable :: ranf_av
 real(mykind), dimension(:,:,:,:), allocatable :: velfluc
 real(mykind) :: ranf_rms
!
 integer :: i,j,k,m,l,ii,iii,jj,kk,kkk
 real(mykind) :: rho,rhofac,rhomean,rhouu,rhovv,rhoww,rhowall
 real(mykind) :: tfluc,tmean,tt
 real(mykind) :: ufluc,uu,uumean
 real(mykind) :: vfluc,vv,vvmean
 real(mykind) :: wfluc,ww,wwmean
!
 nx_1 = nxmax/nblocks(1)
 nx_3 = nxmax/nblocks(3)
 ny_1 = nymax/nblocks(1)
 ny_3 = nymax/nblocks(3)
 nz_1 = nzmax/nblocks(1)
 nz_3 = nzmax/nblocks(3)
!
 allocate(ranf(3,1:nx_1,1-nfmax:nymax+nfmax,1:nz_3))
!
!Random field generation
!
 do k = 1,nz_3
  do j = 1-nfmax,nymax+nfmax
   do i = 1,nx_1
    do m = 1,3
     call gasdev_s(ranf(m,i,j,k))
    enddo
   enddo
  enddo
 enddo
!
! Convolution in y
!
 allocate(ranf_y(3,1:nx_1,1:nymax,1:nz_3))
 ranf_y = 0._mykind
 do k=1,nz_3
  do j=1,nymax
   do i=1,nx_1
    do jj=-nfmax,nfmax
     do m=1,3
      ranf_y(m,i,j,k) = ranf_y(m,i,j,k)+by_df(m,j,jj)*ranf(m,i,j+jj,k)
     enddo
    enddo
   enddo
  enddo
 enddo
!
 deallocate(ranf)
 allocate(ranf_y_trasp(3,1:nxmax,1:ny_1,1:nz_3))
!
! Perform transposition ranf_y ==> ranf_y_trasp
!
 call trasp_xz_to_yz(ranf_y,ranf_y_trasp,1,nx_1,ny_1,nz_3,nblocks(1))
!
!Convolution in x
!
 deallocate(ranf_y)
 allocate(ranf_x(3,1:nxmax,1:ny_1,1:nz_3))
 ranf_x = 0._mykind
 do k=1,nz_3
  do j=1,ny_1
   do i=1,nxmax
    do ii=-nfmax,nfmax
     iii = i+ii
     if (iii<1) iii = nxmax+iii
     if (iii>nxmax) iii = iii-nxmax
     do m=1,3
      ranf_x(m,i,j,k) = ranf_x(m,i,j,k)+bx_df(m,i,ii)*ranf_y_trasp(m,iii,j,k)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(ranf_y_trasp)
 allocate(ranf_x_trasp(3,1:nx_3,1:ny_1,1:nzmax))
!
! Perform transposition ranf_x ==> ranf_x_trasp
!
 call trasp_yz_to_xy(ranf_x,ranf_x_trasp,1,nx_3,ny_1,nz_3,nblocks(3))
!
! Convolution in z
!
 deallocate(ranf_x)
 allocate(ranf_z(3,1:nx_3,1:ny_1,1:nzmax))
 ranf_z  = 0._mykind
 do k=1,nzmax
  do j=1,ny_1
   do i=1,nx_3
    do kk=-nfmax,nfmax
     kkk = k+kk
     if (kkk<1) kkk = nzmax+kkk
     if (kkk>nzmax) kkk = kkk-nzmax
     do m=1,3
      ranf_z(m,i,j,k) = ranf_z(m,i,j,k)+bz_df(m,j,kk)*ranf_x_trasp(m,i,j,kkk)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(ranf_x_trasp)
! 
! Remove mean and transpose back
!
 allocate(ranf_av (2,3,1:nx_3,1:ny_1))
 ranf_av = 0._mykind
 do k=1,nzmax
  do j=1,ny_1
   do i=1,nx_3
    do m=1,3
     ranf_av(1,m,i,j) = ranf_av(1,m,i,j)+ranf_z(m,i,j,k)
     ranf_av(2,m,i,j) = ranf_av(2,m,i,j)+ranf_z(m,i,j,k)**2
    enddo
   enddo
  enddo
 enddo
 ranf_av = ranf_av/nzmax
  do k=1,nzmax
  do j=1,ny_1
   do i=1,nx_3
    do m=1,3
     ranf_rms = sqrt(ranf_av(2,m,i,j)-ranf_av(1,m,i,j)**2)
     ranf_z(m,i,j,k) = (ranf_z(m,i,j,k)-ranf_av(1,m,i,j))/ranf_rms
    enddo
   enddo
  enddo
 enddo
!
 allocate(ranf_yz(3,1:nxmax,1:ny_1,1:nz_3))
 call trasp_yz_to_xy(ranf_yz,ranf_z,-1,nx_3,ny_1,nz_3,nblocks(3))
 deallocate(ranf_z)
 allocate(ranf_xz(3,1:nx_1,1:nymax,1:nz_3))
 call trasp_xz_to_yz(ranf_xz,ranf_yz,-1,nx_1,ny_1,nz_3,nblocks(1))
 deallocate(ranf_yz)
 allocate(ranf_xz_extended(3,1-ngdf:nx_1+ngdf,1:nymax,1:nz_3))
 do k=1,nz_3
  do j=1,nymax
   do i=1,nx_1
    do m=1,3
     ranf_xz_extended(m,i,j,k) = ranf_xz(m,i,j,k)
    enddo
   enddo
  enddo
 enddo
 deallocate(ranf_xz)
!
! Swap
!
 call swap_ranfxz(ranf_xz_extended)
!
 vf_df = ranf_xz_extended(:,1-ngdf,:,:) ! First plane
!
 allocate(velfluc(3,1-ngdf:nx+ngdf,ny,nz))
 velfluc = 0._mykind
 do k=1,nz
  do j=2,ny
   do i=1-ngdf,nx+ngdf
    do m=1,3
     do l=1,3
      velfluc(m,i,j,k) = velfluc(m,i,j,k)+amat_df(m,l,j)*ranf_xz_extended(l,i,j,k)
     enddo
    enddo
   enddo
  enddo
 enddo
!
 do k=1,nz
  do j=1,ny
   do i=1-ngdf,nx+ngdf
!
    rhomean = wmean(1,i,j)
    uumean  = wmean(2,i,j)/rhomean
    vvmean  = wmean(3,i,j)/rhomean
    wwmean  = wmean(4,i,j)/rhomean
    tmean   = p0/rhomean
    ufluc   = velfluc(1,i,j,k)*u0
    vfluc   = velfluc(2,i,j,k)*u0
    wfluc   = velfluc(3,i,j,k)*u0
    tfluc   = -gm1/gamma*ufluc*uumean/tmean
!   rmloc   = uumean/sqrt(gamma*tmean)
!   tfluc   = (-gm1*rmloc**2)*ufluc/uumean ! equal to the previous but divergent at j==1
    tfluc   = tfluc*dftscaling
    tt      = tmean*(1._mykind+tfluc)
    rho     = p0/tt
    uu      = uumean+ufluc
    vv      = vvmean+vfluc
    ww      = wwmean+wfluc
    rhouu   = rho*uu
    rhovv   = rho*vv
    rhoww   = rho*ww
    w(1,i,j,k) = rho
    w(2,i,j,k) = rhouu
    w(3,i,j,k) = rhovv
    w(4,i,j,k) = rhoww
    w(5,i,j,k) = p0*gm+0.5_mykind*(rhouu**2+rhovv**2+rhoww**2)/rho
   enddo
  enddo
 enddo
!
end subroutine initurb
