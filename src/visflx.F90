subroutine visflx
!
! Evaluation of the viscous fluxes
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,l
 real(mykind) :: ccl,clapl,div3l,drmutdt
 real(mykind) :: ri,rmut,rmutx,rmuty,rmutz
 real(mykind) :: sig11,sig12,sig13
 real(mykind) :: sig21,sig22,sig23
 real(mykind) :: sig31,sig32,sig33
 real(mykind) :: sigq,sigx,sigy,sigz,tt
 real(mykind) :: uu,vv,ww
 real(mykind) :: tx,ty,tz
 real(mykind) :: ux,uy,uz
 real(mykind) :: vx,vy,vz
 real(mykind) :: wx,wy,wz
 real(mykind) :: ulap,ulapx,ulapy,ulapz
 real(mykind) :: vlap,vlapx,vlapy,vlapz
 real(mykind) :: wlap,wlapx,wlapy,wlapz
 real(mykind) :: tlap,tlapx,tlapy,tlapz
!
!sqgmr2  = sqgmr*(1._mykind+s2tinf)
!sqgmr2h = sqgmr2*0.5_mykind
!
!Update viscous fluxes
!
 !$cuf kernel do(2) <<<*,*>>>
 do j=1,ny
  do i=1,nx
   do k=1,nz
!
    uu = w_gpu(2,i,j,k)/w_gpu(1,i,j,k)
    vv = w_gpu(3,i,j,k)/w_gpu(1,i,j,k)
    ww = w_gpu(4,i,j,k)/w_gpu(1,i,j,k)
!
    ux = 0._mykind
    vx = 0._mykind
    wx = 0._mykind
    tx = 0._mykind
    uy = 0._mykind
    vy = 0._mykind
    wy = 0._mykind
    ty = 0._mykind
    uz = 0._mykind
    vz = 0._mykind
    wz = 0._mykind
    tz = 0._mykind
    do l=1,ivis/2
     ccl = coeff_deriv1_gpu(l)
     ux = ux+ccl*(w_gpu(2,i+l,j,k)/w_gpu(1,i+l,j,k)-w_gpu(2,i-l,j,k)/w_gpu(1,i-l,j,k))
     vx = vx+ccl*(w_gpu(3,i+l,j,k)/w_gpu(1,i+l,j,k)-w_gpu(3,i-l,j,k)/w_gpu(1,i-l,j,k))
     wx = wx+ccl*(w_gpu(4,i+l,j,k)/w_gpu(1,i+l,j,k)-w_gpu(4,i-l,j,k)/w_gpu(1,i-l,j,k))
     tx = tx+ccl*(temperature_gpu(i+l,j,k)-temperature_gpu(i-l,j,k))
!
     uy = uy+ccl*(w_gpu(2,i,j+l,k)/w_gpu(1,i,j+l,k)-w_gpu(2,i,j-l,k)/w_gpu(1,i,j-l,k))
     vy = vy+ccl*(w_gpu(3,i,j+l,k)/w_gpu(1,i,j+l,k)-w_gpu(3,i,j-l,k)/w_gpu(1,i,j-l,k))
     wy = wy+ccl*(w_gpu(4,i,j+l,k)/w_gpu(1,i,j+l,k)-w_gpu(4,i,j-l,k)/w_gpu(1,i,j-l,k))
     ty = ty+ccl*(temperature_gpu(i,j+l,k)-temperature_gpu(i,j-l,k))
!
     uz = uz+ccl*(w_gpu(2,i,j,k+l)/w_gpu(1,i,j,k+l)-w_gpu(2,i,j,k-l)/w_gpu(1,i,j,k-l))
     vz = vz+ccl*(w_gpu(3,i,j,k+l)/w_gpu(1,i,j,k+l)-w_gpu(3,i,j,k-l)/w_gpu(1,i,j,k-l))
     wz = wz+ccl*(w_gpu(4,i,j,k+l)/w_gpu(1,i,j,k+l)-w_gpu(4,i,j,k-l)/w_gpu(1,i,j,k-l))
     tz = tz+ccl*(temperature_gpu(i,j,k+l)-temperature_gpu(i,j,k-l))
    enddo
    ux = ux*dcsidx_gpu(i)
    vx = vx*dcsidx_gpu(i)
    wx = wx*dcsidx_gpu(i)
    tx = tx*dcsidx_gpu(i)
    uy = uy*detady_gpu(j)
    vy = vy*detady_gpu(j)
    wy = wy*detady_gpu(j)
    ty = ty*detady_gpu(j)
    uz = uz*dzitdz_gpu(k)
    vz = vz*dzitdz_gpu(k)
    wz = wz*dzitdz_gpu(k)
    tz = tz*dzitdz_gpu(k)
!
    ulapx = coeff_clap_gpu(0)*uu
    ulapy = ulapx
    ulapz = ulapx
    vlapx = coeff_clap_gpu(0)*vv
    vlapy = vlapx
    vlapz = vlapx
    wlapx = coeff_clap_gpu(0)*ww
    wlapy = wlapx
    wlapz = wlapx
    tlapx = coeff_clap_gpu(0)*temperature_gpu(i,j,k)
    tlapy = tlapx
    tlapz = tlapx
    do l=1,ivis/2
     clapl = coeff_clap_gpu(l)
     ulapx = ulapx + clapl*(w_gpu(2,i+l,j,k)/w_gpu(1,i+l,j,k)+w_gpu(2,i-l,j,k)/w_gpu(1,i-l,j,k))
     ulapy = ulapy + clapl*(w_gpu(2,i,j+l,k)/w_gpu(1,i,j+l,k)+w_gpu(2,i,j-l,k)/w_gpu(1,i,j-l,k))
     ulapz = ulapz + clapl*(w_gpu(2,i,j,k+l)/w_gpu(1,i,j,k+l)+w_gpu(2,i,j,k-l)/w_gpu(1,i,j,k-l))
     vlapx = vlapx + clapl*(w_gpu(3,i+l,j,k)/w_gpu(1,i+l,j,k)+w_gpu(3,i-l,j,k)/w_gpu(1,i-l,j,k))
     vlapy = vlapy + clapl*(w_gpu(3,i,j+l,k)/w_gpu(1,i,j+l,k)+w_gpu(3,i,j-l,k)/w_gpu(1,i,j-l,k))
     vlapz = vlapz + clapl*(w_gpu(3,i,j,k+l)/w_gpu(1,i,j,k+l)+w_gpu(3,i,j,k-l)/w_gpu(1,i,j,k-l))
     wlapx = wlapx + clapl*(w_gpu(4,i+l,j,k)/w_gpu(1,i+l,j,k)+w_gpu(4,i-l,j,k)/w_gpu(1,i-l,j,k))
     wlapy = wlapy + clapl*(w_gpu(4,i,j+l,k)/w_gpu(1,i,j+l,k)+w_gpu(4,i,j-l,k)/w_gpu(1,i,j-l,k))
     wlapz = wlapz + clapl*(w_gpu(4,i,j,k+l)/w_gpu(1,i,j,k+l)+w_gpu(4,i,j,k-l)/w_gpu(1,i,j,k-l))
     tlapx = tlapx + clapl*(temperature_gpu(i+l,j,k)+temperature_gpu(i-l,j,k))
     tlapy = tlapy + clapl*(temperature_gpu(i,j+l,k)+temperature_gpu(i,j-l,k))
     tlapz = tlapz + clapl*(temperature_gpu(i,j,k+l)+temperature_gpu(i,j,k-l))
    enddo
!
    ulapx = ulapx*dcsidxs_gpu(i)+ux*dcsidx2_gpu(i)
    vlapx = vlapx*dcsidxs_gpu(i)+vx*dcsidx2_gpu(i)
    wlapx = wlapx*dcsidxs_gpu(i)+wx*dcsidx2_gpu(i)
    tlapx = tlapx*dcsidxs_gpu(i)+tx*dcsidx2_gpu(i)
    ulapy = ulapy*detadys_gpu(j)+uy*detady2_gpu(j)
    vlapy = vlapy*detadys_gpu(j)+vy*detady2_gpu(j)
    wlapy = wlapy*detadys_gpu(j)+wy*detady2_gpu(j)
    tlapy = tlapy*detadys_gpu(j)+ty*detady2_gpu(j)
    ulapz = ulapz*dzitdzs_gpu(k)+uz*dzitdz2_gpu(k)
    vlapz = vlapz*dzitdzs_gpu(k)+vz*dzitdz2_gpu(k)
    wlapz = wlapz*dzitdzs_gpu(k)+wz*dzitdz2_gpu(k)
    tlapz = tlapz*dzitdzs_gpu(k)+tz*dzitdz2_gpu(k)
!
    ulap  = ulapx+ulapy+ulapz
    vlap  = vlapx+vlapy+vlapz
    wlap  = wlapx+wlapy+wlapz
    tlap  = tlapx+tlapy+tlapz
!
    div3l   = ux+vy+wz
    div3l   = div3l/3._mykind
    fhat_gpu(6,i,j,k) = div3l
    tt      = temperature_gpu(i,j,k)
!   tt2     = tt*tt
!   sqrtt   = sqrt(tt)
!   sdivt   = s2tinf/tt
!   sdivt1  = 1._mykind+sdivt
!   rmut    = sqgmr2*sqrtt/sdivt1  ! molecular viscosity
!   drmutdt = sqgmr2h/sqrtt+rmut*s2tinf/tt2
!   drmutdt = drmutdt/sdivt1
    rmut    = sqgmr*tt**vtexp
    drmutdt = rmut*vtexp/tt
    rmutx   = drmutdt*tx
    rmuty   = drmutdt*ty
    rmutz   = drmutdt*tz
    sig11 = 2._mykind*(ux-div3l)
    sig12 = uy+vx 
    sig13 = uz+wx
    sig22 = 2._mykind*(vy-div3l)
    sig23 = vz+wy
    sig33 = 2._mykind*(wz-div3l)
    sigx =  rmutx*sig11       +  rmuty*sig12       +  rmutz*sig13       +  rmut*(ulap)
    sigy =  rmutx*sig12       +  rmuty*sig22       +  rmutz*sig23       +  rmut*(vlap)
    sigz =  rmutx*sig13       +  rmuty*sig23       +  rmutz*sig33       +  rmut*(wlap)
!
    sigq = sigx*uu+sigy*vv+sigz*ww+(sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*rmut+&
          (rmutx*tx+rmuty*ty+rmutz*tz+rmut*tlap)*ggmopr 
!
    fl_gpu(2,i,j,k) = fl_gpu(2,i,j,k) - sigx
    fl_gpu(3,i,j,k) = fl_gpu(3,i,j,k) - sigy
    fl_gpu(4,i,j,k) = fl_gpu(4,i,j,k) - sigz
    fl_gpu(5,i,j,k) = fl_gpu(5,i,j,k) - sigq 
!
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
end subroutine visflx
!
subroutine visflx_div
!
! Evaluation of the viscous fluxes (add terms with flow divergence)
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,l
 real(mykind) :: ccl,clapl,div3l,divx3l,divy3l,divz3l
 real(mykind) :: ri,rmut
 real(mykind) :: sigq,sigx,sigy,sigz,tt
 real(mykind) :: uu,vv,ww
!
!sqgmr2  = sqgmr*(1._mykind+s2tinf)
!sqgmr2h = sqgmr2*0.5_mykind
!
!Update viscous fluxes
!
 !$cuf kernel do(3) <<<*,*,stream=stream1>>>
 do k=1,nz
  do j=1,ny
   do i=1,nx
!
    uu = w_gpu(2,i,j,k)/w_gpu(1,i,j,k)
    vv = w_gpu(3,i,j,k)/w_gpu(1,i,j,k)
    ww = w_gpu(4,i,j,k)/w_gpu(1,i,j,k)
!
    divx3l  = 0._mykind
    divy3l  = 0._mykind
    divz3l  = 0._mykind
!
    do l=1,ivis/2
     ccl = coeff_deriv1_gpu(l)
     divx3l = divx3l+ccl*(fhat_gpu(6,i+l,j,k)-fhat_gpu(6,i-l,j,k))
     divy3l = divy3l+ccl*(fhat_gpu(6,i,j+l,k)-fhat_gpu(6,i,j-l,k))
     divz3l = divz3l+ccl*(fhat_gpu(6,i,j,k+l)-fhat_gpu(6,i,j,k-l))
    enddo
    divx3l = divx3l*dcsidx_gpu(i)
    divy3l = divy3l*detady_gpu(j)
    divz3l = divz3l*dzitdz_gpu(k)
!
    tt   = temperature_gpu(i,j,k)
    rmut = sqgmr*tt**vtexp
    sigx = rmut*divx3l
    sigy = rmut*divy3l
    sigz = rmut*divz3l
    sigq = sigx*uu+sigy*vv+sigz*ww
!
    fl_gpu(2,i,j,k) = fl_gpu(2,i,j,k) - sigx
    fl_gpu(3,i,j,k) = fl_gpu(3,i,j,k) - sigy
    fl_gpu(4,i,j,k) = fl_gpu(4,i,j,k) - sigz
    fl_gpu(5,i,j,k) = fl_gpu(5,i,j,k) - sigq 
!
   enddo
  enddo
 enddo
!
end subroutine visflx_div
