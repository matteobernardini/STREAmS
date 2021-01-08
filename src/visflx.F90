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
 real(mykind) :: sigq,sigx,sigy,sigz,tt,sigqt,sigah
 real(mykind) :: uu,vv,ww
 real(mykind) :: tx,ty,tz
 real(mykind) :: ux,uy,uz
 real(mykind) :: vx,vy,vz
 real(mykind) :: wx,wy,wz
 real(mykind) :: ulap,ulapx,ulapy,ulapz
 real(mykind) :: vlap,vlapx,vlapy,vlapz
 real(mykind) :: wlap,wlapx,wlapy,wlapz
 real(mykind) :: tlap,tlapx,tlapy,tlapz
 real(mykind) :: sqgmr2,sqgmr2h,tt2,sqrtt,sdivt,sdivt1
 real(mykind) :: dy
!
 real(mykind) :: ttt,st,et
!
 sqgmr2  = sqgmr*(1._mykind+s2tinf)
 sqgmr2h = sqgmr2*0.5_mykind
!
!Update viscous fluxes
!
!st = mpi_wtime()

 !$cuf kernel do(2) <<<*,*>>>
 do j=1,ny
  do i=1,nx
   do k=1,nz
!
    uu = wv_gpu(i,j,k,2)
    vv = wv_gpu(i,j,k,3)
    ww = wv_gpu(i,j,k,4)
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
     ux = ux+ccl*(wv_gpu(i+l,j,k,2)-wv_gpu(i-l,j,k,2))
     vx = vx+ccl*(wv_gpu(i+l,j,k,3)-wv_gpu(i-l,j,k,3))
     wx = wx+ccl*(wv_gpu(i+l,j,k,4)-wv_gpu(i-l,j,k,4))
     tx = tx+ccl*(temperature_gpu(i+l,j,k)-temperature_gpu(i-l,j,k))
!
     uy = uy+ccl*(wv_gpu(i,j+l,k,2)-wv_gpu(i,j-l,k,2))
     vy = vy+ccl*(wv_gpu(i,j+l,k,3)-wv_gpu(i,j-l,k,3))
     wy = wy+ccl*(wv_gpu(i,j+l,k,4)-wv_gpu(i,j-l,k,4))
     ty = ty+ccl*(temperature_gpu(i,j+l,k)-temperature_gpu(i,j-l,k))
!
     uz = uz+ccl*(wv_gpu(i,j,k+l,2)-wv_gpu(i,j,k-l,2))
     vz = vz+ccl*(wv_gpu(i,j,k+l,3)-wv_gpu(i,j,k-l,3))
     wz = wz+ccl*(wv_gpu(i,j,k+l,4)-wv_gpu(i,j,k-l,4))
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
     ulapx = ulapx + clapl*(wv_gpu(i+l,j,k,2)+wv_gpu(i-l,j,k,2))
     ulapy = ulapy + clapl*(wv_gpu(i,j+l,k,2)+wv_gpu(i,j-l,k,2))
     ulapz = ulapz + clapl*(wv_gpu(i,j,k+l,2)+wv_gpu(i,j,k-l,2))
     vlapx = vlapx + clapl*(wv_gpu(i+l,j,k,3)+wv_gpu(i-l,j,k,3))
     vlapy = vlapy + clapl*(wv_gpu(i,j+l,k,3)+wv_gpu(i,j-l,k,3))
     vlapz = vlapz + clapl*(wv_gpu(i,j,k+l,3)+wv_gpu(i,j,k-l,3))
     wlapx = wlapx + clapl*(wv_gpu(i+l,j,k,4)+wv_gpu(i-l,j,k,4))
     wlapy = wlapy + clapl*(wv_gpu(i,j+l,k,4)+wv_gpu(i,j-l,k,4))
     wlapz = wlapz + clapl*(wv_gpu(i,j,k+l,4)+wv_gpu(i,j,k-l,4))
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
    fhat_gpu(i,j,k,6) = div3l
    tt      = temperature_gpu(i,j,k)
    if (visc_type==1) then
     rmut    = sqgmr*tt**vtexp
     drmutdt = rmut*vtexp/tt
    else
     tt2     = tt*tt
     sqrtt   = sqrt(tt)
     sdivt   = s2tinf/tt
     sdivt1  = 1._mykind+sdivt
     rmut    = sqgmr2*sqrtt/sdivt1  ! molecular viscosity
     drmutdt = sqgmr2h/sqrtt+rmut*s2tinf/tt2
     drmutdt = drmutdt/sdivt1
    endif
    rmutx   = drmutdt*tx
    rmuty   = drmutdt*ty
    rmutz   = drmutdt*tz
    sig11 = 2._mykind*(ux-div3l)
    sig12 = uy+vx 
    sig13 = uz+wx
    sig22 = 2._mykind*(vy-div3l)
    sig23 = vz+wy
    sig33 = 2._mykind*(wz-div3l)
    sigx  = rmutx*sig11       +  rmuty*sig12       +  rmutz*sig13       +  rmut*(ulap)
    sigy  = rmutx*sig12       +  rmuty*sig22       +  rmutz*sig23       +  rmut*(vlap)
    sigz  = rmutx*sig13       +  rmuty*sig23       +  rmutz*sig33       +  rmut*(wlap)
    sigqt = (rmutx*tx+rmuty*ty+rmutz*tz+rmut*tlap)*ggmopr ! Conduction
    if (iflow==0) fhat_gpu(i,j,k,1) = sigqt
    sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*rmut ! Aerodynamic heating
    if (iflow==0) fhat_gpu(i,j,k,2) = sigah
!
    sigq = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt
!
    fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
    fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
    fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
    fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq 
!
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 if (iflow==0.and.mod(icyc-ncyc0,nprint)==0) call heatflux_compute()
!
!et = mpi_wtime()
!ttt = et-st
!if (masterproc) write(error_unit,*) 'Viscous-I time =', ttt
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
 real(mykind) :: sqgmr2,sqgmr2h,tt2,sqrtt,sdivt,sdivt1
!
 real(mykind) :: ttt,st,et
!
 sqgmr2  = sqgmr*(1._mykind+s2tinf)
 sqgmr2h = sqgmr2*0.5_mykind
!
!Update viscous fluxes
!
!st = mpi_wtime()
!
 !$cuf kernel do(3) <<<*,*,stream=stream1>>>
 do k=1,nz
  do j=1,ny
   do i=1,nx
!
    uu = wv_gpu(i,j,k,2)
    vv = wv_gpu(i,j,k,3)
    ww = wv_gpu(i,j,k,4)
!
    divx3l  = 0._mykind
    divy3l  = 0._mykind
    divz3l  = 0._mykind
!
    do l=1,ivis/2
     ccl = coeff_deriv1_gpu(l)
     divx3l = divx3l+ccl*(fhat_gpu(i+l,j,k,6)-fhat_gpu(i-l,j,k,6))
     divy3l = divy3l+ccl*(fhat_gpu(i,j+l,k,6)-fhat_gpu(i,j-l,k,6))
     divz3l = divz3l+ccl*(fhat_gpu(i,j,k+l,6)-fhat_gpu(i,j,k-l,6))
    enddo
    divx3l = divx3l*dcsidx_gpu(i)
    divy3l = divy3l*detady_gpu(j)
    divz3l = divz3l*dzitdz_gpu(k)
!
    tt   = temperature_gpu(i,j,k)
    if (visc_type==1) then
     rmut    = sqgmr*tt**vtexp
    else
     tt2     = tt*tt
     sqrtt   = sqrt(tt)
     sdivt   = s2tinf/tt
     sdivt1  = 1._mykind+sdivt
     rmut    = sqgmr2*sqrtt/sdivt1  ! molecular viscosity
    endif
    sigx = rmut*divx3l
    sigy = rmut*divy3l
    sigz = rmut*divz3l
    sigq = sigx*uu+sigy*vv+sigz*ww
!
    fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
    fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
    fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
    fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq 
!
   enddo
  enddo
 enddo

 !!!!@cuf iercuda=cudaDeviceSynchronize()
!et = mpi_wtime()
!ttt = et-st
!if (masterproc) write(error_unit,*) 'Viscous-II time =', ttt
!
end subroutine visflx_div
