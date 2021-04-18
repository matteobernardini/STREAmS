subroutine visflx_stag
!
! Evaluation of the viscous fluxes
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,l,iv
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
 real(mykind) :: sqgmr2,sqgmr2h,tt2,sqrtt,sdivt,sdivt1
 real(mykind) :: rmuf
 real(mykind) :: sigqq,sigxx,sigyy,sigzz
 real(mykind) :: dxhl,dyhl,dzhl
 real(mykind) :: fl2,fl3,fl4,fl5
 real(mykind) :: fl2o,fl3o,fl4o,fl5o
 real(mykind) :: cleft,cright
!
 real(mykind) :: ttt,st,et
!
! st = mpi_wtime()
!
 sqgmr2  = sqgmr*(1._mykind+s2tinf)
 sqgmr2h = sqgmr2*0.5_mykind
!
! ddx
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ny
   do i=0,nx
    tt   = 0._mykind
!   rmuf = 0._mykind
    sigx = 0._mykind
    sigy = 0._mykind
    sigz = 0._mykind
    sigq = 0._mykind
    do l=1,ivis/2
     ccl  = coeff_deriv1s_gpu(l)
     sigx = sigx+ccl*(wv_trans_gpu(j,i+l,k,2)-wv_trans_gpu(j,i-l+1,k,2))
     sigy = sigy+ccl*(wv_trans_gpu(j,i+l,k,3)-wv_trans_gpu(j,i-l+1,k,3))
     sigz = sigz+ccl*(wv_trans_gpu(j,i+l,k,4)-wv_trans_gpu(j,i-l+1,k,4))
     sigq = sigq+ccl*(temperature_trans_gpu(j,i+l,k)-temperature_trans_gpu(j,i-l+1,k))
     cleft  = cx_midpi_gpu(i,1-l)
     cright = cx_midpi_gpu(i,l)
!    rmuf = rmuf+ccl*(wv_gpu(i+l,j,k,1)+wv_gpu(i-l+1,j,k,1))
     tt   = tt+cright*temperature_trans_gpu(j,i+l,k)+cleft*temperature_trans_gpu(j,i-l+1,k)
    enddo
    if (visc_type==1) then
     rmuf    = sqgmr*tt**vtexp
    else
     tt2     = tt*tt
     sqrtt   = sqrt(tt)
     sdivt   = s2tinf/tt
     sdivt1  = 1._mykind+sdivt
     rmuf    = sqgmr2*sqrtt/sdivt1  ! molecular viscosity
    endif
!
    rmuf = rmuf*dcsidxh_gpu(i)
    sigx = sigx*rmuf
    sigy = sigy*rmuf
    sigz = sigz*rmuf
    sigq = sigq*rmuf*ggmopr
!
    if (i>0) then
     dxhl  = 1._mykind/(xh_gpu(i)-xh_gpu(i-1))
!    uu    = wv_trans_gpu(j,i,k,2)
!    vv    = wv_trans_gpu(j,i,k,3)
!    ww    = wv_trans_gpu(j,i,k,4)
     sigxx = sigx*dxhl
     sigyy = sigy*dxhl
     sigzz = sigz*dxhl
     sigqq = sigq*dxhl
     sigqq = sigqq+uu*sigxx+vv*sigyy+ww*sigzz
!    fl_trans_gpu(j,i,k,2) = fl_trans_gpu(j,i,k,2) - sigxx
!    fl_trans_gpu(j,i,k,3) = fl_trans_gpu(j,i,k,3) - sigyy
!    fl_trans_gpu(j,i,k,4) = fl_trans_gpu(j,i,k,4) - sigzz
!    fl_trans_gpu(j,i,k,5) = fl_trans_gpu(j,i,k,5) - sigqq
     fl2 = fl2o-sigxx
     fl3 = fl3o-sigyy
     fl4 = fl4o-sigzz
     fl5 = fl5o-sigqq
    endif
    if (i<nx) then
     dxhl  = 1._mykind/(xh_gpu(i+1)-xh_gpu(i))
     uu    = wv_trans_gpu(j,i+1,k,2)
     vv    = wv_trans_gpu(j,i+1,k,3)
     ww    = wv_trans_gpu(j,i+1,k,4)
     sigxx = sigx*dxhl
     sigyy = sigy*dxhl
     sigzz = sigz*dxhl
     sigqq = sigq*dxhl
     sigqq = sigqq+uu*sigxx+vv*sigyy+ww*sigzz
!    fl_trans_gpu(j,i+1,k,2) = fl_trans_gpu(j,i+1,k,2) + sigxx
!    fl_trans_gpu(j,i+1,k,3) = fl_trans_gpu(j,i+1,k,3) + sigyy
!    fl_trans_gpu(j,i+1,k,4) = fl_trans_gpu(j,i+1,k,4) + sigzz
!    fl_trans_gpu(j,i+1,k,5) = fl_trans_gpu(j,i+1,k,5) + sigqq
     fl2o = sigxx 
     fl3o = sigyy 
     fl4o = sigzz 
     fl5o = sigqq
    endif
    if (i>0) then
     fl_trans_gpu(j,i,k,2) = fl2
     fl_trans_gpu(j,i,k,3) = fl3
     fl_trans_gpu(j,i,k,4) = fl4
     fl_trans_gpu(j,i,k,5) = fl5
    endif
!
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 !$cuf kernel do(3) <<<*,*,stream=stream1>>>
  do k=1,nz
   do j=1,ny
    do i=1,nx
     do iv=2,nv
      fl_gpu(i,j,k,iv) = fl_gpu(i,j,k,iv)+fl_trans_gpu(j,i,k,iv)
     enddo
    enddo
   enddo
  enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
! ddy
 !$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do i=1,nx
   do j=0,ny
    tt   = 0._mykind
!   rmuf = 0._mykind
    sigx = 0._mykind
    sigy = 0._mykind
    sigz = 0._mykind
    sigq = 0._mykind
    do l=1,ivis/2
     ccl  = coeff_deriv1s_gpu(l)
     sigx = sigx+ccl*(wv_gpu(i,j+l,k,2)-wv_gpu(i,j-l+1,k,2))
     sigy = sigy+ccl*(wv_gpu(i,j+l,k,3)-wv_gpu(i,j-l+1,k,3))
     sigz = sigz+ccl*(wv_gpu(i,j+l,k,4)-wv_gpu(i,j-l+1,k,4))
     sigq = sigq+ccl*(temperature_gpu(i,j+l,k)-temperature_gpu(i,j-l+1,k))
     cleft  = cy_midpi_gpu(j,1-l)
     cright = cy_midpi_gpu(j,l)
!    rmuf = rmuf+ccl*(wv_gpu(i,j+l,k,1)+wv_gpu(i,j-l+1,k,1))
     tt   = tt+cright*temperature_gpu(i,j+l,k)+cleft*temperature_gpu(i,j-l+1,k)
    enddo
    if (visc_type==1) then
     rmuf    = sqgmr*tt**vtexp
    else
     tt2     = tt*tt
     sqrtt   = sqrt(tt)
     sdivt   = s2tinf/tt
     sdivt1  = 1._mykind+sdivt
     rmuf    = sqgmr2*sqrtt/sdivt1  ! molecular viscosity
    endif
!
    rmuf = rmuf*detadyh_gpu(j)
    sigx = sigx*rmuf
    sigy = sigy*rmuf
    sigz = sigz*rmuf
    sigq = sigq*rmuf*ggmopr
!
    if (j>0) then ! left
     dyhl  = 1._mykind/(yh_gpu(j)-yh_gpu(j-1))
!    uu    = wv_gpu(i,j,k,2)
!    vv    = wv_gpu(i,j,k,3)
!    ww    = wv_gpu(i,j,k,4)
     sigxx = sigx*dyhl
     sigyy = sigy*dyhl
     sigzz = sigz*dyhl
     sigqq = sigq*dyhl
     sigqq = sigqq+uu*sigxx+vv*sigyy+ww*sigzz
!    fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigxx
!    fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigyy
!    fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigzz
!    fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigqq
     fl2 = fl2o-sigxx
     fl3 = fl3o-sigyy
     fl4 = fl4o-sigzz
     fl5 = fl5o-sigqq
    endif
    if (j<ny) then ! right
     dyhl  = 1._mykind/(yh_gpu(j+1)-yh_gpu(j))
     uu   = wv_gpu(i,j+1,k,2)
     vv   = wv_gpu(i,j+1,k,3)
     ww   = wv_gpu(i,j+1,k,4)
     sigxx = sigx*dyhl
     sigyy = sigy*dyhl
     sigzz = sigz*dyhl
     sigqq = sigq*dyhl
     sigqq = sigqq+uu*sigxx+vv*sigyy+ww*sigzz
!    fl_gpu(i,j+1,k,2) = fl_gpu(i,j+1,k,2) + sigxx
!    fl_gpu(i,j+1,k,3) = fl_gpu(i,j+1,k,3) + sigyy
!    fl_gpu(i,j+1,k,4) = fl_gpu(i,j+1,k,4) + sigzz
!    fl_gpu(i,j+1,k,5) = fl_gpu(i,j+1,k,5) + sigqq
     fl2o = sigxx 
     fl3o = sigyy 
     fl4o = sigzz 
     fl5o = sigqq
    endif
    if (j>0) then
     fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) + fl2
     fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) + fl3
     fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) + fl4
     fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) + fl5
    endif
!
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
! ddz
 !$cuf kernel do(2) <<<*,*>>>
 do j=1,ny
  do i=1,nx
   do k=0,nz
    tt   = 0._mykind
!   rmuf = 0._mykind
    sigx = 0._mykind
    sigy = 0._mykind
    sigz = 0._mykind
    sigq = 0._mykind
    do l=1,ivis/2
     ccl  = coeff_deriv1s_gpu(l)
     sigx = sigx+ccl*(wv_gpu(i,j,k+l,2)-wv_gpu(i,j,k-l+1,2))
     sigy = sigy+ccl*(wv_gpu(i,j,k+l,3)-wv_gpu(i,j,k-l+1,3))
     sigz = sigz+ccl*(wv_gpu(i,j,k+l,4)-wv_gpu(i,j,k-l+1,4))
     sigq = sigq+ccl*(temperature_gpu(i,j,k+l)-temperature_gpu(i,j,k-l+1))
     cleft  = cz_midpi_gpu(k,1-l)
     cright = cz_midpi_gpu(k,l)
!    rmuf = rmuf+ccl*(wv_gpu(i,j,k+l,1)+wv_gpu(i,j,k-l+1,1))
     tt   = tt+cright*temperature_gpu(i,j,k+l)+cleft*temperature_gpu(i,j,k-l+1)
    enddo
    if (visc_type==1) then
     rmuf    = sqgmr*tt**vtexp
    else
     tt2     = tt*tt
     sqrtt   = sqrt(tt)
     sdivt   = s2tinf/tt
     sdivt1  = 1._mykind+sdivt
     rmuf    = sqgmr2*sqrtt/sdivt1  ! molecular viscosity
    endif
!
    rmuf = rmuf*dzitdzh_gpu(k)
    sigx = sigx*rmuf
    sigy = sigy*rmuf
    sigz = sigz*rmuf
    sigq = sigq*rmuf*ggmopr
!
    if (k>0) then
     dzhl  = 1._mykind/(zh_gpu(k)-zh_gpu(k-1))
!    uu    = wv_gpu(i,j,k,2)
!    vv    = wv_gpu(i,j,k,3)
!    ww    = wv_gpu(i,j,k,4)
     sigxx = sigx*dzhl
     sigyy = sigy*dzhl
     sigzz = sigz*dzhl
     sigqq = sigq*dzhl
     sigqq = sigqq+uu*sigxx+vv*sigyy+ww*sigzz
!    fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigxx
!    fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigyy
!    fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigzz
!    fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigqq
     fl2 = fl2o-sigxx
     fl3 = fl3o-sigyy
     fl4 = fl4o-sigzz
     fl5 = fl5o-sigqq
    endif
    if (k<nz) then
     dzhl  = 1._mykind/(zh_gpu(k+1)-zh_gpu(k))
     uu   = wv_gpu(i,j,k+1,2)
     vv   = wv_gpu(i,j,k+1,3)
     ww   = wv_gpu(i,j,k+1,4)
     sigxx = sigx*dzhl
     sigyy = sigy*dzhl
     sigzz = sigz*dzhl
     sigqq = sigq*dzhl
     sigqq = sigqq+uu*sigxx+vv*sigyy+ww*sigzz
!    fl_gpu(i,j,k+1,2) = fl_gpu(i,j,k+1,2) + sigxx
!    fl_gpu(i,j,k+1,3) = fl_gpu(i,j,k+1,3) + sigyy
!    fl_gpu(i,j,k+1,4) = fl_gpu(i,j,k+1,4) + sigzz
!    fl_gpu(i,j,k+1,5) = fl_gpu(i,j,k+1,5) + sigqq
     fl2o = sigxx 
     fl3o = sigyy 
     fl4o = sigzz 
     fl5o = sigqq
    endif
    if (k>0) then
     fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) + fl2
     fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) + fl3
     fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) + fl4
     fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) + fl5
    endif
!
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
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
!
    do l=1,ivis/2
!
     ccl = coeff_deriv1_gpu(l)
!
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
!
    enddo
!
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
    div3l   = ux+vy+wz
    div3l   = div3l/3._mykind
    fhat_gpu(i,j,k,6) = div3l
    tt      = temperature_gpu(i,j,k)
    if (visc_type==1) then
     rmut    = sqgmr*tt**vtexp
     drmutdt = rmut*vtexp/tt
    else
     sqgmr2  = sqgmr*(1._mykind+s2tinf)
     sqgmr2h = sqgmr2*0.5_mykind
     tt2     = tt*tt
     sqrtt   = sqrt(tt)
     sdivt   = s2tinf/tt
     sdivt1  = 1._mykind+sdivt
     rmut    = sqgmr2*sqrtt/sdivt1  ! molecular viscosity
     drmutdt = sqgmr2h/sqrtt+rmut*s2tinf/tt2
     drmutdt = drmutdt/sdivt1
    endif
    rmutx = drmutdt*tx
    rmuty = drmutdt*ty
    rmutz = drmutdt*tz
    sig11 = ux-2._mykind*div3l
    sig12 = vx
    sig13 = wx
    sig21 = uy
    sig22 = vy-2._mykind*div3l
    sig23 = wy
    sig31 = uz
    sig32 = vz
    sig33 = wz-2._mykind*div3l
    sigx  = rmutx*sig11 + rmuty*sig12 + rmutz*sig13
    sigy  = rmutx*sig21 + rmuty*sig22 + rmutz*sig23
    sigz  = rmutx*sig31 + rmuty*sig32 + rmutz*sig33
    sigq  = sigx*uu+sigy*vv+sigz*ww
!
    fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) - sigx
    fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) - sigy
    fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) - sigz
    fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq 
!
    sig11 = 2._mykind*(ux-div3l)
    sig12 = uy+vx
    sig13 = uz+wx
    sig22 = 2._mykind*(vy-div3l)
    sig23 = vz+wy
    sig33 = 2._mykind*(wz-div3l)
    sigq  = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*rmut
!
    fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) - sigq 
!
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()

!et = mpi_wtime()
!ttt = et-st
!if (masterproc) write(error_unit,*) 'Viscous-I time =', ttt
!
end subroutine visflx_stag
