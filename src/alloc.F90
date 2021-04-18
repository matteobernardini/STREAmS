subroutine allocate_vars()
!
! Allocate variables for the computation
!
 use mod_streams
 implicit none
!
#ifdef USE_CUDA
 allocate(fl_trans_gpu(1:ny,1:nx,1:nz,nv))
 allocate(temperature_trans_gpu(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng))
 allocate(fhat_trans_gpu(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,5))

 allocate(wv_trans_gpu(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,nv))

 allocate(w_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,nv))
 allocate(fl_gpu(nx,ny,nz,nv))
 allocate(fln_gpu(nx,ny,nz,nv))
 allocate(temperature_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
 allocate(ducros_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
 allocate(wmean_gpu(nv,1-ng:nx+1+ng,ny))
 allocate(dcsidx_gpu(nx),dcsidx2_gpu(nx),dcsidxs_gpu(nx))
 allocate(detady_gpu(ny),detady2_gpu(ny),detadys_gpu(ny))
 allocate(dzitdz_gpu(nz),dzitdz2_gpu(nz),dzitdzs_gpu(nz))
 allocate(dcsidxh_gpu(0:nx))
 allocate(detadyh_gpu(0:ny))
 allocate(dzitdzh_gpu(0:nz))
 allocate(x_gpu(1-ng:nx+ng))
 allocate(y_gpu(1-ng:ny+ng))
 allocate(z_gpu(1-ng:nz+ng))
 allocate(xg_gpu(1-ng:nxmax+ng+1))
 allocate(coeff_deriv1_gpu(3))
 allocate(coeff_deriv1s_gpu(3))
 allocate(coeff_clap_gpu(0:3))
 allocate(coeff_midpi_gpu(3))
 allocate(cx_midpi_gpu(0:nx,-2:3))
 allocate(cy_midpi_gpu(0:ny,-2:3))
 allocate(cz_midpi_gpu(0:nz,-2:3))
 allocate(fhat_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,6))
 allocate(ibcnr_gpu(6))
 allocate(dcoe_gpu(4,4))
 allocate(winf_gpu(nv),winf1_gpu(nv))
 allocate(rf_gpu(3,1-nfmax:ny+nfmax,1-nfmax:nz+nfmax))
 allocate(rfy_gpu(3,ny,1-nfmax:nz+nfmax))
 allocate(vf_df_gpu(3,ny,nz))
 allocate(by_df_gpu(3,ny,-nfmax:nfmax))
 allocate(bz_df_gpu(3,ny,-nfmax:nfmax))
 allocate(amat_df_gpu(3,3,ny))
!
 allocate(yplus_inflow_gpu(1-ng:ny+ng))
 allocate(yplus_recyc_gpu(1-ng:ny+ng))
 allocate(eta_inflow_gpu(1-ng:ny+ng))
 allocate(eta_recyc_gpu(1-ng:ny+ng))
 allocate(map_j_inn_gpu(ny))
 allocate(map_j_out_gpu(ny))
 allocate(weta_inflow_gpu(ny))
!
 allocate(gplus_x (5,2*iweno,ny,nz))
 allocate(gminus_x(5,2*iweno,ny,nz))
 allocate(gplus_y (5,2*iweno,nx,nz))
 allocate(gminus_y(5,2*iweno,nx,nz))
 allocate(gplus_z (5,2*iweno,nx,ny))
 allocate(gminus_z(5,2*iweno,nx,ny))
#endif
!
 allocate(wv_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,nv))
 allocate(w_order(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,nv))

 allocate(wallpfield_gpu(nx,nz))
 allocate(slicexy_gpu(nv,nx,ny))
 allocate(vf_df_old(3,ny,nz) )
 allocate(uf(3,ny,nz) )
 allocate(evmax_mat_yz(ny,nz) )
 allocate(evmax_mat_y(ny) )
 allocate(bulk5g_gpu(5))
 allocate(rtrms_ib_gpu(ny,nz))
 allocate(rtrms_ib_1d_gpu(ny))
!
 allocate(wbuf1s_gpu(ng,ny,nz,nv))  
 allocate(wbuf2s_gpu(ng,ny,nz,nv))  
 allocate(wbuf3s_gpu(nx,ng,nz,nv))  
 allocate(wbuf4s_gpu(nx,ng,nz,nv))  
 allocate(wbuf5s_gpu(nx,ny,ng,nv))  
 allocate(wbuf6s_gpu(nx,ny,ng,nv))  
 allocate(wbuf1r_gpu(ng,ny,nz,nv))  
 allocate(wbuf2r_gpu(ng,ny,nz,nv))  
 allocate(wbuf3r_gpu(nx,ng,nz,nv))  
 allocate(wbuf4r_gpu(nx,ng,nz,nv))  
 allocate(wbuf5r_gpu(nx,ny,ng,nv))  
 allocate(wbuf6r_gpu(nx,ny,ng,nv))  
 allocate(divbuf1s_gpu(ng,ny,nz))  
 allocate(divbuf2s_gpu(ng,ny,nz))  
 allocate(divbuf3s_gpu(nx,ng,nz))  
 allocate(divbuf4s_gpu(nx,ng,nz))  
 allocate(divbuf5s_gpu(nx,ny,ng))  
 allocate(divbuf6s_gpu(nx,ny,ng))  
 allocate(divbuf1r_gpu(ng,ny,nz))  
 allocate(divbuf2r_gpu(ng,ny,nz))  
 allocate(divbuf3r_gpu(nx,ng,nz))  
 allocate(divbuf4r_gpu(nx,ng,nz))  
 allocate(divbuf5r_gpu(nx,ny,ng))  
 allocate(divbuf6r_gpu(nx,ny,ng))  
 allocate(ducbuf1s_gpu(ng,ny,nz))
 allocate(ducbuf2s_gpu(ng,ny,nz))
 allocate(ducbuf3s_gpu(nx,ng,nz))
 allocate(ducbuf4s_gpu(nx,ng,nz))
 allocate(ducbuf5s_gpu(nx,ny,ng))
 allocate(ducbuf6s_gpu(nx,ny,ng))
 allocate(ducbuf1r_gpu(ng,ny,nz))
 allocate(ducbuf2r_gpu(ng,ny,nz))
 allocate(ducbuf3r_gpu(nx,ng,nz))
 allocate(ducbuf4r_gpu(nx,ng,nz))
 allocate(ducbuf5r_gpu(nx,ny,ng))
 allocate(ducbuf6r_gpu(nx,ny,ng))
!
 allocate(w(nv,1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
 allocate(fl(nx,ny,nz,nv))
 allocate(fln(nx,ny,nz,nv))
 allocate(temperature(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
 allocate(ducros(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
 allocate(wmean(nv,1-ng:nx+1+ng,ny))
 allocate(dcsidx(nx),dcsidx2(nx),dcsidxs(nx))
 allocate(detady(ny),detady2(ny),detadys(ny))
 allocate(dzitdz(nz),dzitdz2(nz),dzitdzs(nz))
 allocate(dcsidxh(0:nx))
 allocate(detadyh(0:ny))
 allocate(dzitdzh(0:nz))
 allocate(x(1-ng:nx+ng))
 allocate(y(1-ng:ny+ng))
 allocate(yn(1:ny+1))    ! Both yn and yn_gpu must be always allocated
 allocate(yn_gpu(1:ny+1))
 allocate(z(1-ng:nz+ng))
 allocate(xg(1-ng:nxmax+ng+1))
 allocate(coeff_deriv1(3))
 allocate(coeff_deriv1s(3))
 allocate(coeff_clap(0:3))
 allocate(coeff_midpi(3))
 allocate(cx_midpi(0:nx,-2:3))
 allocate(cy_midpi(0:ny,-2:3))
 allocate(cz_midpi(0:nz,-2:3))
 allocate(fhat(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,6)) ! Size increased to exchange divergence
 allocate(ibcnr(6))
 allocate(dcoe(4,4))
 allocate(winf(nv),winf1(nv))
 allocate(rf(3,1-nfmax:ny+nfmax,1-nfmax:nz+nfmax))
 allocate(rfy(3,ny,1-nfmax:nz+nfmax))
 allocate(vf_df(3,ny,nz))
 allocate(by_df(3,ny,-nfmax:nfmax))
 allocate(bz_df(3,ny,-nfmax:nfmax))
 allocate(amat_df(3,3,ny))
 allocate(wallpfield(nx,nz))
 allocate(slicexy(nv,nx,ny))
 allocate(xh(0:nx))
 allocate(yh(0:ny))
 allocate(zh(0:nz))
 allocate(xgh(0:nxmax))
 allocate(ygh(0:nymax))
 allocate(zgh(0:nzmax))
 allocate(yplus_inflow(1-ng:ny+ng))
 allocate(yplus_recyc(1-ng:ny+ng))
 allocate(eta_inflow(1-ng:ny+ng))
 allocate(eta_recyc(1-ng:ny+ng))
 allocate(map_j_inn(ny))
 allocate(map_j_out(ny))
 allocate(weta_inflow(ny))
!
 allocate(ibc(6))
 allocate(dxg(nxmax))
 allocate(dyg(nymax))
 allocate(dzg(nzmax))
 allocate(w_av(nvmean,nx,ny))
 allocate(w_avzg(nvmean,nx,ny))
 allocate(w_av_1d(nvmean,ny))
 allocate(w_avxzg(nvmean,ny))
 allocate(bx_df(3,nxmax,-nfmax:nfmax))
 allocate(wbuf1s(ng,ny,nz,nv))  
 allocate(wbuf2s(ng,ny,nz,nv))  
 allocate(wbuf3s(nx,ng,nz,nv))  
 allocate(wbuf4s(nx,ng,nz,nv))  
 allocate(wbuf5s(nx,ny,ng,nv))  
 allocate(wbuf6s(nx,ny,ng,nv))  
 allocate(wbuf1r(ng,ny,nz,nv))  
 allocate(wbuf2r(ng,ny,nz,nv))  
 allocate(wbuf3r(nx,ng,nz,nv))  
 allocate(wbuf4r(nx,ng,nz,nv))  
 allocate(wbuf5r(nx,ny,ng,nv))  
 allocate(wbuf6r(nx,ny,ng,nv))  
 allocate(divbuf1s(ng,ny,nz))  
 allocate(divbuf2s(ng,ny,nz))  
 allocate(divbuf3s(nx,ng,nz))  
 allocate(divbuf4s(nx,ng,nz))  
 allocate(divbuf5s(nx,ny,ng))  
 allocate(divbuf6s(nx,ny,ng))  
 allocate(divbuf1r(ng,ny,nz))  
 allocate(divbuf2r(ng,ny,nz))  
 allocate(divbuf3r(nx,ng,nz))  
 allocate(divbuf4r(nx,ng,nz))  
 allocate(divbuf5r(nx,ny,ng))  
 allocate(divbuf6r(nx,ny,ng))  
 allocate(ducbuf1s(ng,ny,nz))
 allocate(ducbuf2s(ng,ny,nz))
 allocate(ducbuf3s(nx,ng,nz))
 allocate(ducbuf4s(nx,ng,nz))
 allocate(ducbuf5s(nx,ny,ng))
 allocate(ducbuf6s(nx,ny,ng))
 allocate(ducbuf1r(ng,ny,nz))
 allocate(ducbuf2r(ng,ny,nz))
 allocate(ducbuf3r(nx,ng,nz))
 allocate(ducbuf4r(nx,ng,nz))
 allocate(ducbuf5r(nx,ny,ng))
 allocate(ducbuf6r(nx,ny,ng))
 allocate(yg(1-ng:nymax+ng  ))
 allocate(zg(1-ng:nzmax+ng  ))
!
 allocate(wrecyc_gpu(ng,ny,nz,nv))
 allocate(wrecycav_gpu(ng,ny,nv))
!
endsubroutine allocate_vars

subroutine copy_cpu_to_gpu()
!
 use mod_streams
 implicit none
 integer :: i,j,k,iv
!
 do iv=1,5
  do k=1-ng,nz+ng
   do j=1-ng,ny+ng
    do i=1-ng,nx+ng
     w_order(i,j,k,iv) = w(iv,i,j,k)
    enddo
   enddo
  enddo
 enddo
!
 yn_gpu           = yn
!
#ifdef USE_CUDA
 w_gpu            = w_order
!fl_gpu           = fl
!fln_gpu          = fln
!temperature_gpu  = temperature
 ducros_gpu       = ducros
 wmean_gpu        = wmean
 dcsidx_gpu       = dcsidx      
 dcsidx2_gpu      = dcsidx2      
 dcsidxs_gpu      = dcsidxs     
 detady_gpu       = detady      
 detady2_gpu      = detady2     
 detadys_gpu      = detadys     
 dzitdz_gpu       = dzitdz      
 dzitdz2_gpu      = dzitdz2     
 dzitdzs_gpu      = dzitdzs     
 dcsidxh_gpu      = dcsidxh
 detadyh_gpu      = detadyh
 dzitdzh_gpu      = dzitdzh
 x_gpu            = x
 y_gpu            = y
 z_gpu            = z
 xh_gpu           = xh
 yh_gpu           = yh
 zh_gpu           = zh
 xg_gpu           = xg
 coeff_deriv1_gpu = coeff_deriv1
 coeff_deriv1s_gpu= coeff_deriv1s
 coeff_clap_gpu   = coeff_clap
 coeff_midpi_gpu  = coeff_midpi
 cx_midpi_gpu     = cx_midpi
 cy_midpi_gpu     = cy_midpi
 cz_midpi_gpu     = cz_midpi
!fhat_gpu         = fhat
 ibcnr_gpu        = ibcnr
 dcoe_gpu         = dcoe
 winf_gpu         = winf
 winf1_gpu        = winf1
!rf_gpu           = rf
!rfy_gpu          = rfy ! not needed
 vf_df_gpu        = vf_df
 by_df_gpu        = by_df
 bz_df_gpu        = bz_df
 amat_df_gpu      = amat_df
 yplus_inflow_gpu = yplus_inflow
 yplus_recyc_gpu  = yplus_recyc
 eta_inflow_gpu   = eta_inflow
 eta_recyc_gpu    = eta_recyc
 map_j_inn_gpu    = map_j_inn
 map_j_out_gpu    = map_j_out
 weta_inflow_gpu  = weta_inflow
#else
 call move_alloc( w_order      , w_gpu           )
 call move_alloc( fl           , fl_gpu          )
 call move_alloc( fln          , fln_gpu         )
 call move_alloc( temperature  , temperature_gpu )
 call move_alloc( ducros       , ducros_gpu      )
 call move_alloc( wmean        , wmean_gpu       )
 call move_alloc( dcsidx       , dcsidx_gpu      )
 call move_alloc( dcsidx2      , dcsidx2_gpu     )
 call move_alloc( dcsidxs      , dcsidxs_gpu     )
 call move_alloc( detady       , detady_gpu      )
 call move_alloc( detady2      , detady2_gpu     )
 call move_alloc( detadys      , detadys_gpu     )
 call move_alloc( dzitdz       , dzitdz_gpu      )
 call move_alloc( dzitdz2      , dzitdz2_gpu     )
 call move_alloc( dzitdzs      , dzitdzs_gpu     )
 call move_alloc( dcsidxh      , dcsidxh_gpu     )
 call move_alloc( detadyh      , detadyh_gpu     )
 call move_alloc( dzitdzh      , dzitdzh_gpu     )
 call move_alloc( x            , x_gpu           )
 call move_alloc( y            , y_gpu           )
 call move_alloc( z            , z_gpu           )
 call move_alloc( xh           , xh_gpu          )
 call move_alloc( yh           , yh_gpu          )
 call move_alloc( zh           , zh_gpu          )
 call move_alloc( xg           , xg_gpu          )
 call move_alloc( coeff_deriv1 , coeff_deriv1_gpu)
 call move_alloc( coeff_deriv1s, coeff_deriv1s_gpu)
 call move_alloc( coeff_clap   , coeff_clap_gpu  )
 call move_alloc( coeff_midpi  , coeff_midpi_gpu )
 call move_alloc( cx_midpi     , cx_midpi_gpu    )
 call move_alloc( cy_midpi     , cy_midpi_gpu    )
 call move_alloc( cz_midpi     , cz_midpi_gpu    )
 call move_alloc( fhat         , fhat_gpu        )
 call move_alloc( ibcnr        , ibcnr_gpu       )
 call move_alloc( dcoe         , dcoe_gpu        )
 call move_alloc( winf         , winf_gpu        )
 call move_alloc( winf1        , winf1_gpu       )
 call move_alloc( rf           , rf_gpu          )
 call move_alloc( rfy          , rfy_gpu         )
 call move_alloc( vf_df        , vf_df_gpu       )
 call move_alloc( by_df        , by_df_gpu       )
 call move_alloc( bz_df        , bz_df_gpu       )
 call move_alloc( amat_df      , amat_df_gpu     )
 call move_alloc( yplus_inflow , yplus_inflow_gpu)
 call move_alloc( yplus_recyc  , yplus_recyc_gpu )
 call move_alloc( eta_inflow   , eta_inflow_gpu  )
 call move_alloc( eta_recyc    , eta_recyc_gpu   )
 call move_alloc( map_j_inn    , map_j_inn_gpu   )
 call move_alloc( map_j_out    , map_j_out_gpu   )
 call move_alloc( weta_inflow  , weta_inflow_gpu )
#endif
end subroutine copy_cpu_to_gpu

subroutine copy_gpu_to_cpu
!
 use mod_streams
 implicit none
 integer :: i,j,k,iv
!
#ifdef USE_CUDA
 w_order = w_gpu
 temperature = temperature_gpu
 vf_df = vf_df_gpu
#else
 call move_alloc(w_gpu, w_order)
 call move_alloc(temperature_gpu , temperature)
 call move_alloc(vf_df_gpu , vf_df)
 call move_alloc(x_gpu , x)
 call move_alloc(xg_gpu,xg)
 call move_alloc(y_gpu , y)
 call move_alloc(z_gpu , z)
 call move_alloc(xh_gpu,xh)
 call move_alloc(yh_gpu,yh)
 call move_alloc(zh_gpu,zh)
 call move_alloc(fl_gpu, fl)
#endif
 do iv=1,5
  do k=1-ng,nz+ng
   do j=1-ng,ny+ng
    do i=1-ng,nx+ng
     w(iv,i,j,k) = w_order(i,j,k,iv)
    enddo
   enddo
  enddo
 enddo
end subroutine copy_gpu_to_cpu

subroutine reset_cpu_gpu
!
 use mod_streams
 implicit none
!
#ifdef USE_CUDA
#else
 call move_alloc(w_order, w_gpu)
 call move_alloc(temperature, temperature_gpu)
 call move_alloc(vf_df, vf_df_gpu)
 call move_alloc(x , x_gpu)
 call move_alloc(xg,xg_gpu)
 call move_alloc(y , y_gpu)
 call move_alloc(z , z_gpu)
 call move_alloc(xh,xh_gpu)
 call move_alloc(yh,yh_gpu)
 call move_alloc(zh,zh_gpu)
 call move_alloc(fl, fl_gpu)
#endif
end subroutine reset_cpu_gpu
