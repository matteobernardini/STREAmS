module mod_streams
 use mpi
#ifdef USE_CUDA
  use cudafor
#endif
 use, intrinsic :: iso_fortran_env, only : error_unit
 implicit none
 save
!
 integer, parameter :: singtype = selected_real_kind(6,37)    ! single precision
 integer, parameter :: doubtype = selected_real_kind(15,307)  ! double precision
!integer, parameter :: quadtype = selected_real_kind(33,4931) ! extended precision
#ifdef SINGLE_PRECISION
 integer, parameter :: mykind    = singtype
 integer, parameter :: mpi_prec = mpi_real4
 real(mykind), parameter :: tol_iter = 0.00001_mykind
#else
 integer, parameter :: mykind    = doubtype
 integer, parameter :: mpi_prec = mpi_real8
 real(mykind), parameter :: tol_iter = 0.000000001_mykind
#endif
!
 integer, parameter :: nv      =  5   ! Physical variables
 integer, parameter :: nsmoo   = 10   ! Number of smoothing iterations for wall-normal mesh
 integer, parameter :: nvmean  = 20
 integer, parameter :: nsolmax = 999999
 integer, parameter :: itmax   = 100000
 real(mykind), parameter :: pi = 4._mykind*atan(1._mykind)
!
! MPI related parameters 
!
 integer, parameter :: ndims = 3
 integer, dimension(:), allocatable  :: nblocks
 logical, dimension(:), allocatable  :: pbc
 integer, dimension(mpi_status_size) :: istatus
!
 integer, dimension(:), allocatable :: ncoords
 integer :: mp_cart,mp_cartx,mp_carty,mp_cartz
 integer :: nrank,nproc,nrank_x, nrank_y, nrank_z
 integer :: ileftx,irightx,ilefty,irighty,ileftz,irightz
 integer :: iermpi, iercuda
!
 integer :: nxmax
 integer :: nymax,nymaxwr
 integer :: nzmax
 integer :: nx
 integer :: ny
 integer :: nz
 integer :: ng   ! Number of ghost nodes
 integer :: ngdf ! Number of ghost nodes for digital filtering
 integer :: io_type
!
 integer :: enable_plot3d, enable_vtk
!
! Useful code variables
!
 real(mykind), parameter :: gamma  = 1.4_mykind
 real(mykind), parameter :: pr     = 0.72_mykind
!real(mykind), parameter :: gamma  = 1.1312_mykind
!real(mykind), parameter :: pr     = 0.455_mykind
 real(mykind), parameter :: gm1    = gamma-1._mykind
 real(mykind), parameter :: gm     = 1._mykind/gm1
 real(mykind), parameter :: ggmopr = gamma*gm/pr
 real(mykind), parameter :: vtexp  = 3._mykind/4._mykind
 real(mykind), parameter :: rfac   = 0.89_mykind !pr**(1._mykind/3._mykind)
 real(mykind) :: rm,re,taw,sqgmr,s2tinf,retauinflow,trat,twall,tref_dimensional
 real(mykind) :: rtrms
 real(mykind), dimension(0:nsolmax) :: tsol, tsol_restart
 real(mykind) :: dtsave, dtsave_restart
 integer :: iflow
 integer :: idiski, ndim
 integer :: istore, istore_restart 
 integer :: iorder,iweno
 integer :: visc_type
 real(mykind) :: tresduc
 real(mykind), dimension(:,:), allocatable :: dcoe
 real(mykind), dimension(:,:), allocatable :: dcoe_gpu
 real(mykind), dimension(:), allocatable  :: winf,winf1
 real(mykind), dimension(:), allocatable  :: winf_gpu,winf1_gpu
 real(mykind) :: rho0,t0,p0,u0,v0,w0,s0
 real(mykind) :: ximp,thetas,deflec,tanhsfac,xsh
 real(mykind) :: pgradf
 real(mykind) :: heatflux,aeroheat,bulkcooling
 integer :: ivis
 logical :: masterproc
 logical :: dfupdated
!
 character(3) :: chx,chy,chz
 character(6) :: stat_io
!
! Vector of conservative variables and fluxes
 real(mykind), dimension(:,:,:,:), allocatable :: w,fl,fln
 real(mykind), dimension(:,:,:,:), allocatable :: w_order
 real(mykind), dimension(:,:,:,:), allocatable :: w_gpu,fl_gpu,fln_gpu

 real(mykind), dimension(:,:,:,:), allocatable :: fhat_trans_gpu, fl_trans_gpu
 real(mykind), dimension(:,:,:,:), allocatable :: wv_gpu, wv_trans_gpu
 real(mykind), dimension(:,:,:), allocatable :: temperature_trans_gpu

 real(mykind), dimension(:,:,:), allocatable :: temperature
 real(mykind), dimension(:,:,:), allocatable :: temperature_gpu
 logical, dimension(:,:,:), allocatable :: ducros,ducros_gpu
! 
! RK data
 real(mykind), dimension(3) :: gamvec,rhovec
 real(mykind) :: dtglobal,cfl,dtmin,alpdt,telaps,telaps0,alpdtold
 integer :: icyc,ncyc,ncyc0,nstep,nprint
!
! Coordinates and metric related quantities 
 integer :: jbgrid ! Parameter for grid stretching
 real(mykind) :: rlx,rly,rlz,rlywr,dyp_target
 real(mykind), dimension(:), allocatable :: x,y,z,yn
 real(mykind), dimension(:), allocatable :: x_gpu,y_gpu,z_gpu,yn_gpu
 real(mykind), dimension(:), allocatable :: xg,yg,zg
 real(mykind), dimension(:), allocatable :: xg_gpu
 real(mykind), dimension(:), allocatable :: dcsidx,dcsidx2,dcsidxs
 real(mykind), dimension(:), allocatable :: detady,detady2,detadys
 real(mykind), dimension(:), allocatable :: dzitdz,dzitdz2,dzitdzs
 real(mykind), dimension(:), allocatable :: dcsidx_gpu,dcsidx2_gpu,dcsidxs_gpu
 real(mykind), dimension(:), allocatable :: detady_gpu,detady2_gpu,detadys_gpu
 real(mykind), dimension(:), allocatable :: dzitdz_gpu,dzitdz2_gpu,dzitdzs_gpu
 real(mykind), dimension(:), allocatable :: dxg,dyg,dzg
!
 real(mykind), dimension(:), allocatable :: xgh,ygh,zgh
 real(mykind), dimension(:), allocatable :: xh,yh,zh
 real(mykind), dimension(:), allocatable :: xh_gpu,yh_gpu,zh_gpu
 real(mykind), dimension(:), allocatable :: dcsidxh,detadyh,dzitdzh
 real(mykind), dimension(:), allocatable :: dcsidxh_gpu,detadyh_gpu,dzitdzh_gpu
!
 real(mykind) :: dpdx,rhobulk,ubulk,tbulk
 real(mykind) :: target_tbulk
!
! BC data
 integer, dimension(:), allocatable :: ibc,ibcnr
 integer, dimension(:), allocatable :: ibcnr_gpu
! Digital filtering
 integer, parameter :: nfmax = 64
 integer :: rand_start 
 real(mykind) :: dftscaling
 real(mykind), dimension(3) :: xlen_df
 real(mykind), dimension(:,:,:), allocatable :: vf_df,vf_df_gpu
 real(mykind), dimension(:,:,:), allocatable :: rf,rf_gpu,rfy,rfy_gpu
 real(mykind), dimension(:,:,:), allocatable :: bx_df,by_df,bz_df
 real(mykind), dimension(:,:,:), allocatable :: by_df_gpu,bz_df_gpu
 real(mykind), dimension(:,:,:), allocatable :: amat_df,amat_df_gpu
! Recycling-rescaling 
 real(mykind) :: xrecyc,betarecyc
 integer :: ibrecyc ! x-coord of block with recycling station
 integer :: irecyc  ! i stat of recycling station in block ibrecyc
 real(mykind), dimension(:,:,:,:), allocatable :: wrecyc_gpu
 real(mykind), dimension(:,:,:), allocatable :: wrecycav_gpu
 real(mykind), dimension(:), allocatable :: yplus_inflow,yplus_recyc
 real(mykind), dimension(:), allocatable :: yplus_inflow_gpu,yplus_recyc_gpu
 real(mykind), dimension(:), allocatable :: eta_inflow,eta_recyc
 real(mykind), dimension(:), allocatable :: eta_inflow_gpu,eta_recyc_gpu
 real(mykind), dimension(:), allocatable :: weta_inflow,weta_inflow_gpu
 integer, dimension(:), allocatable :: map_j_inn,map_j_out
 integer, dimension(:), allocatable :: map_j_inn_gpu,map_j_out_gpu
!
! Statistical quantities
 integer :: istat,itav,nstat,nstatloc
 real(mykind), dimension(:,:,:), allocatable :: w_av,w_avzg
 real(mykind), dimension(:,:  ), allocatable :: w_av_1d,w_avxzg
 real(mykind), dimension(:), allocatable     ::  xstat
 integer     , dimension(:), allocatable     :: ixstat
 integer     , dimension(:), allocatable     :: igxstat
!
! Mean field for initialization
 real(mykind), dimension(:,:,:), allocatable :: wmean,wmean_gpu
!
 real(mykind), dimension(:), allocatable :: coeff_deriv1
 real(mykind), dimension(:), allocatable :: coeff_deriv1_gpu
 real(mykind), dimension(:), allocatable :: coeff_clap
 real(mykind), dimension(:), allocatable :: coeff_clap_gpu
 real(mykind), dimension(:), allocatable :: coeff_deriv1s
 real(mykind), dimension(:), allocatable :: coeff_deriv1s_gpu
 real(mykind), dimension(:), allocatable :: coeff_midpi
 real(mykind), dimension(:), allocatable :: coeff_midpi_gpu
 real(mykind), dimension(:,:), allocatable :: cx_midpi,cy_midpi,cz_midpi
 real(mykind), dimension(:,:), allocatable :: cx_midpi_gpu,cy_midpi_gpu,cz_midpi_gpu
!
 real(mykind), dimension(:,:,:,:), allocatable :: wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
 real(mykind), dimension(:,:,:,:), allocatable :: wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
 real(mykind), dimension(:,:,:), allocatable :: divbuf1s_gpu, divbuf2s_gpu, divbuf3s_gpu, divbuf4s_gpu, divbuf5s_gpu, divbuf6s_gpu
 real(mykind), dimension(:,:,:), allocatable :: divbuf1r_gpu, divbuf2r_gpu, divbuf3r_gpu, divbuf4r_gpu, divbuf5r_gpu, divbuf6r_gpu
 logical, dimension(:,:,:), allocatable :: ducbuf1s_gpu, ducbuf2s_gpu, ducbuf3s_gpu, ducbuf4s_gpu, ducbuf5s_gpu, ducbuf6s_gpu
 logical, dimension(:,:,:), allocatable :: ducbuf1r_gpu, ducbuf2r_gpu, ducbuf3r_gpu, ducbuf4r_gpu, ducbuf5r_gpu, ducbuf6r_gpu
!
 real(mykind), dimension(:,:,:,:), allocatable :: wbuf1s, wbuf2s, wbuf3s, wbuf4s, wbuf5s, wbuf6s
 real(mykind), dimension(:,:,:,:), allocatable :: wbuf1r, wbuf2r, wbuf3r, wbuf4r, wbuf5r, wbuf6r
 real(mykind), dimension(:,:,:), allocatable :: divbuf1s, divbuf2s, divbuf3s, divbuf4s, divbuf5s, divbuf6s
 real(mykind), dimension(:,:,:), allocatable :: divbuf1r, divbuf2r, divbuf3r, divbuf4r, divbuf5r, divbuf6r
 logical, dimension(:,:,:), allocatable :: ducbuf1s, ducbuf2s, ducbuf3s, ducbuf4s, ducbuf5s, ducbuf6s
 logical, dimension(:,:,:), allocatable :: ducbuf1r, ducbuf2r, ducbuf3r, ducbuf4r, ducbuf5r, ducbuf6r
!
 real(mykind), dimension(:,:), allocatable :: wallpfield,wallpfield_gpu
 real(mykind), dimension(:,:,:), allocatable :: slicexy,slicexy_gpu
 real(mykind), dimension(:,:,:,:), allocatable :: fhat,fhat_gpu
!
 real(mykind), allocatable, dimension(:,:,:) :: vf_df_old
 real(mykind), allocatable, dimension(:,:,:) :: uf
 real(mykind), allocatable, dimension(:,:) :: evmax_mat_yz
 real(mykind), allocatable, dimension(:) :: evmax_mat_y
 real(mykind), allocatable, dimension(:) :: bulk5g_gpu
 real(mykind), allocatable, dimension(:,:) :: rtrms_ib_gpu
 real(mykind), allocatable, dimension(:) :: rtrms_ib_1d_gpu
!
 real(mykind), dimension(:,:,:,:), allocatable :: gplus_x,gminus_x
 real(mykind), dimension(:,:,:,:), allocatable :: gplus_y,gminus_y
 real(mykind), dimension(:,:,:,:), allocatable :: gplus_z,gminus_z
!
#ifdef USE_CUDA
 attributes(device) :: fhat_trans_gpu, fl_trans_gpu
 attributes(device) :: temperature_trans_gpu
 attributes(device) :: wv_gpu, wv_trans_gpu
!
 integer :: local_comm, mydev
 attributes(device) :: w_gpu,fl_gpu,fln_gpu
 attributes(device) :: temperature_gpu,ducros_gpu
 attributes(device) :: dcsidx_gpu,dcsidx2_gpu,dcsidxs_gpu
 attributes(device) :: detady_gpu,detady2_gpu,detadys_gpu
 attributes(device) :: dzitdz_gpu,dzitdz2_gpu,dzitdzs_gpu
 attributes(device) :: dcsidxh_gpu,detadyh_gpu,dzitdzh_gpu
 attributes(device) :: x_gpu,y_gpu,yn_gpu,z_gpu
 attributes(device) :: xh_gpu,yh_gpu,zh_gpu
 attributes(device) :: xg_gpu
 attributes(device) :: coeff_deriv1_gpu
 attributes(device) :: coeff_clap_gpu
 attributes(device) :: coeff_deriv1s_gpu
 attributes(device) :: coeff_midpi_gpu
 attributes(device) :: cx_midpi_gpu,cy_midpi_gpu,cz_midpi_gpu
 attributes(device) :: ibcnr_gpu
 attributes(device) :: dcoe_gpu
 attributes(device) :: wmean_gpu
 attributes(device) :: winf_gpu,winf1_gpu
 attributes(device) :: rf_gpu,rfy_gpu
 attributes(device) :: vf_df_gpu
 attributes(device) :: by_df_gpu
 attributes(device) :: bz_df_gpu
 attributes(device) :: amat_df_gpu
 attributes(device) :: fhat_gpu
 attributes(device) :: wbuf1s_gpu, wbuf2s_gpu, wbuf3s_gpu, wbuf4s_gpu, wbuf5s_gpu, wbuf6s_gpu
 attributes(device) :: wbuf1r_gpu, wbuf2r_gpu, wbuf3r_gpu, wbuf4r_gpu, wbuf5r_gpu, wbuf6r_gpu
 attributes(device) :: divbuf1s_gpu, divbuf2s_gpu, divbuf3s_gpu, divbuf4s_gpu, divbuf5s_gpu, divbuf6s_gpu
 attributes(device) :: divbuf1r_gpu, divbuf2r_gpu, divbuf3r_gpu, divbuf4r_gpu, divbuf5r_gpu, divbuf6r_gpu
 attributes(device) :: ducbuf1s_gpu, ducbuf2s_gpu, ducbuf3s_gpu, ducbuf4s_gpu, ducbuf5s_gpu, ducbuf6s_gpu
 attributes(device) :: ducbuf1r_gpu, ducbuf2r_gpu, ducbuf3r_gpu, ducbuf4r_gpu, ducbuf5r_gpu, ducbuf6r_gpu
 attributes(device) :: wrecyc_gpu,wrecycav_gpu
 attributes(device) :: yplus_inflow_gpu,yplus_recyc_gpu
 attributes(device) :: eta_inflow_gpu,eta_recyc_gpu
 attributes(device) :: map_j_inn_gpu,map_j_out_gpu
 attributes(device) :: weta_inflow_gpu

 attributes(pinned) :: wbuf1s, wbuf2s, wbuf3s, wbuf4s, wbuf5s, wbuf6s
 attributes(pinned) :: wbuf1r, wbuf2r, wbuf3r, wbuf4r, wbuf5r, wbuf6r
 attributes(pinned) :: divbuf1s, divbuf2s, divbuf3s, divbuf4s, divbuf5s, divbuf6s
 attributes(pinned) :: divbuf1r, divbuf2r, divbuf3r, divbuf4r, divbuf5r, divbuf6r
 attributes(pinned) :: ducbuf1s, ducbuf2s, ducbuf3s, ducbuf4s, ducbuf5s, ducbuf6s
 attributes(pinned) :: ducbuf1r, ducbuf2r, ducbuf3r, ducbuf4r, ducbuf5r, ducbuf6r

 integer(kind=cuda_stream_kind) :: stream1, stream2
 attributes(device) :: vf_df_old,uf
 attributes(device) :: evmax_mat_yz,evmax_mat_y
 attributes(device) :: bulk5g_gpu
 attributes(device) :: rtrms_ib_gpu,rtrms_ib_1d_gpu
 attributes(device) :: wallpfield_gpu
 attributes(device) :: slicexy_gpu

 attributes(device) :: gplus_x, gminus_x
 attributes(device) :: gplus_y, gminus_y
 attributes(device) :: gplus_z, gminus_z
#endif

end module mod_streams
