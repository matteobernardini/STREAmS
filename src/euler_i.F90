subroutine euler_i(istart, iend)
!
! Evaluation of convective terms in the x direction
!
 use mod_streams
 implicit none
!
 integer :: istart, iend
!
 real(mykind) :: evmax1,evmax2,evmax3,evmax4,evmax5
 real(mykind) :: gplus1_1,gplus1_2,gplus1_3,gplus1_4
 real(mykind) :: gplus2_1,gplus2_2,gplus2_3,gplus2_4
 real(mykind) :: gplus3_1,gplus3_2,gplus3_3,gplus3_4
 real(mykind) :: gplus4_1,gplus4_2,gplus4_3,gplus4_4
 real(mykind) :: gplus5_1,gplus5_2,gplus5_3,gplus5_4
 real(mykind) :: gminus1_1,gminus1_2,gminus1_3,gminus1_4
 real(mykind) :: gminus2_1,gminus2_2,gminus2_3,gminus2_4
 real(mykind) :: gminus3_1,gminus3_2,gminus3_3,gminus3_4
 real(mykind) :: gminus4_1,gminus4_2,gminus4_3,gminus4_4
 real(mykind) :: gminus5_1,gminus5_2,gminus5_3,gminus5_4
!
 real(mykind) :: el_1_1,el_1_2,el_1_3,el_1_4,el_1_5
 real(mykind) :: el_2_1,el_2_2,el_2_3,el_2_4,el_2_5
 real(mykind) :: el_3_1,el_3_2,el_3_3,el_3_4,el_3_5
 real(mykind) :: el_4_1,el_4_2,el_4_3,el_4_4,el_4_5
 real(mykind) :: el_5_1,el_5_2,el_5_3,el_5_4,el_5_5
 real(mykind) :: er_1_1,er_1_2,er_1_3,er_1_4,er_1_5
 real(mykind) :: er_2_1,er_2_2,er_2_3,er_2_4,er_2_5
 real(mykind) :: er_3_1,er_3_2,er_3_3,er_3_4,er_3_5
 real(mykind) :: er_4_1,er_4_2,er_4_3,er_4_4,er_4_5
 real(mykind) :: er_5_1,er_5_2,er_5_3,er_5_4,er_5_5
 real(mykind) :: gr1,gr2,gr3,gr4,gr5
 real(mykind) :: gl1,gl2,gl3,gl4,gl5
 real(mykind) :: ghat1,ghat2,ghat3,ghat4,ghat5
 real(mykind) :: fi1,fi2,fi3,fi4,fi5
 real(mykind) :: ev1,ev2,ev3,ev4,ev5
!
!real(mykind),dimension(6,1-ng:nx+ng) :: u,v
!real(mykind),dimension(6,ng,1-ng:nx+ng) :: uv
 real(mykind) :: uui,uuip,uuip2
 real(mykind) :: vvi,vvip,vvip2
 real(mykind) :: wwi,wwip,wwip2
 real(mykind) :: ppi,ppip,ppip2
 real(mykind) :: enti,entip,entip2
 real(mykind) :: uv_1_1_i,uv_1_2_i
 real(mykind) :: uv_2_1_i,uv_2_2_i
 real(mykind) :: uv_3_1_i,uv_3_2_i
 real(mykind) :: uv_4_1_i,uv_4_2_i
 real(mykind) :: uv_5_1_i,uv_5_2_i
 real(mykind) :: uv_6_1_i,uv_6_2_i
 real(mykind) :: uv_1_1_im,uv_1_2_im
 real(mykind) :: uv_2_1_im,uv_2_2_im
 real(mykind) :: uv_3_1_im,uv_3_2_im
 real(mykind) :: uv_4_1_im,uv_4_2_im
 real(mykind) :: uv_5_1_im,uv_5_2_im
 real(mykind) :: uv_6_1_im,uv_6_2_im
!real(mykind),dimension(5) :: fh
!real(mykind),dimension(6) :: ft,uvs
 real(mykind) :: fh1,fh2,fh3,fh4,fh5
 real(mykind) :: ft1,ft2,ft3,ft4,ft5,ft6
 real(mykind) :: uvs1,uvs2,uvs3,uvs4,uvs5,uvs6
 integer :: endi,endj,endk
!
 integer :: i,j,k,l,m,n
 integer :: im,ip,ll,lmax,mm,ng2
 real(mykind) :: b1,b2,c,cc,ci,ducm,ent,gc,gpr,h,hp,pp,qq,qqp
 real(mykind) :: r,rho,rhoe,rhom,rhou,rhov,rhow,ri,rip,rp1,up
 real(mykind) :: uu,vp,vv,wc,wp,ww
 real(mykind) :: df
!----------------------------------------------------------------------------------------
! from wenorec
!----------------------------------------------------------------------------------------
 integer :: ii
 integer, parameter :: nvar = 5
!real(mykind),dimension(nvar,2*ng) :: vm_weno,vp_weno
!real(mykind),dimension(nvar) :: vminus,vplus
 real(mykind) :: vm_weno1_1,vm_weno1_2,vm_weno1_3,vm_weno1_4
 real(mykind) :: vm_weno2_1,vm_weno2_2,vm_weno2_3,vm_weno2_4
 real(mykind) :: vm_weno3_1,vm_weno3_2,vm_weno3_3,vm_weno3_4
 real(mykind) :: vm_weno4_1,vm_weno4_2,vm_weno4_3,vm_weno4_4
 real(mykind) :: vm_weno5_1,vm_weno5_2,vm_weno5_3,vm_weno5_4
 real(mykind) :: vp_weno1_1,vp_weno1_2,vp_weno1_3,vp_weno1_4
 real(mykind) :: vp_weno2_1,vp_weno2_2,vp_weno2_3,vp_weno2_4
 real(mykind) :: vp_weno3_1,vp_weno3_2,vp_weno3_3,vp_weno3_4
 real(mykind) :: vp_weno4_1,vp_weno4_2,vp_weno4_3,vp_weno4_4
 real(mykind) :: vp_weno5_1,vp_weno5_2,vp_weno5_3,vp_weno5_4
 real(mykind) :: vminus1,vplus1,vminus2,vplus2,vminus3,vplus3,vminus4,vplus4,vminus5,vplus5
!real(mykind),dimension(-1:4) :: dwe           ! linear weights
!real(mykind),dimension(-1:4) :: alfp,alfm     ! alpha_l
!real(mykind),dimension(-1:4) :: betap,betam   ! beta_l
!real(mykind),dimension(-1:4) :: omp,omm       ! WENO weights
 real(mykind) :: dwe0           ! linear weights
 real(mykind) :: alfp0,alfm0     ! alpha_l
 real(mykind) :: betap0,betam0   ! beta_l
 real(mykind) :: omp0,omm0       ! WENO weights
 real(mykind) :: dwe1           ! linear weights
 real(mykind) :: alfp1,alfm1     ! alpha_l
 real(mykind) :: betap1,betam1   ! beta_l
 real(mykind) :: omp1,omm1       ! WENO weights
 real(mykind) :: sump, summ, temp
 real(mykind) :: eps
 logical :: ishk
!----------------------------------------------------------------------------------------
!
 endi = nx
 endj = ny
 endk = nz
 ng2  = iorder
 lmax = iorder/2 ! max stencil width
 if (ibc(2)==4.or.ibc(2)==8) then
  endi = nx-1
 endif
!
 eps = 0.000001_mykind
!
 !$cuf kernel do(2) <<<*,*,stream=stream1>>> 
 do k=1,endk ! loop in the z-direction
  do j=1,endj ! loop in the y-direction
!
!  l=1
   i = istart-1
   rhom     = w_gpu(1,i,j,k)  +w_gpu(1,i+1,j,k)
   uui      = w_gpu(2,i,j,k)  /w_gpu(1,i,j,k)
   uuip     = w_gpu(2,i+1,j,k)/w_gpu(1,i+1,j,k)
   vvi      = w_gpu(3,i,j,k)  /w_gpu(1,i,j,k)
   vvip     = w_gpu(3,i+1,j,k)/w_gpu(1,i+1,j,k)
   wwi      = w_gpu(4,i,j,k)  /w_gpu(1,i,j,k)
   wwip     = w_gpu(4,i+1,j,k)/w_gpu(1,i+1,j,k)
   ppi      = w_gpu(1,i,j,k)  *temperature_gpu(i,j,k)
   ppip     = w_gpu(1,i+1,j,k)*temperature_gpu(i+1,j,k)
   enti     = (w_gpu(5,i,j,k)+ppi)/w_gpu(1,i,j,k)
   entip    = (w_gpu(5,i+1,j,k)+ppip)/w_gpu(1,i+1,j,k)
   uv_1_1_i = (uui+uuip)*(             2._mykind)*rhom
   uv_2_1_i = (uui+uuip)*(uui+uuip)*rhom
   uv_3_1_i = (uui+uuip)*(vvi+vvip)*rhom
   uv_4_1_i = (uui+uuip)*(wwi+wwip)*rhom
   uv_5_1_i = (uui+uuip)*(enti+entip)*rhom
   uv_6_1_i = (   +2._mykind  )*(ppi+ppip)
!  l=2
   rhom     = w_gpu(1,i,j,k)+w_gpu(1,i+2,j,k)
   uuip2    = w_gpu(2,i+2,j,k)/w_gpu(1,i+2,j,k)
   vvip2    = w_gpu(3,i+2,j,k)/w_gpu(1,i+2,j,k)
   wwip2    = w_gpu(4,i+2,j,k)/w_gpu(1,i+2,j,k)
   ppip2    = w_gpu(1,i+2,j,k)*temperature_gpu(i+2,j,k)
   entip2   = (w_gpu(5,i+2,j,k)+ppip2)/w_gpu(1,i+2,j,k)
   uv_1_2_i = (uui+uuip2)*(             2._mykind)*rhom
   uv_2_2_i = (uui+uuip2)*(uui+uuip2)*rhom
   uv_3_2_i = (uui+uuip2)*(vvi+vvip2)*rhom
   uv_4_2_i = (uui+uuip2)*(wwi+wwip2)*rhom
   uv_5_2_i = (uui+uuip2)*(enti+entip2)*rhom
   uv_6_2_i = (       2._mykind)*(ppi+ppip2)
!
!  Evaluation of numerical fluxes
!
   do i=istart,iend ! i is the index of intermediate nodes 
!
    uv_1_1_im = uv_1_1_i
    uv_2_1_im = uv_2_1_i
    uv_3_1_im = uv_3_1_i
    uv_4_1_im = uv_4_1_i
    uv_5_1_im = uv_5_1_i
    uv_6_1_im = uv_6_1_i
    uv_1_2_im = uv_1_2_i
    uv_2_2_im = uv_2_2_i
    uv_3_2_im = uv_3_2_i
    uv_4_2_im = uv_4_2_i
    uv_5_2_im = uv_5_2_i
    uv_6_2_im = uv_6_2_i
!   l=1
    rhom     = w_gpu(1,i,j,k)  +w_gpu(1,i+1,j,k)
    uui      = w_gpu(2,i,j,k)  /w_gpu(1,i,j,k)
    uuip     = w_gpu(2,i+1,j,k)/w_gpu(1,i+1,j,k)
    vvi      = w_gpu(3,i,j,k)  /w_gpu(1,i,j,k)
    vvip     = w_gpu(3,i+1,j,k)/w_gpu(1,i+1,j,k)
    wwi      = w_gpu(4,i,j,k)  /w_gpu(1,i,j,k)
    wwip     = w_gpu(4,i+1,j,k)/w_gpu(1,i+1,j,k)
    ppi      = w_gpu(1,i,j,k)  *temperature_gpu(i,j,k)
    ppip     = w_gpu(1,i+1,j,k)*temperature_gpu(i+1,j,k)
    enti     = (w_gpu(5,i,j,k)+ppi)/w_gpu(1,i,j,k)
    entip    = (w_gpu(5,i+1,j,k)+ppip)/w_gpu(1,i+1,j,k)
    uv_1_1_i = (uui+uuip)*(             2._mykind)*rhom
    uv_2_1_i = (uui+uuip)*(uui+uuip)*rhom
    uv_3_1_i = (uui+uuip)*(vvi+vvip)*rhom
    uv_4_1_i = (uui+uuip)*(wwi+wwip)*rhom
    uv_5_1_i = (uui+uuip)*(enti+entip)*rhom
    uv_6_1_i = (   +2._mykind  )*(ppi+ppip)
!   l=2
    rhom     = w_gpu(1,i,j,k)+w_gpu(1,i+2,j,k)
    uuip2    = w_gpu(2,i+2,j,k)/w_gpu(1,i+2,j,k)
    vvip2    = w_gpu(3,i+2,j,k)/w_gpu(1,i+2,j,k)
    wwip2    = w_gpu(4,i+2,j,k)/w_gpu(1,i+2,j,k)
    ppip2    = w_gpu(1,i+2,j,k)*temperature_gpu(i+2,j,k)
    entip2   = (w_gpu(5,i+2,j,k)+ppip2)/w_gpu(1,i+2,j,k)
    uv_1_2_i = (uui+uuip2)*(             2._mykind)*rhom
    uv_2_2_i = (uui+uuip2)*(uui+uuip2)*rhom
    uv_3_2_i = (uui+uuip2)*(vvi+vvip2)*rhom
    uv_4_2_i = (uui+uuip2)*(wwi+wwip2)*rhom
    uv_5_2_i = (uui+uuip2)*(enti+entip2)*rhom
    uv_6_2_i = (       2._mykind)*(ppi+ppip2)
!
    ishk = .false.
    do ii=i-iweno+1,i+iweno
     if (ducros_gpu(ii,j,k)) ishk = .true.
    enddo
!
    if (.not.ishk) then
!
     ft1 = 0._mykind
     ft2 = 0._mykind
     ft3 = 0._mykind
     ft4 = 0._mykind
     ft5 = 0._mykind
     ft6 = 0._mykind
! l=1
     uvs1 = 0._mykind
     uvs2 = 0._mykind
     uvs3 = 0._mykind
     uvs4 = 0._mykind
     uvs5 = 0._mykind
     uvs6 = 0._mykind
     uvs1 = uvs1 + uv_1_1_i
     uvs2 = uvs2 + uv_2_1_i
     uvs3 = uvs3 + uv_3_1_i
     uvs4 = uvs4 + uv_4_1_i
     uvs5 = uvs5 + uv_5_1_i
     uvs6 = uvs6 + uv_6_1_i
     ft1 = ft1 + dcoe_gpu(1,lmax)*uvs1
     ft2 = ft2 + dcoe_gpu(1,lmax)*uvs2
     ft3 = ft3 + dcoe_gpu(1,lmax)*uvs3
     ft4 = ft4 + dcoe_gpu(1,lmax)*uvs4
     ft5 = ft5 + dcoe_gpu(1,lmax)*uvs5
     ft6 = ft6 + dcoe_gpu(1,lmax)*uvs6
! l=2
     uvs1 = 0._mykind
     uvs2 = 0._mykind
     uvs3 = 0._mykind
     uvs4 = 0._mykind
     uvs5 = 0._mykind
     uvs6 = 0._mykind
     uvs1 = uvs1 + uv_1_2_i
     uvs2 = uvs2 + uv_2_2_i
     uvs3 = uvs3 + uv_3_2_i
     uvs4 = uvs4 + uv_4_2_i
     uvs5 = uvs5 + uv_5_2_i
     uvs6 = uvs6 + uv_6_2_i
     uvs1 = uvs1 + uv_1_2_im
     uvs2 = uvs2 + uv_2_2_im
     uvs3 = uvs3 + uv_3_2_im
     uvs4 = uvs4 + uv_4_2_im
     uvs5 = uvs5 + uv_5_2_im
     uvs6 = uvs6 + uv_6_2_im
     ft1 = ft1 + dcoe_gpu(2,lmax)*uvs1
     ft2 = ft2 + dcoe_gpu(2,lmax)*uvs2
     ft3 = ft3 + dcoe_gpu(2,lmax)*uvs3
     ft4 = ft4 + dcoe_gpu(2,lmax)*uvs4
     ft5 = ft5 + dcoe_gpu(2,lmax)*uvs5
     ft6 = ft6 + dcoe_gpu(2,lmax)*uvs6
!
     fh1 = 0.25_mykind*ft1
     fh2 = 0.25_mykind*ft2
     fh3 = 0.25_mykind*ft3
     fh4 = 0.25_mykind*ft4
     fh5 = 0.25_mykind*ft5
     fh2 = fh2 + 0.5_mykind*ft6
!
     fhat_gpu(1,i,j,k) = fh1
     fhat_gpu(2,i,j,k) = fh2
     fhat_gpu(3,i,j,k) = fh3
     fhat_gpu(4,i,j,k) = fh4
     fhat_gpu(5,i,j,k) = fh5

    else
!
     ip = i + 1
!
!    Compute Roe average
!
!    Left state (node i)
     ri        =  1._mykind/w_gpu(1,i,j,k)
     uu        =  w_gpu(2,i,j,k) * ri
     vv        =  w_gpu(3,i,j,k) * ri
     ww        =  w_gpu(4,i,j,k) * ri
     qq        =  0.5_mykind * (uu*uu  +vv*vv + ww*ww)
!    pp        =  gm1 * (w_gpu(5,i,j,k) - w_gpu(1,i,j,k) * qq)
!    h         =  (w_gpu(5,i,j,k)  + pp) * ri
     h         =  gamma*w_gpu(5,i,j,k)*ri-gm1*qq
!    Right state (node i+1)
     rip       =  1._mykind/w_gpu(1,ip,j,k)
     up        =  w_gpu(2,ip,j,k) * rip
     vp        =  w_gpu(3,ip,j,k) * rip
     wp        =  w_gpu(4,ip,j,k) * rip
     qqp       =  0.5_mykind * (up*up  +vp*vp +wp*wp)
!    ppp       =  gm1 * (w_gpu(5,ip,j,k) - w_gpu(1,ip,j,k) * qqp)
!    hp        =  (w_gpu(5,ip,j,k)  + ppp) * rip
     hp        =  gamma*w_gpu(5,ip,j,k)*rip-gm1*qqp
!    average state
     r         =  w_gpu(1,ip,j,k)*ri
     r         =  sqrt(r)
     rp1       =  1._mykind/(r  +1._mykind)
     uu        =  (r*up  +uu)*rp1
     vv        =  (r*vp  +vv)*rp1
     ww        =  (r*wp  +ww)*rp1
     h         =  (r*hp  +h)*rp1
     qq        =  0.5_mykind * (uu*uu  +vv*vv + ww*ww)
     cc        =  gm1 * (h - qq)
     c         =  sqrt(cc) 
     ci        =  1._mykind/c 
!
!    left eigenvectors matrix (at Roe state)
!
     b2        =   gm1/cc
     b1        =   b2 * qq
     el_1_1  =   0.5_mykind * (b1     + uu * ci)
     el_2_1  =  -0.5_mykind * (b2 * uu +     ci)
     el_3_1  =  -0.5_mykind * (b2 * vv         )
     el_4_1  =  -0.5_mykind * (b2 * ww         )
     el_5_1  =   0.5_mykind * b2
     el_1_2  =   1._mykind - b1
     el_2_2  =   b2*uu
     el_3_2  =   b2*vv
     el_4_2  =   b2*ww
     el_5_2  =  -b2
     el_1_3  =   0.5_mykind * (b1     - uu * ci)
     el_2_3  =  -0.5_mykind * (b2 * uu -     ci)
     el_3_3  =  -0.5_mykind * (b2 * vv         )
     el_4_3  =  -0.5_mykind * (b2 * ww         )
     el_5_3  =   0.5_mykind * b2
     el_1_4  =  -vv
     el_2_4  =   0._mykind
     el_3_4  =   1._mykind
     el_4_4  =   0._mykind
     el_5_4  =   0._mykind
     el_1_5  =  -ww
     el_2_5  =   0._mykind
     el_3_5  =   0._mykind
     el_4_5  =   1._mykind
     el_5_5  =   0._mykind
!
!    right eigenvectors matrix (at Roe state)
!
     er_1_1  =  1._mykind
     er_2_1  =  1._mykind
     er_3_1  =  1._mykind
     er_4_1  =  0._mykind
     er_5_1  =  0._mykind
     er_1_2  =  uu -  c
     er_2_2  =  uu
     er_3_2  =  uu +  c
     er_4_2  =  0._mykind
     er_5_2  =  0._mykind
     er_1_3  =  vv
     er_2_3  =  vv
     er_3_3  =  vv
     er_4_3  =  1._mykind
     er_5_3  =  0._mykind
     er_1_4  =  ww
     er_2_4  =  ww
     er_3_4  =  ww
     er_4_4  =  0._mykind
     er_5_4  =  1._mykind
     er_1_5  =  h  - uu * c
     er_2_5  =  qq
     er_3_5  =  h  + uu * c
     er_4_5  =  vv
     er_5_5  =  ww
!
     evmax1 = -1._mykind
     evmax2 = -1._mykind
     evmax3 = -1._mykind
     evmax4 = -1._mykind
     evmax5 = -1._mykind
     do l=1,ng2 ! LLF
      ll = i + l - iorder/2
!
      rho  = w_gpu(1,ll,j,k)
      rhou = w_gpu(2,ll,j,k)
      ri   = 1._mykind/rho
      uu   = rhou*ri
      pp   = rho*temperature_gpu(ll,j,k)
      gpr  = gamma * pp * ri
      c    = sqrt (gpr)
      ev1 = abs(uu-c)
      ev2 = abs(uu)
      ev3 = abs(uu+c)
      ev4 = ev2
      ev5 = ev2
!
      evmax1 = max(ev1,evmax1)
      evmax2 = max(ev2,evmax2)
      evmax3 = max(ev3,evmax3)
      evmax4 = max(ev4,evmax4)
      evmax5 = max(ev5,evmax5)
     enddo
!
!--------------------------
     l  = 1
     ll = i + l - iorder/2
!
     rho  = w_gpu(1,ll,j,k)
     rhou = w_gpu(2,ll,j,k)
     ri   = 1._mykind/rho
     uu   = rhou*ri
     pp   = rho*temperature_gpu(ll,j,k)
     fi1  =       w_gpu(2,ll,j,k)
     fi2  = uu *  w_gpu(2,ll,j,k)  + pp
     fi3  = uu *  w_gpu(3,ll,j,k)
     fi4  = uu *  w_gpu(4,ll,j,k)
     fi5  = uu * (w_gpu(5,ll,j,k)  + pp)
!
     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_1 * w_gpu(1,ll,j,k)
     gc = gc + el_1_1 * fi1
     wc = wc + el_2_1 * w_gpu(2,ll,j,k)
     gc = gc + el_2_1 * fi2
     wc = wc + el_3_1 * w_gpu(3,ll,j,k)
     gc = gc + el_3_1 * fi3
     wc = wc + el_4_1 * w_gpu(4,ll,j,k)
     gc = gc + el_4_1 * fi4
     wc = wc + el_5_1 * w_gpu(5,ll,j,k)
     gc = gc + el_5_1 * fi5
     gplus1_1  = 0.5_mykind * (gc + evmax1 * wc)
     gminus1_1 = gc - gplus1_1

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_2 * w_gpu(1,ll,j,k)
     gc = gc + el_1_2 * fi1
     wc = wc + el_2_2 * w_gpu(2,ll,j,k)
     gc = gc + el_2_2 * fi2
     wc = wc + el_3_2 * w_gpu(3,ll,j,k)
     gc = gc + el_3_2 * fi3
     wc = wc + el_4_2 * w_gpu(4,ll,j,k)
     gc = gc + el_4_2 * fi4
     wc = wc + el_5_2 * w_gpu(5,ll,j,k)
     gc = gc + el_5_2 * fi5
     gplus2_1  = 0.5_mykind * (gc + evmax2 * wc)
     gminus2_1 = gc - gplus2_1

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_3 * w_gpu(1,ll,j,k)
     gc = gc + el_1_3 * fi1
     wc = wc + el_2_3 * w_gpu(2,ll,j,k)
     gc = gc + el_2_3 * fi2
     wc = wc + el_3_3 * w_gpu(3,ll,j,k)
     gc = gc + el_3_3 * fi3
     wc = wc + el_4_3 * w_gpu(4,ll,j,k)
     gc = gc + el_4_3 * fi4
     wc = wc + el_5_3 * w_gpu(5,ll,j,k)
     gc = gc + el_5_3 * fi5
     gplus3_1  = 0.5_mykind * (gc + evmax3 * wc)
     gminus3_1 = gc - gplus3_1

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_4 * w_gpu(1,ll,j,k)
     gc = gc + el_1_4 * fi1
     wc = wc + el_2_4 * w_gpu(2,ll,j,k)
     gc = gc + el_2_4 * fi2
     wc = wc + el_3_4 * w_gpu(3,ll,j,k)
     gc = gc + el_3_4 * fi3
     wc = wc + el_4_4 * w_gpu(4,ll,j,k)
     gc = gc + el_4_4 * fi4
     wc = wc + el_5_4 * w_gpu(5,ll,j,k)
     gc = gc + el_5_4 * fi5
     gplus4_1  = 0.5_mykind * (gc + evmax4 * wc)
     gminus4_1 = gc - gplus4_1

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_5 * w_gpu(1,ll,j,k)
     gc = gc + el_1_5 * fi1
     wc = wc + el_2_5 * w_gpu(2,ll,j,k)
     gc = gc + el_2_5 * fi2
     wc = wc + el_3_5 * w_gpu(3,ll,j,k)
     gc = gc + el_3_5 * fi3
     wc = wc + el_4_5 * w_gpu(4,ll,j,k)
     gc = gc + el_4_5 * fi4
     wc = wc + el_5_5 * w_gpu(5,ll,j,k)
     gc = gc + el_5_5 * fi5
     gplus5_1  = 0.5_mykind * (gc + evmax5 * wc)
     gminus5_1 = gc - gplus5_1
!---------------------------------------------
     l  = 2
     ll = i + l - iorder/2
!
     rho  = w_gpu(1,ll,j,k)
     rhou = w_gpu(2,ll,j,k)
     ri   = 1._mykind/rho
     uu   = rhou*ri
     pp   = rho*temperature_gpu(ll,j,k)
     fi1  =       w_gpu(2,ll,j,k)
     fi2  = uu *  w_gpu(2,ll,j,k)  + pp
     fi3  = uu *  w_gpu(3,ll,j,k)
     fi4  = uu *  w_gpu(4,ll,j,k)
     fi5  = uu * (w_gpu(5,ll,j,k)  + pp)
!
     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_1 * w_gpu(1,ll,j,k)
     gc = gc + el_1_1 * fi1
     wc = wc + el_2_1 * w_gpu(2,ll,j,k)
     gc = gc + el_2_1 * fi2
     wc = wc + el_3_1 * w_gpu(3,ll,j,k)
     gc = gc + el_3_1 * fi3
     wc = wc + el_4_1 * w_gpu(4,ll,j,k)
     gc = gc + el_4_1 * fi4
     wc = wc + el_5_1 * w_gpu(5,ll,j,k)
     gc = gc + el_5_1 * fi5
     gplus1_2  = 0.5_mykind * (gc + evmax1 * wc)
     gminus1_2 = gc - gplus1_2

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_2 * w_gpu(1,ll,j,k)
     gc = gc + el_1_2 * fi1
     wc = wc + el_2_2 * w_gpu(2,ll,j,k)
     gc = gc + el_2_2 * fi2
     wc = wc + el_3_2 * w_gpu(3,ll,j,k)
     gc = gc + el_3_2 * fi3
     wc = wc + el_4_2 * w_gpu(4,ll,j,k)
     gc = gc + el_4_2 * fi4
     wc = wc + el_5_2 * w_gpu(5,ll,j,k)
     gc = gc + el_5_2 * fi5
     gplus2_2  = 0.5_mykind * (gc + evmax2 * wc)
     gminus2_2 = gc - gplus2_2

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_3 * w_gpu(1,ll,j,k)
     gc = gc + el_1_3 * fi1
     wc = wc + el_2_3 * w_gpu(2,ll,j,k)
     gc = gc + el_2_3 * fi2
     wc = wc + el_3_3 * w_gpu(3,ll,j,k)
     gc = gc + el_3_3 * fi3
     wc = wc + el_4_3 * w_gpu(4,ll,j,k)
     gc = gc + el_4_3 * fi4
     wc = wc + el_5_3 * w_gpu(5,ll,j,k)
     gc = gc + el_5_3 * fi5
     gplus3_2  = 0.5_mykind * (gc + evmax3 * wc)
     gminus3_2 = gc - gplus3_2

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_4 * w_gpu(1,ll,j,k)
     gc = gc + el_1_4 * fi1
     wc = wc + el_2_4 * w_gpu(2,ll,j,k)
     gc = gc + el_2_4 * fi2
     wc = wc + el_3_4 * w_gpu(3,ll,j,k)
     gc = gc + el_3_4 * fi3
     wc = wc + el_4_4 * w_gpu(4,ll,j,k)
     gc = gc + el_4_4 * fi4
     wc = wc + el_5_4 * w_gpu(5,ll,j,k)
     gc = gc + el_5_4 * fi5
     gplus4_2  = 0.5_mykind * (gc + evmax4 * wc)
     gminus4_2 = gc - gplus4_2

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_5 * w_gpu(1,ll,j,k)
     gc = gc + el_1_5 * fi1
     wc = wc + el_2_5 * w_gpu(2,ll,j,k)
     gc = gc + el_2_5 * fi2
     wc = wc + el_3_5 * w_gpu(3,ll,j,k)
     gc = gc + el_3_5 * fi3
     wc = wc + el_4_5 * w_gpu(4,ll,j,k)
     gc = gc + el_4_5 * fi4
     wc = wc + el_5_5 * w_gpu(5,ll,j,k)
     gc = gc + el_5_5 * fi5
     gplus5_2  = 0.5_mykind * (gc + evmax5 * wc)
     gminus5_2 = gc - gplus5_2
!---------------------------------------------
     l  = 3
     ll = i + l - iorder/2
!
     rho  = w_gpu(1,ll,j,k)
     rhou = w_gpu(2,ll,j,k)
     ri   = 1._mykind/rho
     uu   = rhou*ri
     pp   = rho*temperature_gpu(ll,j,k)
     fi1  =       w_gpu(2,ll,j,k)
     fi2  = uu *  w_gpu(2,ll,j,k)  + pp
     fi3  = uu *  w_gpu(3,ll,j,k)
     fi4  = uu *  w_gpu(4,ll,j,k)
     fi5  = uu * (w_gpu(5,ll,j,k)  + pp)
!
     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_1 * w_gpu(1,ll,j,k)
     gc = gc + el_1_1 * fi1
     wc = wc + el_2_1 * w_gpu(2,ll,j,k)
     gc = gc + el_2_1 * fi2
     wc = wc + el_3_1 * w_gpu(3,ll,j,k)
     gc = gc + el_3_1 * fi3
     wc = wc + el_4_1 * w_gpu(4,ll,j,k)
     gc = gc + el_4_1 * fi4
     wc = wc + el_5_1 * w_gpu(5,ll,j,k)
     gc = gc + el_5_1 * fi5
     gplus1_3  = 0.5_mykind * (gc + evmax1 * wc)
     gminus1_3 = gc - gplus1_3

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_2 * w_gpu(1,ll,j,k)
     gc = gc + el_1_2 * fi1
     wc = wc + el_2_2 * w_gpu(2,ll,j,k)
     gc = gc + el_2_2 * fi2
     wc = wc + el_3_2 * w_gpu(3,ll,j,k)
     gc = gc + el_3_2 * fi3
     wc = wc + el_4_2 * w_gpu(4,ll,j,k)
     gc = gc + el_4_2 * fi4
     wc = wc + el_5_2 * w_gpu(5,ll,j,k)
     gc = gc + el_5_2 * fi5
     gplus2_3  = 0.5_mykind * (gc + evmax2 * wc)
     gminus2_3 = gc - gplus2_3

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_3 * w_gpu(1,ll,j,k)
     gc = gc + el_1_3 * fi1
     wc = wc + el_2_3 * w_gpu(2,ll,j,k)
     gc = gc + el_2_3 * fi2
     wc = wc + el_3_3 * w_gpu(3,ll,j,k)
     gc = gc + el_3_3 * fi3
     wc = wc + el_4_3 * w_gpu(4,ll,j,k)
     gc = gc + el_4_3 * fi4
     wc = wc + el_5_3 * w_gpu(5,ll,j,k)
     gc = gc + el_5_3 * fi5
     gplus3_3  = 0.5_mykind * (gc + evmax3 * wc)
     gminus3_3 = gc - gplus3_3

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_4 * w_gpu(1,ll,j,k)
     gc = gc + el_1_4 * fi1
     wc = wc + el_2_4 * w_gpu(2,ll,j,k)
     gc = gc + el_2_4 * fi2
     wc = wc + el_3_4 * w_gpu(3,ll,j,k)
     gc = gc + el_3_4 * fi3
     wc = wc + el_4_4 * w_gpu(4,ll,j,k)
     gc = gc + el_4_4 * fi4
     wc = wc + el_5_4 * w_gpu(5,ll,j,k)
     gc = gc + el_5_4 * fi5
     gplus4_3  = 0.5_mykind * (gc + evmax4 * wc)
     gminus4_3 = gc - gplus4_3

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_5 * w_gpu(1,ll,j,k)
     gc = gc + el_1_5 * fi1
     wc = wc + el_2_5 * w_gpu(2,ll,j,k)
     gc = gc + el_2_5 * fi2
     wc = wc + el_3_5 * w_gpu(3,ll,j,k)
     gc = gc + el_3_5 * fi3
     wc = wc + el_4_5 * w_gpu(4,ll,j,k)
     gc = gc + el_4_5 * fi4
     wc = wc + el_5_5 * w_gpu(5,ll,j,k)
     gc = gc + el_5_5 * fi5
     gplus5_3  = 0.5_mykind * (gc + evmax5 * wc)
     gminus5_3 = gc - gplus5_3
!---------------------------------------------
     l  = 4
     ll = i + l - iorder/2
!
     rho  = w_gpu(1,ll,j,k)
     rhou = w_gpu(2,ll,j,k)
     ri   = 1._mykind/rho
     uu   = rhou*ri
     pp   = rho*temperature_gpu(ll,j,k)
     fi1  =       w_gpu(2,ll,j,k)
     fi2  = uu *  w_gpu(2,ll,j,k)  + pp
     fi3  = uu *  w_gpu(3,ll,j,k)
     fi4  = uu *  w_gpu(4,ll,j,k)
     fi5  = uu * (w_gpu(5,ll,j,k)  + pp)
!
     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_1 * w_gpu(1,ll,j,k)
     gc = gc + el_1_1 * fi1
     wc = wc + el_2_1 * w_gpu(2,ll,j,k)
     gc = gc + el_2_1 * fi2
     wc = wc + el_3_1 * w_gpu(3,ll,j,k)
     gc = gc + el_3_1 * fi3
     wc = wc + el_4_1 * w_gpu(4,ll,j,k)
     gc = gc + el_4_1 * fi4
     wc = wc + el_5_1 * w_gpu(5,ll,j,k)
     gc = gc + el_5_1 * fi5
     gplus1_4  = 0.5_mykind * (gc + evmax1 * wc)
     gminus1_4 = gc - gplus1_4

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_2 * w_gpu(1,ll,j,k)
     gc = gc + el_1_2 * fi1
     wc = wc + el_2_2 * w_gpu(2,ll,j,k)
     gc = gc + el_2_2 * fi2
     wc = wc + el_3_2 * w_gpu(3,ll,j,k)
     gc = gc + el_3_2 * fi3
     wc = wc + el_4_2 * w_gpu(4,ll,j,k)
     gc = gc + el_4_2 * fi4
     wc = wc + el_5_2 * w_gpu(5,ll,j,k)
     gc = gc + el_5_2 * fi5
     gplus2_4  = 0.5_mykind * (gc + evmax2 * wc)
     gminus2_4 = gc - gplus2_4

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_3 * w_gpu(1,ll,j,k)
     gc = gc + el_1_3 * fi1
     wc = wc + el_2_3 * w_gpu(2,ll,j,k)
     gc = gc + el_2_3 * fi2
     wc = wc + el_3_3 * w_gpu(3,ll,j,k)
     gc = gc + el_3_3 * fi3
     wc = wc + el_4_3 * w_gpu(4,ll,j,k)
     gc = gc + el_4_3 * fi4
     wc = wc + el_5_3 * w_gpu(5,ll,j,k)
     gc = gc + el_5_3 * fi5
     gplus3_4  = 0.5_mykind * (gc + evmax3 * wc)
     gminus3_4 = gc - gplus3_4

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_4 * w_gpu(1,ll,j,k)
     gc = gc + el_1_4 * fi1
     wc = wc + el_2_4 * w_gpu(2,ll,j,k)
     gc = gc + el_2_4 * fi2
     wc = wc + el_3_4 * w_gpu(3,ll,j,k)
     gc = gc + el_3_4 * fi3
     wc = wc + el_4_4 * w_gpu(4,ll,j,k)
     gc = gc + el_4_4 * fi4
     wc = wc + el_5_4 * w_gpu(5,ll,j,k)
     gc = gc + el_5_4 * fi5
     gplus4_4  = 0.5_mykind * (gc + evmax4 * wc)
     gminus4_4 = gc - gplus4_4

     wc = 0._mykind
     gc = 0._mykind
     wc = wc + el_1_5 * w_gpu(1,ll,j,k)
     gc = gc + el_1_5 * fi1
     wc = wc + el_2_5 * w_gpu(2,ll,j,k)
     gc = gc + el_2_5 * fi2
     wc = wc + el_3_5 * w_gpu(3,ll,j,k)
     gc = gc + el_3_5 * fi3
     wc = wc + el_4_5 * w_gpu(4,ll,j,k)
     gc = gc + el_4_5 * fi4
     wc = wc + el_5_5 * w_gpu(5,ll,j,k)
     gc = gc + el_5_5 * fi5
     gplus5_4  = 0.5_mykind * (gc + evmax5 * wc)
     gminus5_4 = gc - gplus5_4
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!
!    Reconstruction of the '+' and '-' fluxes
!
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
     vp_weno1_1 = gplus1_1
     vp_weno2_1 = gplus2_1
     vp_weno3_1 = gplus3_1
     vp_weno4_1 = gplus4_1
     vp_weno5_1 = gplus5_1
     vp_weno1_2 = gplus1_2
     vp_weno2_2 = gplus2_2
     vp_weno3_2 = gplus3_2
     vp_weno4_2 = gplus4_2
     vp_weno5_2 = gplus5_2
     vp_weno1_3 = gplus1_3
     vp_weno2_3 = gplus2_3
     vp_weno3_3 = gplus3_3
     vp_weno4_3 = gplus4_3
     vp_weno5_3 = gplus5_3
     vp_weno1_4 = gplus1_4
     vp_weno2_4 = gplus2_4
     vp_weno3_4 = gplus3_4
     vp_weno4_4 = gplus4_4
     vp_weno5_4 = gplus5_4
     vm_weno1_1 = gminus1_1
     vm_weno2_1 = gminus2_1
     vm_weno3_1 = gminus3_1
     vm_weno4_1 = gminus4_1
     vm_weno5_1 = gminus5_1
     vm_weno1_2 = gminus1_2
     vm_weno2_2 = gminus2_2
     vm_weno3_2 = gminus3_2
     vm_weno4_2 = gminus4_2
     vm_weno5_2 = gminus5_2
     vm_weno1_3 = gminus1_3
     vm_weno2_3 = gminus2_3
     vm_weno3_3 = gminus3_3
     vm_weno4_3 = gminus4_3
     vm_weno5_3 = gminus5_3
     vm_weno1_4 = gminus1_4
     vm_weno2_4 = gminus2_4
     vm_weno3_4 = gminus3_4
     vm_weno4_4 = gminus4_4
     vm_weno5_4 = gminus5_4
!
!-------------------------------------------------------------------------------------------------------------
     dwe1   = 2._mykind/3._mykind
     dwe0   = 1._mykind/3._mykind
!
! m=1
     betap0  = (vp_weno1_2-vp_weno1_1)**2
     betap1  = (vp_weno1_3-vp_weno1_2)**2
     betam0  = (vm_weno1_4-vm_weno1_3)**2
     betam1  = (vm_weno1_3-vm_weno1_2)**2
!
     sump = 0._mykind
     summ = 0._mykind
     alfp0 = dwe0/(eps+betap0)**2
     alfm0 = dwe0/(eps+betam0)**2
     sump = sump + alfp0
     summ = summ + alfm0
     alfp1 = dwe1/(eps+betap1)**2
     alfm1 = dwe1/(eps+betam1)**2
     sump = sump + alfp1
     summ = summ + alfm1
     omp0 = alfp0/sump
     omm0 = alfm0/summ
     omp1 = alfp1/sump
     omm1 = alfm1/summ
!
     vminus1 = omp0 *(-vp_weno1_1+3._mykind*vp_weno1_2) + omp1 *( vp_weno1_2+ vp_weno1_3)
     vplus1  = omm0 *(-vm_weno1_4+3._mykind*vm_weno1_3) + omm1 *( vm_weno1_2+ vm_weno1_3)
! m=2
     betap0  = (vp_weno2_2-vp_weno2_1)**2
     betap1  = (vp_weno2_3-vp_weno2_2)**2
     betam0  = (vm_weno2_4-vm_weno2_3)**2
     betam1  = (vm_weno2_3-vm_weno2_2)**2
!
     sump = 0._mykind
     summ = 0._mykind
     alfp0 = dwe0/(eps+betap0)**2
     alfm0 = dwe0/(eps+betam0)**2
     sump = sump + alfp0
     summ = summ + alfm0
     alfp1 = dwe1/(eps+betap1)**2
     alfm1 = dwe1/(eps+betam1)**2
     sump = sump + alfp1
     summ = summ + alfm1
     omp0 = alfp0/sump
     omm0 = alfm0/summ
     omp1 = alfp1/sump
     omm1 = alfm1/summ
!
     vminus2 = omp0 *(-vp_weno2_1+3._mykind*vp_weno2_2) + omp1 *( vp_weno2_2+ vp_weno2_3)
     vplus2  = omm0 *(-vm_weno2_4+3._mykind*vm_weno2_3) + omm1 *( vm_weno2_2+ vm_weno2_3)
! m=3
     betap0  = (vp_weno3_2-vp_weno3_1)**2
     betap1  = (vp_weno3_3-vp_weno3_2)**2
     betam0  = (vm_weno3_4-vm_weno3_3)**2
     betam1  = (vm_weno3_3-vm_weno3_2)**2
!
     sump = 0._mykind
     summ = 0._mykind
     alfp0 = dwe0/(eps+betap0)**2
     alfm0 = dwe0/(eps+betam0)**2
     sump = sump + alfp0
     summ = summ + alfm0
     alfp1 = dwe1/(eps+betap1)**2
     alfm1 = dwe1/(eps+betam1)**2
     sump = sump + alfp1
     summ = summ + alfm1
     omp0 = alfp0/sump
     omm0 = alfm0/summ
     omp1 = alfp1/sump
     omm1 = alfm1/summ
!
     vminus3 = omp0 *(-vp_weno3_1+3._mykind*vp_weno3_2) + omp1 *( vp_weno3_2+ vp_weno3_3)
     vplus3  = omm0 *(-vm_weno3_4+3._mykind*vm_weno3_3) + omm1 *( vm_weno3_2+ vm_weno3_3)
! m=4
     betap0  = (vp_weno4_2-vp_weno4_1)**2
     betap1  = (vp_weno4_3-vp_weno4_2)**2
     betam0  = (vm_weno4_4-vm_weno4_3)**2
     betam1  = (vm_weno4_3-vm_weno4_2)**2
!
     sump = 0._mykind
     summ = 0._mykind
     alfp0 = dwe0/(eps+betap0)**2
     alfm0 = dwe0/(eps+betam0)**2
     sump = sump + alfp0
     summ = summ + alfm0
     alfp1 = dwe1/(eps+betap1)**2
     alfm1 = dwe1/(eps+betam1)**2
     sump = sump + alfp1
     summ = summ + alfm1
     omp0 = alfp0/sump
     omm0 = alfm0/summ
     omp1 = alfp1/sump
     omm1 = alfm1/summ
!
     vminus4 = omp0 *(-vp_weno4_1+3._mykind*vp_weno4_2) + omp1 *( vp_weno4_2+ vp_weno4_3)
     vplus4  = omm0 *(-vm_weno4_4+3._mykind*vm_weno4_3) + omm1 *( vm_weno4_2+ vm_weno4_3)
! m=5
     betap0  = (vp_weno5_2-vp_weno5_1)**2
     betap1  = (vp_weno5_3-vp_weno5_2)**2
     betam0  = (vm_weno5_4-vm_weno5_3)**2
     betam1  = (vm_weno5_3-vm_weno5_2)**2
!
     sump = 0._mykind
     summ = 0._mykind
     alfp0 = dwe0/(eps+betap0)**2
     alfm0 = dwe0/(eps+betam0)**2
     sump = sump + alfp0
     summ = summ + alfm0
     alfp1 = dwe1/(eps+betap1)**2
     alfm1 = dwe1/(eps+betam1)**2
     sump = sump + alfp1
     summ = summ + alfm1
     omp0 = alfp0/sump
     omm0 = alfm0/summ
     omp1 = alfp1/sump
     omm1 = alfm1/summ
!
     vminus5 = omp0 *(-vp_weno5_1+3._mykind*vp_weno5_2) + omp1 *( vp_weno5_2+ vp_weno5_3)
     vplus5  = omm0 *(-vm_weno5_4+3._mykind*vm_weno5_3) + omm1 *( vm_weno5_2+ vm_weno5_3)
!
     vminus1 = 0.5_mykind*vminus1
     vplus1  = 0.5_mykind*vplus1
     vminus2 = 0.5_mykind*vminus2
     vplus2  = 0.5_mykind*vplus2
     vminus3 = 0.5_mykind*vminus3
     vplus3  = 0.5_mykind*vplus3
     vminus4 = 0.5_mykind*vminus4
     vplus4  = 0.5_mykind*vplus4
     vminus5 = 0.5_mykind*vminus5
     vplus5  = 0.5_mykind*vplus5
!-------------------------------------------------------------------------------------------------------------
     gl1 = vminus1 
     gr1 = vplus1
     gl2 = vminus2 
     gr2 = vplus2
     gl3 = vminus3 
     gr3 = vplus3
     gl4 = vminus4 
     gr4 = vplus4
     gl5 = vminus5 
     gr5 = vplus5
!-------------------------------------------------------------------------------------------------------------
!
     ghat1 = gl1 + gr1 ! char._mykind flux
     ghat2 = gl2 + gr2 ! char._mykind flux
     ghat3 = gl3 + gr3 ! char._mykind flux
     ghat4 = gl4 + gr4 ! char._mykind flux
     ghat5 = gl5 + gr5 ! char._mykind flux
!
!   Return to conservative fluxes
     fhat_gpu(1,i,j,k) = 0._mykind
     fhat_gpu(1,i,j,k) = fhat_gpu(1,i,j,k) + er_1_1 * ghat1
     fhat_gpu(1,i,j,k) = fhat_gpu(1,i,j,k) + er_2_1 * ghat2
     fhat_gpu(1,i,j,k) = fhat_gpu(1,i,j,k) + er_3_1 * ghat3
     fhat_gpu(1,i,j,k) = fhat_gpu(1,i,j,k) + er_4_1 * ghat4
     fhat_gpu(1,i,j,k) = fhat_gpu(1,i,j,k) + er_5_1 * ghat5
     fhat_gpu(2,i,j,k) = 0._mykind
     fhat_gpu(2,i,j,k) = fhat_gpu(2,i,j,k) + er_1_2 * ghat1
     fhat_gpu(2,i,j,k) = fhat_gpu(2,i,j,k) + er_2_2 * ghat2
     fhat_gpu(2,i,j,k) = fhat_gpu(2,i,j,k) + er_3_2 * ghat3
     fhat_gpu(2,i,j,k) = fhat_gpu(2,i,j,k) + er_4_2 * ghat4
     fhat_gpu(2,i,j,k) = fhat_gpu(2,i,j,k) + er_5_2 * ghat5
     fhat_gpu(3,i,j,k) = 0._mykind
     fhat_gpu(3,i,j,k) = fhat_gpu(3,i,j,k) + er_1_3 * ghat1
     fhat_gpu(3,i,j,k) = fhat_gpu(3,i,j,k) + er_2_3 * ghat2
     fhat_gpu(3,i,j,k) = fhat_gpu(3,i,j,k) + er_3_3 * ghat3
     fhat_gpu(3,i,j,k) = fhat_gpu(3,i,j,k) + er_4_3 * ghat4
     fhat_gpu(3,i,j,k) = fhat_gpu(3,i,j,k) + er_5_3 * ghat5
     fhat_gpu(4,i,j,k) = 0._mykind
     fhat_gpu(4,i,j,k) = fhat_gpu(4,i,j,k) + er_1_4 * ghat1
     fhat_gpu(4,i,j,k) = fhat_gpu(4,i,j,k) + er_2_4 * ghat2
     fhat_gpu(4,i,j,k) = fhat_gpu(4,i,j,k) + er_3_4 * ghat3
     fhat_gpu(4,i,j,k) = fhat_gpu(4,i,j,k) + er_4_4 * ghat4
     fhat_gpu(4,i,j,k) = fhat_gpu(4,i,j,k) + er_5_4 * ghat5
     fhat_gpu(5,i,j,k) = 0._mykind
     fhat_gpu(5,i,j,k) = fhat_gpu(5,i,j,k) + er_1_5 * ghat1
     fhat_gpu(5,i,j,k) = fhat_gpu(5,i,j,k) + er_2_5 * ghat2
     fhat_gpu(5,i,j,k) = fhat_gpu(5,i,j,k) + er_3_5 * ghat3
     fhat_gpu(5,i,j,k) = fhat_gpu(5,i,j,k) + er_4_5 * ghat4
     fhat_gpu(5,i,j,k) = fhat_gpu(5,i,j,k) + er_5_5 * ghat5
!
    endif
   enddo ! end of loop on the cell faces
!
!  Update net flux 
   if (iend == endi) then
    do i=1,endi ! loop on the inner nodes
     im = i-1
     do m=1,5
      df = (fhat_gpu(m,i,j,k)-fhat_gpu(m,im,j,k))*dcsidx_gpu(i)
      fl_gpu(m,i,j,k) = fl_gpu(m,i,j,k) + df
     enddo
    enddo
   endif
  enddo ! end of j-loop
 enddo ! end of k-loop
!
#ifdef USE_CUDA
 if (iend == (iorder/2-1)) then
   !@cuf iercuda=cudaDeviceSynchronize()
 endif
#endif
!
end subroutine euler_i
