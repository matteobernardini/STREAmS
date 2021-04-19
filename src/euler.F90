#include "cuda_definitions.h"

module mod_euler

contains

subroutine euler_i(istart,iend,num_call,last_call,i1,i2)
!
! Evaluation of convective terms in the x direction
!
 use mod_streams
 implicit none
!
 integer :: istart,iend,num_call,last_call,i1,i2
 integer :: lmax,ng2
 real(mykind) :: tt,st,et
#ifdef USE_CUDA
 type(dim3) :: grid, tBlock
 integer :: i, j, k, iv, n_th_x, n_th_y
#else
 real(mykind),dimension(5,1-ng:nx+ng) :: ui,fi,ev
 real(mykind),dimension(5) :: evmax
 real(mykind),dimension(5,5) :: el,er
 real(mykind),dimension(5) :: gr,gl,ghat
 real(mykind),dimension(5,2*iweno) :: gplus,gminus
 integer,dimension(1-ng:nx+ng) :: ishk
!
 real(mykind),dimension(1-ng:nx+ng) :: den
 real(mykind),dimension(6,1-ng:nx+ng) :: u,v
 real(mykind),dimension(6,ng,1-ng:nx+ng) :: uv
 real(mykind),dimension(5,nx) :: df
 real(mykind),dimension(5) :: fh
 real(mykind),dimension(6) :: ft,uvs
!
 integer :: i,j,k,l,m,n
 integer :: im,imax,imin,ip,ll,mm
 real(mykind) :: b1,b2,c,cc,ci,ducm,ent,gc,gpr,h,hp,pp,qq,qqp
 real(mykind) :: r,rho,rhoe,rhom,rhou,rhov,rhow,ri,rip,rp1,up
 real(mykind) :: uu,vp,vv,wc,wp,ww
#endif
!----------------------------------------------------------------------------------------
!
 ng2  = 2*iweno
 lmax = iorder/2 ! max stencil width
!
#ifdef USE_CUDA
#ifdef CUDA_ASYNC
  select case (num_call)
  case (1) ! Inner nodes
   !$cuf kernel do(3) <<<*,*,stream=stream1>>>
   do k=1,nz
    do i=1,nx
     do j=1,ny
      do iv=1,nv
       wv_trans_gpu(j,i,k,iv) = wv_gpu(i,j,k,iv)
      enddo
     enddo
    enddo
   enddo
   !$cuf kernel do(3) <<<*,*,stream=stream1>>>
   do k=1,nz
    do i=1,nx
     do j=1,ny
       temperature_trans_gpu(j,i,k) = temperature_gpu(i,j,k)
     enddo
    enddo
   enddo
  case(2) ! First nodes 
   !$cuf kernel do(3) <<<*,*,stream=stream1>>>
   do k=1,nz
    do i=1-ng,0
     do j=1,ny
      do iv=1,nv
       wv_trans_gpu(j,i,k,iv) = wv_gpu(i,j,k,iv)
      enddo
     enddo
    enddo
   enddo
   !$cuf kernel do(3) <<<*,*,stream=stream1>>>
   do k=1,nz
    do i=1-ng,0
     do j=1,ny
       temperature_trans_gpu(j,i,k) = temperature_gpu(i,j,k)
     enddo
    enddo
   enddo
   case(3)
   !$cuf kernel do(3) <<<*,*,stream=stream1>>>
   do k=1,nz
    do i=nx+1,nx+ng
     do j=1,ny
      do iv=1,nv
       wv_trans_gpu(j,i,k,iv) = wv_gpu(i,j,k,iv)
      enddo
     enddo
    enddo
   enddo
  !$cuf kernel do(3) <<<*,*,stream=stream1>>>
   do k=1,nz
    do i=nx+1,nx+ng
     do j=1,ny
       temperature_trans_gpu(j,i,k) = temperature_gpu(i,j,k)
     enddo
    enddo
   enddo
  end select
#else
!$cuf kernel do(3) <<<*,*,stream=stream1>>>
  do k=1,nz
   do i=1-ng,nx+ng
    do j=1,ny
     do iv=1,nv
      wv_trans_gpu(j,i,k,iv) = wv_gpu(i,j,k,iv)
     enddo
    enddo
   enddo
  enddo
!$cuf kernel do(3) <<<*,*,stream=stream1>>>
  do k=1,nz
   do i=1-ng,nx+ng
    do j=1,ny
      temperature_trans_gpu(j,i,k) = temperature_gpu(i,j,k)
    enddo
   enddo
  enddo
#endif

  if (tresduc<1._mykind) then
   n_th_x = EULERWENO_THREADS_X
   n_th_y = EULERWENO_THREADS_Y
   tBlock = dim3(n_th_x,n_th_y,1)
   grid = dim3(ceiling(real(ny)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
   call euler_i_kernel<<<grid, tBlock, 0, stream1>>>(nv, nx, ny, nz, ng, ng2, gamma, gm1, &
     iweno, istart, iend, num_call, last_call, i1, i2, lmax, iorder, iflow, &
     w_gpu, temperature_gpu, ducros_gpu, &
     fhat_trans_gpu, temperature_trans_gpu, fl_trans_gpu, & 
     fhat_gpu, fl_gpu, dcoe_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, gplus_x, gminus_x, wv_trans_gpu )
  else
   n_th_x = EULERCENTRAL_THREADS_X
   n_th_y = EULERCENTRAL_THREADS_Y
   tBlock = dim3(n_th_x,n_th_y,1)
   grid = dim3(ceiling(real(ny)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
   call euler_i_kernel_central<<<grid, tBlock, 0, stream1>>>(nv, nx, ny, nz, ng, &
    istart, iend, num_call, last_call, i1, i2, lmax, &
    fhat_trans_gpu, temperature_trans_gpu, fl_trans_gpu, dcoe_gpu, dcsidx_gpu, wv_trans_gpu)
  endif

 !iercuda = cudaGetLastError()
 !if (iercuda /= cudaSuccess) then
 ! call fail_input("CUDA ERROR! Try to reduce the number of Euler threads in cuda_definitions.h: "//cudaGetErrorString(iercuda))
 !endif
!
 if (num_call==last_call) then
  !$cuf kernel do(3) <<<*,*,stream=stream1>>>
  do k=1,nz
   do j=1,ny
    do i=i1,i2
     do iv=1,nv
      fl_gpu(i,j,k,iv) = fl_trans_gpu(j,i,k,iv)
     enddo
    enddo
   enddo
  enddo
  !@cuf iercuda=cudaDeviceSynchronize()
 endif
!
#else
!
 imin = istart - lmax
 imax = iend   + lmax
 do k=1,nz ! loop in the z-direction
  do j=1,ny ! loop in the y-direction
!  Sweep in the x-direction
   do i=imin,imax ! loop on all the nodes
    do m=1,5
     ui(m,i) = w_gpu(i,j,k,m)
    enddo
    rho  =  w_gpu(i,j,k,1)
    rhou =  w_gpu(i,j,k,2)
    rhov =  w_gpu(i,j,k,3)
    rhow =  w_gpu(i,j,k,4)
    rhoe =  w_gpu(i,j,k,5)
    ri   =  1._mykind/rho
    uu   =  rhou*ri
    vv   =  rhov*ri
    ww   =  rhow*ri
    qq   =  0.5_mykind*(uu*uu+vv*vv+ww*ww)
    pp   =  rho*temperature_gpu(i,j,k)
    ent  =  (rhoe+pp)*ri
    do m=1,5
     u(m,i) = uu
    enddo
    u(6,i) = 1._mykind
    v(1,i) = 1._mykind
    v(2,i) = uu
    v(3,i) = vv
    v(4,i) = ww
    v(5,i) = ent
    v(6,i) = pp
    den(i) = rho
    gpr = gamma * pp * ri
    c  =  sqrt (gpr)
    fi(1,i) =       ui(2,i)
    fi(2,i) = uu *  ui(2,i)  + pp
    fi(3,i) = uu *  ui(3,i)
    fi(4,i) = uu *  ui(4,i)
    fi(5,i) = uu * (ui(5,i)  + pp)
    ev(1,i) = abs(uu-c)
    ev(2,i) = abs(uu)
    ev(3,i) = abs(uu+c)
    ev(4,i) = ev(2,i)
    ev(5,i) = ev(2,i)
   enddo ! end of loop on the nodes
!
!  Compute u*v averages
!
   do i=imin,iend
    do l=1,lmax
     rhom = den(i)+den(i+l)
     do m=1,5
      uv(m,l,i) = (u(m,i)+u(m,i+l))*(v(m,i)+v(m,i+l))*rhom
     enddo
     uv(6,l,i) = (u(6,i)+u(6,i+l))*(v(6,i)+v(6,i+l))
    enddo
   enddo
!
!  Smoothness properties
   ishk = 0
   do i=istart-1,iend
    if (any(ducros_gpu(i-iweno+1:i+iweno,j,k))) then
     ishk(i) = 1
    endif
   enddo
!
!  Evaluation of numerical fluxes
!
   do i=istart-1,iend ! i is the index of intermediate nodes
!
    if (ishk(i)==0) then
     ft = 0._mykind
     do l=1,lmax
      uvs = 0._mykind
      do n=0,l-1
       uvs = uvs + uv(:,l,i-n)
      enddo
      ft = ft + dcoe_gpu(l,lmax)*uvs
     enddo
     do m=1,5
      fh(m) = 0.25_mykind*ft(m)
     enddo
     fh(2) = fh(2) + 0.5_mykind*ft(6)
!
     fhat_gpu(i,j,k,1:5) = fh(:)
!
    else

     ip = i + 1
!
!    Compute Roe average
!
!    Left state (node i)
     ri        =  1._mykind/ui(1,i)
     uu        =  ui(2,i) * ri
     vv        =  ui(3,i) * ri
     ww        =  ui(4,i) * ri
     qq        =  0.5_mykind * (uu*uu  +vv*vv + ww*ww)
!    pp        =  gm1 * (ui(5,i) - ui(1,i) * qq)
!    h         =  (ui(5,i)  + pp) * ri
     h         =  gamma*ui(5,i)*ri-gm1*qq
!    Right state (node i+1)
     rip       =  1._mykind/ui(1,ip)
     up        =  ui(2,ip) * rip
     vp        =  ui(3,ip) * rip
     wp        =  ui(4,ip) * rip
     qqp       =  0.5_mykind * (up*up  +vp*vp +wp*wp)
!    ppp       =  gm1 * (ui(5,ip) - ui(1,ip) * qqp)
!    hp        =  (ui(5,ip)  + ppp) * rip
     hp        =  gamma*ui(5,ip)*rip-gm1*qqp
!    average state
     r         =  ui(1,ip)*ri
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
     el(1,1)   =   0.5_mykind * (b1     + uu * ci)
     el(2,1)   =  -0.5_mykind * (b2 * uu +     ci)
     el(3,1)   =  -0.5_mykind * (b2 * vv         )
     el(4,1)   =  -0.5_mykind * (b2 * ww         )
     el(5,1)   =   0.5_mykind * b2
     el(1,2)   =   1._mykind - b1
     el(2,2)   =   b2*uu
     el(3,2)   =   b2*vv
     el(4,2)   =   b2*ww
     el(5,2)   =  -b2
     el(1,3)   =   0.5_mykind * (b1     - uu * ci)
     el(2,3)   =  -0.5_mykind * (b2 * uu -     ci)
     el(3,3)   =  -0.5_mykind * (b2 * vv         )
     el(4,3)   =  -0.5_mykind * (b2 * ww         )
     el(5,3)   =   0.5_mykind * b2
     el(1,4)   =  -vv
     el(2,4)   =   0._mykind
     el(3,4)   =   1._mykind
     el(4,4)   =   0._mykind
     el(5,4)   =   0._mykind
     el(1,5)   =  -ww
     el(2,5)   =   0._mykind
     el(3,5)   =   0._mykind
     el(4,5)   =   1._mykind
     el(5,5)   =   0._mykind
!
!    right eigenvectors matrix (at Roe state)
!
     er(1,1)   =  1._mykind
     er(2,1)   =  1._mykind
     er(3,1)   =  1._mykind
     er(4,1)   =  0._mykind
     er(5,1)   =  0._mykind
     er(1,2)   =  uu -  c
     er(2,2)   =  uu
     er(3,2)   =  uu +  c
     er(4,2)   =  0._mykind
     er(5,2)   =  0._mykind
     er(1,3)   =  vv
     er(2,3)   =  vv
     er(3,3)   =  vv
     er(4,3)   =  1._mykind
     er(5,3)   =  0._mykind
     er(1,4)   =  ww
     er(2,4)   =  ww
     er(3,4)   =  ww
     er(4,4)   =  0._mykind
     er(5,4)   =  1._mykind
     er(1,5)   =  h  - uu * c
     er(2,5)   =  qq
     er(3,5)   =  h  + uu * c
     er(4,5)   =  vv
     er(5,5)   =  ww

     do  m=1,5 ! loop on characteristic fields
      evmax(m) = -1._mykind
     enddo
     do l=1,ng2 ! LLF
      ll = i + l - iweno
      do  m=1,5
       evmax(m) = max(ev(m,ll),evmax(m))
      enddo
     enddo
     do l=1,ng2 ! loop over the stencil centered at face i
      ll = i + l - iweno
      do m=1,5
       wc = 0._mykind
       gc = 0._mykind
       do mm=1,5
        wc = wc + el(mm,m) * ui(mm,ll)
        gc = gc + el(mm,m) * fi(mm,ll)
       enddo
       gplus (m,l)   = 0.5_mykind * (gc + evmax(m) * wc)
       gminus(m,l)   = gc - gplus(m,l)
      enddo
     enddo
!
!    Reconstruction of the '+' and '-' fluxes
!
     call wenorec(5,gplus,gminus,gl,gr,iweno)
!
     do m=1,5
      ghat(m) = gl(m) + gr(m) ! char. flux
     enddo
!
!   Return to conservative fluxes
     do m=1,5
      fhat_gpu(i,j,k,m) = 0._mykind
      do mm=1,5
       fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * ghat(mm)
      enddo
     enddo
!
    endif
   enddo ! end of loop on the cell faces
!
!  Update net flux
   do i=istart,iend ! loop on the inner nodes
    do m=1,5
     fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + (fhat_gpu(i,j,k,m)-fhat_gpu(i-1,j,k,m))*dcsidx_gpu(i)
    enddo
   enddo
  enddo ! end of j-loop
 enddo ! end of k-loop
#endif
!
#ifdef USE_CUDA
 if (num_call/=1) then
   !@cuf iercuda=cudaDeviceSynchronize()
 endif
#endif
!
end subroutine euler_i

subroutine euler_j(jstart,jend)
!
! Evaluation of convective terms in the y direction
!
 use mod_streams
 implicit none
!
 integer :: jstart,jend
 integer :: lmax,ng2
!
#ifdef USE_CUDA
 type(dim3) :: grid, tBlock
 real(mykind) :: tt,st,et
 integer :: n_th_x, n_th_y
#else
 real(mykind),dimension(5,1-ng:ny+ng) :: uj,fj,ev
 real(mykind),dimension(5) :: evmax
 real(mykind),dimension(5,5) :: el,er
 real(mykind),dimension(5) :: gr,gl,ghat
 real(mykind),dimension(5,2*ng) :: gplus,gminus
 integer,dimension(1-ng:ny+ng) :: ishk
!
 real(mykind),dimension(1-ng:ny+ng) :: den
 real(mykind),dimension(6,1-ng:ny+ng) :: u,v
 real(mykind),dimension(6,ng,1-ng:ny+ng) :: uv
 real(mykind),dimension(5,ny) :: df
 real(mykind),dimension(5) :: fh
 real(mykind),dimension(6) :: ft,uvs
!
 integer :: i,j,k,l,m,n
 integer :: jm,jmax,jmin,jp,ll,mm
 real(mykind) :: b1,b2,c,cc,ci,ducm,ent,gc,gpr,h,hp,pp,qq,qqp
 real(mykind) :: r,rho,rhoe,rhom,rhou,rhov,rhow,ri,rip,rp1,up
 real(mykind) :: uu,vp,vv,wc,wp,ww
#endif
!----------------------------------------------------------------------------------------
!
 ng2  = 2*iweno
 lmax = iorder/2 ! max stencil width
!
#ifdef USE_CUDA
    if (tresduc<1._mykind) then
     n_th_x = EULERWENO_THREADS_X
     n_th_y = EULERWENO_THREADS_Y
     tBlock = dim3(n_th_x,n_th_y,1)
     grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
     call euler_j_kernel<<<grid, tBlock, 0, stream1>>>(nv, nx, ny, nz, ng, ng2, gamma, gm1, &
         iweno, jstart, jend, lmax, iorder, iflow, & 
         w_gpu, temperature_gpu, ducros_gpu, &
         fhat_gpu, fl_gpu, dcoe_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, gplus_y, gminus_y, wv_gpu )
    else
     n_th_x = EULERCENTRAL_THREADS_X
     n_th_y = EULERCENTRAL_THREADS_Y
     tBlock = dim3(n_th_x,n_th_y,1)
     grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(nz)/tBlock%y),1)
     call euler_j_kernel_central<<<grid, tBlock, 0, stream1>>>(nv, nx, ny, nz, ng, &
         jstart, jend, lmax, iflow, &
         temperature_gpu, fhat_gpu, fl_gpu, dcoe_gpu, detady_gpu, wv_gpu)
    endif
!
#else
!
 jmin = jstart  - lmax
 jmax = jend + lmax
 do k=1,nz ! loop in the z-direction
  do i=1,nx ! loop in the x-direction
!  Sweep in the y-direction
   do j=jmin,jmax ! loop on all nodes needed for updating
    do m=1,5
     uj(m,j) = w_gpu(i,j,k,m)
    enddo
    rho  =  w_gpu(i,j,k,1)
    rhou =  w_gpu(i,j,k,2)
    rhov =  w_gpu(i,j,k,3)
    rhow =  w_gpu(i,j,k,4)
    rhoe =  w_gpu(i,j,k,5)
    ri   =  1._mykind/rho
    uu   =  rhou*ri
    vv   =  rhov*ri
    ww   =  rhow*ri
    qq   =  0.5_mykind*(uu*uu+vv*vv+ww*ww)
    pp   =  rho*temperature_gpu(i,j,k)
    ent  =  (rhoe+pp)*ri
    do m=1,5
     u(m,j) = vv
    enddo
    u(6,j) = 1._mykind
    v(1,j) = 1._mykind
    v(2,j) = uu
    v(3,j) = vv
    v(4,j) = ww
    v(5,j) = ent
    v(6,j) = pp
    den(j) = rho
    gpr = gamma * pp * ri
    c  =  sqrt (gpr)
    fj(1,j) =       uj(3,j)
    fj(2,j) = vv *  uj(2,j)
    fj(3,j) = vv *  uj(3,j) + pp
    fj(4,j) = vv *  uj(4,j)
    fj(5,j) = vv * (uj(5,j)  + pp)
    ev(1,j) = abs(vv-c)
    ev(2,j) = abs(vv)
    ev(3,j) = abs(vv+c)
    ev(4,j) = ev(2,j)
    ev(5,j) = ev(2,j)
   enddo ! end of loop on the nodes
!
!  Compute u*v averages
!
   do j=jmin,jend
    do l=1,lmax
     rhom = den(j)+den(j+l)
     do m=1,5
      uv(m,l,j) = (u(m,j)+u(m,j+l))*(v(m,j)+v(m,j+l))*rhom
     enddo
     uv(6,l,j) = (u(6,j)+u(6,j+l))*(v(6,j)+v(6,j+l))
    enddo
   enddo
!
!  Smoothness properties
   ishk = 0
   do j=jstart-1,jend
    if (any(ducros_gpu(i,j-iweno+1:j+iweno,k))) then
     ishk(j) = 1
    endif
   enddo
!
!  Evaluation of numerical fluxes
!
   do j=jstart-1,jend ! j is the index of intermediate nodes (faces)
!
    if (ishk(j)==0) then
     ft = 0._mykind
     do l=1,lmax
      uvs = 0._mykind
      do n=0,l-1
       uvs = uvs + uv(:,l,j-n)
      enddo
      ft = ft + dcoe_gpu(l,lmax)*uvs
     enddo
     do m=1,5
      fh(m) = 0.25_mykind*ft(m)
     enddo
     if (iflow==0) then
      if (j==0.or.j==ny) then
       fh(:) = 0._mykind
      endif
     endif
     fh(3) = fh(3) + 0.5_mykind*ft(6)
!
     fhat_gpu(i,j,k,1:5) = fh(:)
!
    else
!
     jp = j + 1
!
!   Compute Roe average
!
!    Left state (node j)
     ri        =  1._mykind/uj(1,j)
     uu        =  uj(2,j) * ri
     vv        =  uj(3,j) * ri
     ww        =  uj(4,j) * ri
     qq        =  0.5_mykind * (uu*uu  +vv*vv + ww*ww)
!    pp        =  gm1 * (uj(5,j) - uj(1,j) * qq)
!    h         =  (uj(5,j)  + pp) * ri
     h         =  gamma*uj(5,j)*ri-gm1*qq
!    Right state (node j+1)
     rip       =  1._mykind/uj(1,jp)
     up        =  uj(2,jp) * rip
     vp        =  uj(3,jp) * rip
     wp        =  uj(4,jp) * rip
     qqp       =  0.5_mykind * (up*up  +vp*vp + wp*wp)
!    ppp       =  gm1 * (uj(5,jp) - uj(1,jp) * qqp)
!    hp        =  (uj(5,jp)  + ppp) * rip
     hp        =  gamma*uj(5,jp)*rip-gm1*qqp
!    average state
     r         =  uj(1,jp)*ri
     r         =  sqrt(r)
     rp1       =  1._mykind/(r  +1._mykind)
     uu        =  (r*up  +uu)*rp1
     vv        =  (r*vp  +vv)*rp1
     ww        =  (r*wp  +ww)*rp1
     h         =  (r*hp  +h)*rp1
     qq        =  0.5_mykind * (uu*uu  + vv*vv + ww*ww)
     cc        =  gm1 * (h - qq)
     c         =  sqrt(cc)
     ci        =  1._mykind/c
!
!    left eigenvectors matrix (at Roe state)
!
     b2        =   gm1/cc
     b1        =   b2 * qq
     el(1,1)   =   0.5_mykind * (b1     + vv * ci)
     el(2,1)   =  -0.5_mykind * (b2 * uu         )
     el(3,1)   =  -0.5_mykind * (b2 * vv +     ci)
     el(4,1)   =  -0.5_mykind * (b2 * ww         )
     el(5,1)   =   0.5_mykind * b2
     el(1,2)   =   1._mykind - b1
     el(2,2)   =   b2 * uu
     el(3,2)   =   b2 * vv
     el(4,2)   =   b2 * ww
     el(5,2)   =  -b2
     el(1,3)   =   0.5_mykind * (b1     - vv * ci)
     el(2,3)   =  -0.5_mykind * (b2 * uu         )
     el(3,3)   =  -0.5_mykind * (b2 * vv -     ci)
     el(4,3)   =  -0.5_mykind * (b2 * ww         )
     el(5,3)   =   0.5_mykind * b2
     el(1,4)   =  -ww
     el(2,4)   =   0._mykind
     el(3,4)   =   0._mykind
     el(4,4)   =   1._mykind
     el(5,4)   =   0._mykind
     el(1,5)   =   uu
     el(2,5)   =  -1._mykind
     el(3,5)   =   0._mykind
     el(4,5)   =   0._mykind
     el(5,5)   =   0._mykind
!
!    right eigenvectors matrix (at Roe state)
!
     er(1,1)   =  1._mykind
     er(2,1)   =  1._mykind
     er(3,1)   =  1._mykind
     er(4,1)   =  0._mykind
     er(5,1)   =  0._mykind
     er(1,2)   =  uu
     er(2,2)   =  uu
     er(3,2)   =  uu
     er(4,2)   =  0._mykind
     er(5,2)   = -1._mykind
     er(1,3)   =  vv - c
     er(2,3)   =  vv
     er(3,3)   =  vv + c
     er(4,3)   =  0._mykind
     er(5,3)   =  0._mykind
     er(1,4)   =  ww
     er(2,4)   =  ww
     er(3,4)   =  ww
     er(4,4)   =  1._mykind
     er(5,4)   =  0._mykind
     er(1,5)   =  h  - vv * c
     er(2,5)   =  qq
     er(3,5)   =  h  + vv * c
     er(4,5)   =  ww
     er(5,5)   = -uu

     do  m=1,5 ! loop on characteristic fields
      evmax(m) = -1._mykind
     enddo
     do l=1,ng2 ! LLF
      ll = j + l - iweno
      do  m=1,5
       evmax(m) = max(ev(m,ll),evmax(m))
      enddo
     enddo
     do l=1,ng2 ! loop over the stencil centered at face j
      ll = j + l - iweno
      do m=1,5
       wc = 0._mykind
       gc = 0._mykind
       do mm=1,5
        wc = wc + el(mm,m) * uj(mm,ll)
        gc = gc + el(mm,m) * fj(mm,ll)
       enddo
       gplus (m,l)   = 0.5_mykind * (gc + evmax(m) * wc)
       gminus(m,l)   = gc - gplus(m,l)
      enddo
     enddo
!
!    Reconstruction of the '+' and '-' fluxes
!
     call wenorec (5,gplus,gminus,gl,gr,iweno)
!
     do m=1,5
      ghat(m) = gl(m) + gr(m) ! char. flux
     enddo
!
!    Return to conservative fluxes
     do m=1,5
      fhat_gpu(i,j,k,m) = 0._mykind
      do mm=1,5
       fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * ghat(mm)
      enddo
     enddo
!
    endif
   enddo ! end of loop on the cell faces
!
!  Update net flux
   do j=jstart,jend ! loop on the inner nodes
    do m=1,5
     fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + (fhat_gpu(i,j,k,m)-fhat_gpu(i,j-1,k,m))*detady_gpu(j)
    enddo
   enddo
  enddo ! end of i-loop
 enddo ! end of k-loop
#endif
!
end subroutine euler_j

subroutine euler_k(kstart,kend)
!
! Evaluation of convective terms in the z direction
!
 use mod_streams
 implicit none
!
 integer :: kstart,kend
 integer :: lmax,ng2
#ifdef USE_CUDA
 type(dim3) :: grid, tBlock
 real(mykind) :: tt,st,et
 integer :: n_th_x, n_th_y
#else
 real(mykind),dimension(5,1-ng:nz+ng) :: uk,fk,ev
 real(mykind),dimension(5) :: evmax
 real(mykind),dimension(5,5) :: el,er
 real(mykind),dimension(5) :: gr,gl,ghat
 real(mykind),dimension(5,2*ng) :: gplus,gminus
 integer,dimension(1-ng:nz+ng) :: ishk
!
 real(mykind),dimension(1-ng:nz+ng) :: den
 real(mykind),dimension(6,1-ng:nz+ng) :: u,v
 real(mykind),dimension(6,ng,1-ng:nz+ng) :: uv
 real(mykind),dimension(5,nz) :: df
 real(mykind),dimension(5) :: fh
 real(mykind),dimension(6) :: ft,uvs
!
 integer :: i,j,k,l,m,n
 integer :: km,kmax,kmin,kp,ll,mm
 real(mykind) :: b1,b2,c,cc,ci,ducm,ent,gc,gpr,h,hp,pp,qq,qqp
 real(mykind) :: r,rho,rhoe,rhom,rhou,rhov,rhow,ri,rip,rp1,up
 real(mykind) :: uu,vp,vv,wc,wp,ww
#endif
!----------------------------------------------------------------------------------------
!
 ng2  = 2*iweno
 lmax = iorder/2 ! max stencil width
!
#ifdef USE_CUDA
    if (tresduc<1._mykind) then
     n_th_x = EULERWENO_THREADS_X
     n_th_y = EULERWENO_THREADS_Y
     tBlock = dim3(n_th_x,n_th_y,1)
     grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(ny)/tBlock%y),1)
     call euler_k_kernel<<<grid, tBlock>>>(nv, nx, ny, nz, ng, ng2, gamma, gm1, &
         iweno, kstart, kend, lmax, iorder, iflow, &
         w_gpu, temperature_gpu, ducros_gpu, &
         fhat_gpu, fl_gpu, dcoe_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, gplus_z, gminus_z, wv_gpu )
    else
     n_th_x = EULERCENTRAL_THREADS_X
     n_th_y = EULERCENTRAL_THREADS_Y
     tBlock = dim3(n_th_x,n_th_y,1)
     grid = dim3(ceiling(real(nx)/tBlock%x),ceiling(real(ny)/tBlock%y),1)
     call euler_k_kernel_central<<<grid, tBlock, 0>>>(nv, nx, ny, nz, ng, &
         kstart, kend, lmax, &
         temperature_gpu, fhat_gpu, fl_gpu, dcoe_gpu, dzitdz_gpu, wv_gpu)
    endif
!
#else
!
 kmin = kstart - lmax
 kmax = kend   + lmax
!
 do j=1,ny ! loop in the y-direction
  do i=1,nx ! loop in the x-direction
!  Sweep in the z-direction
   do k=kmin,kmax ! loop on all the nodes
    do m=1,5
     uk(m,k) = w_gpu(i,j,k,m)
    enddo
    rho  =  w_gpu(i,j,k,1)
    rhou =  w_gpu(i,j,k,2)
    rhov =  w_gpu(i,j,k,3)
    rhow =  w_gpu(i,j,k,4)
    rhoe =  w_gpu(i,j,k,5)
    ri   =  1._mykind/rho
    uu   =  rhou*ri
    vv   =  rhov*ri
    ww   =  rhow*ri
    qq   =  0.5_mykind*(uu*uu+vv*vv+ww*ww)
    pp   =  rho*temperature_gpu(i,j,k)
    ent  =  (rhoe+pp)*ri
    do m=1,5
     u(m,k) = ww
    enddo
    u(6,k) = 1._mykind
    v(1,k) = 1._mykind
    v(2,k) = uu
    v(3,k) = vv
    v(4,k) = ww
    v(5,k) = ent
    v(6,k) = pp
    den(k) = rho
    gpr = gamma * pp * ri
    c  =  sqrt (gpr)
    fk(1,k) =       uk(4,k)
    fk(2,k) = ww *  uk(2,k)
    fk(3,k) = ww *  uk(3,k)
    fk(4,k) = ww *  uk(4,k) + pp
    fk(5,k) = ww * (uk(5,k) + pp)
    ev(1,k) = abs(ww-c)
    ev(2,k) = abs(ww)
    ev(3,k) = abs(ww+c)
    ev(4,k) = ev(2,k)
    ev(5,k) = ev(2,k)
   enddo ! end of loop on the nodes
!
!  Compute u*v averages
!
   do k=kmin,kend
    do l=1,lmax
     rhom = den(k)+den(k+l)
     do m=1,5
      uv(m,l,k) = (u(m,k)+u(m,k+l))*(v(m,k)+v(m,k+l))*rhom
     enddo
     uv(6,l,k) = (u(6,k)+u(6,k+l))*(v(6,k)+v(6,k+l))
    enddo
   enddo
!
!  Smoothness properties
   ishk = 0
   do k=kstart-1,kend
    if (any(ducros_gpu(i,j,k-iweno+1:k+iweno))) then
     ishk(k) = 1
    endif
   enddo
!
!  Evaluation of numerical fluxes
!
   do k=kstart-1,kend ! k is the index of intermediate nodes
!
    if (ishk(k)==0) then
     ft = 0._mykind
     do l=1,lmax
      uvs = 0._mykind
      do n=0,l-1
       uvs = uvs + uv(:,l,k-n)
      enddo
      ft = ft + dcoe_gpu(l,lmax)*uvs
     enddo
     do m=1,5
      fh(m) = 0.25_mykind*ft(m)
     enddo
     fh(4) = fh(4) + 0.5_mykind*ft(6)
!
     fhat_gpu(i,j,k,1:5) = fh(:)
!
    else

     kp = k + 1
!
!    Compute Roe average
!
!    Left state (node k)
     ri        =  1._mykind/uk(1,k)
     uu        =  uk(2,k) * ri
     vv        =  uk(3,k) * ri
     ww        =  uk(4,k) * ri
     qq        =  0.5_mykind * (uu*uu  +vv*vv + ww*ww)
!    pp        =  gm1 * (uk(5,k) - uk(1,k) * qq)
!    h         =  (uk(5,k)  + pp) * ri
     h         =  gamma*uk(5,k)*ri-gm1*qq
!    Right state (node k+1)
     rip       =  1._mykind/uk(1,kp)
     up        =  uk(2,kp) * rip
     vp        =  uk(3,kp) * rip
     wp        =  uk(4,kp) * rip
     qqp       =  0.5_mykind * (up*up  +vp*vp + wp*wp)
!    ppp       =  gm1 * (uk(5,kp) - uk(1,kp) * qqp)
!    hp        =  (uk(5,kp)  + ppp) * rip
     hp        =  gamma*uk(5,kp)*rip-gm1*qqp
!    average state
     r         =  uk(1,kp)*ri
     r         =  sqrt(r)
     rp1       =  1._mykind/(r  +1._mykind)
     uu        =  (r*up  +uu)*rp1
     vv        =  (r*vp  +vv)*rp1
     ww        =  (r*wp  +ww)*rp1
     h         =  (r*hp  +h)*rp1
     qq        =  0.5_mykind * (uu*uu  + vv*vv + ww*ww)
     cc        =  gm1 * (h - qq)
     c         =  sqrt(cc)
     ci        =  1._mykind/c
!
!    left eigenvectors matrix (at Roe state)
!
     b2        =   gm1/cc
     b1        =   b2 * qq
     el(1,1)   =   0.5_mykind * (b1     + ww * ci)
     el(2,1)   =  -0.5_mykind * (b2 * uu         )
     el(3,1)   =  -0.5_mykind * (b2 * vv         )
     el(4,1)   =  -0.5_mykind * (b2 * ww +     ci)
     el(5,1)   =   0.5_mykind * b2
     el(1,2)   =   1._mykind - b1
     el(2,2)   =   b2 * uu
     el(3,2)   =   b2 * vv
     el(4,2)   =   b2 * ww
     el(5,2)   =  -b2
     el(1,3)   =   0.5_mykind * (b1     - ww * ci)
     el(2,3)   =  -0.5_mykind * (b2 * uu         )
     el(3,3)   =  -0.5_mykind * (b2 * vv         )
     el(4,3)   =  -0.5_mykind * (b2 * ww -     ci)
     el(5,3)   =   0.5_mykind * b2
     el(1,4)   =  -uu
     el(2,4)   =   1._mykind
     el(3,4)   =   0._mykind
     el(4,4)   =   0._mykind
     el(5,4)   =   0._mykind
     el(1,5)   =  -vv
     el(2,5)   =   0._mykind
     el(3,5)   =   1._mykind
     el(4,5)   =   0._mykind
     el(5,5)   =   0._mykind
!
!   right eigenvectors matrix (at Roe state)
!
     er(1,1)   =  1._mykind
     er(2,1)   =  1._mykind
     er(3,1)   =  1._mykind
     er(4,1)   =  0._mykind
     er(5,1)   =  0._mykind
     er(1,2)   =  uu
     er(2,2)   =  uu
     er(3,2)   =  uu
     er(4,2)   =  1._mykind
     er(5,2)   =  0._mykind
     er(1,3)   =  vv
     er(2,3)   =  vv
     er(3,3)   =  vv
     er(4,3)   =  0._mykind
     er(5,3)   =  1._mykind
     er(1,4)   =  ww - c
     er(2,4)   =  ww
     er(3,4)   =  ww + c
     er(4,4)   =  0._mykind
     er(5,4)   =  0._mykind
     er(1,5)   =  h  - ww * c
     er(2,5)   =  qq
     er(3,5)   =  h  + ww * c
     er(4,5)   =  uu
     er(5,5)   =  vv

     do  m=1,5 ! loop on characteristic fields
      evmax(m) = -1._mykind
     enddo
     do l=1,ng2 ! LLF
      ll = k + l - iweno
      do  m=1,5
       evmax(m) = max(ev(m,ll),evmax(m))
      enddo
     enddo
     do l=1,ng2 ! loop over the stencil centered at face k
      ll = k + l - iweno
      do m=1,5
       wc = 0._mykind
       gc = 0._mykind
       do mm=1,5
        wc = wc + el(mm,m) * uk(mm,ll)
        gc = gc + el(mm,m) * fk(mm,ll)
       enddo
       gplus (m,l)   = 0.5_mykind * (gc + evmax(m) * wc)
       gminus(m,l)   = gc - gplus(m,l)
      enddo
     enddo
!
!   Reconstruction of the '+' and '-' fluxes
!
     call wenorec (5,gplus,gminus,gl,gr,iweno)
!
     do m=1,5
      ghat(m) = gl(m) + gr(m) ! char. flux
     enddo
!
!   Return to conservative fluxes
     do m=1,5
      fhat_gpu(i,j,k,m) = 0._mykind
      do mm=1,5
       fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * ghat(mm)
      enddo
     enddo
!
    endif
   enddo ! end of loop on the cell faces
!
!  Update net flux
   do k=kstart,kend ! loop on the inner nodes
    do m=1,5
     fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + (fhat_gpu(i,j,k,m)-fhat_gpu(i,j,k-1,m)) * dzitdz_gpu(k)
    enddo
   enddo
  enddo ! end of i-loop
 enddo ! end of j-loop
#endif

 if (ibc(1)==9) then
#ifdef USE_CUDA
#else
  call move_alloc(rf_gpu, rf)
#endif
  call generateinflowrand
 endif
 !@cuf iercuda=cudaDeviceSynchronize()
 if (ibc(1)==9) then
#ifdef USE_CUDA
  rf_gpu = rf
#else
  call move_alloc(rf, rf_gpu)
#endif
 endif
!
end subroutine euler_k

#ifdef USE_CUDA
    attributes(global) subroutine euler_i_kernel_central(nv, nx, ny, nz, ng, &
        istart, iend, num_call, last_call, i1, i2, lmax, &
        fhat_trans_gpu, temperature_trans_gpu, fl_trans_gpu, dcoe_gpu, dcsidx_gpu, wv_trans_gpu)
        implicit none
        integer, parameter :: mykind = MYKIND
        ! Passed arguments
        integer, value :: nv, nx, ny, nz, ng
        integer, value :: istart, iend, lmax, num_call, last_call, i1, i2
        real(mykind) :: dcoe_gpu(4,4)
        real(mykind) :: dcsidx_gpu(nx)

        real(mykind) :: wv_trans_gpu(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,nv)
        real(mykind) :: fl_trans_gpu(1:ny,1:nx,1:nz,nv)
        real(mykind) :: fhat_trans_gpu(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,6)
        real(mykind) :: temperature_trans_gpu(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng)

        ! Local variables
        integer :: i, j, k, m, l
        real(mykind) :: rhom, uui, vvi, wwi, ppi, enti
        real(mykind) :: uuip, vvip, wwip, ppip, entip
        real(mykind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(mykind) :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uv_part
 
        j = blockDim%x * (blockIdx%x - 1) + threadIdx%x 
        k = blockDim%y * (blockIdx%y - 1) + threadIdx%y
        if (j > ny .or. k > nz) return

        do i=istart-1,iend

            ft1 = 0._mykind
            ft2 = 0._mykind
            ft3 = 0._mykind
            ft4 = 0._mykind
            ft5 = 0._mykind
            ft6 = 0._mykind
            do l=1,lmax
                uvs1 = 0._mykind
                uvs2 = 0._mykind
                uvs3 = 0._mykind
                uvs4 = 0._mykind
                uvs5 = 0._mykind
                uvs6 = 0._mykind
                do m=0,l-1
                    rhom    =  wv_trans_gpu(j,i-m,k,1) + wv_trans_gpu(j,i-m+l,k,1)

                    uui   = wv_trans_gpu(j,i-m,k,2)
                    vvi   = wv_trans_gpu(j,i-m,k,3)
                    wwi   = wv_trans_gpu(j,i-m,k,4)
                    ppi   = wv_trans_gpu(j,i-m,k,1)*temperature_trans_gpu(j,i-m,k)
                    enti  = wv_trans_gpu(j,i-m,k,5)

                    uuip   = wv_trans_gpu(j,i-m+l,k,2)
                    vvip   = wv_trans_gpu(j,i-m+l,k,3)
                    wwip   = wv_trans_gpu(j,i-m+l,k,4)
                    ppip   = wv_trans_gpu(j,i-m+l,k,1)*temperature_trans_gpu(j,i-m+l,k)
                    entip  = wv_trans_gpu(j,i-m+l,k,5)

                    uv_part = (uui+uuip) * rhom
                    uvs1 = uvs1 + uv_part * (2._mykind)
                    uvs2 = uvs2 + uv_part * (uui+uuip)
                    uvs3 = uvs3 + uv_part * (vvi+vvip)
                    uvs4 = uvs4 + uv_part * (wwi+wwip)
                    uvs5 = uvs5 + uv_part * (enti+entip)
                    uvs6 = uvs6 + (2._mykind)*(ppi+ppip)
                enddo
                ft1 = ft1 + dcoe_gpu(l,lmax)*uvs1
                ft2 = ft2 + dcoe_gpu(l,lmax)*uvs2
                ft3 = ft3 + dcoe_gpu(l,lmax)*uvs3
                ft4 = ft4 + dcoe_gpu(l,lmax)*uvs4
                ft5 = ft5 + dcoe_gpu(l,lmax)*uvs5
                ft6 = ft6 + dcoe_gpu(l,lmax)*uvs6
            enddo
!
            fhat_trans_gpu(j,i,k,1) = 0.25_mykind*ft1
            fhat_trans_gpu(j,i,k,2) = 0.25_mykind*ft2 + 0.5_mykind*ft6
            fhat_trans_gpu(j,i,k,3) = 0.25_mykind*ft3
            fhat_trans_gpu(j,i,k,4) = 0.25_mykind*ft4
            fhat_trans_gpu(j,i,k,5) = 0.25_mykind*ft5
        enddo
!
!       Update net flux 
        if (num_call == last_call) then
            do i=i1,i2 ! loop on the inner nodes
                do m=1,5
                    fl_trans_gpu(j,i,k,m) = (fhat_trans_gpu(j,i,k,m)-fhat_trans_gpu(j,i-1,k,m))*dcsidx_gpu(i)
                enddo
            enddo
        endif
!
    endsubroutine euler_i_kernel_central

    attributes(global) subroutine euler_i_kernel(nv, nx, ny, nz, ng, ng2, gamma, gm1, &
        iweno, istart, iend, num_call, last_call, i1, i2, lmax, iorder, iflow, &
        w_gpu, temperature_gpu, ducros_gpu, &
        fhat_trans_gpu, temperature_trans_gpu, fl_trans_gpu, & 
        fhat_gpu, fl_gpu, dcoe_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, gplus, gminus, wv_trans_gpu )
        implicit none
        integer, parameter :: mykind = MYKIND
        ! Passed arguments
        integer, value :: nv, nx, ny, nz, ng, ng2, iorder, iflow
        real(mykind), value :: gamma, gm1
        integer, value :: iweno, istart, iend, lmax, num_call, last_call, i1, i2
        real(mykind) :: wv_trans_gpu(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,nv)
        real(mykind) :: w_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,nv)
        real(mykind) :: temperature_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng)
        logical :: ducros_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng)
        real(mykind) :: fhat_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,6)
        real(mykind) :: fl_gpu(nx,ny,nz,nv)
        real(mykind) :: dcoe_gpu(4,4)
        real(mykind) :: dcsidx_gpu(nx), detady_gpu(ny), dzitdz_gpu(nz)
        real(mykind) :: gplus(5,2*iweno,ny,nz)
        real(mykind) :: gminus(5,2*iweno,ny,nz)

        real(mykind) :: fl_trans_gpu(1:ny,1:nx,1:nz,nv)
        real(mykind) :: fhat_trans_gpu(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng,6)
        real(mykind) :: temperature_trans_gpu(1-ng:ny+ng,1-ng:nx+ng,1-ng:nz+ng)

        ! Local variables
        integer :: i, j, k, im, m, ii, l, ip, ll, mm
        logical :: ishk
        real(mykind) :: fh1, fh2, fh3, fh4, fh5
        real(mykind) :: rhom, uui, vvi, wwi, ppi, enti
        real(mykind) :: uuip, vvip, wwip, ppip, entip
        real(mykind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(mykind) :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6
        real(mykind) :: ri, uu, vv, ww, qq, h
        real(mykind) :: rip, up, vp, wp, qqp, hp
        real(mykind) :: r, rp1, cc, c, ci, b1, b2
        real(mykind) :: gpr, pp
        real(mykind) :: df
        ! for weno
        real(mykind), dimension(5) :: ev, evmax, fi, ghat, gl, gr
        real(mykind), dimension(5,5) :: el, er
        real(mykind) :: rho, rhou, wc, gc
 
        j = blockDim%x * (blockIdx%x - 1) + threadIdx%x 
        k = blockDim%y * (blockIdx%y - 1) + threadIdx%y
        if (j > ny .or. k > nz) return

        do i=istart-1,iend

            ishk = .false.
            do ii=i-iweno+1,i+iweno
                if (ducros_gpu(ii,j,k)) ishk = .true.
            enddo
!
            if (.not.ishk) then
                ft1 = 0._mykind
                ft2 = 0._mykind
                ft3 = 0._mykind
                ft4 = 0._mykind
                ft5 = 0._mykind
                ft6 = 0._mykind
                do l=1,lmax
                    uvs1 = 0._mykind
                    uvs2 = 0._mykind
                    uvs3 = 0._mykind
                    uvs4 = 0._mykind
                    uvs5 = 0._mykind
                    uvs6 = 0._mykind
                    do m=0,l-1
                        rhom  = wv_trans_gpu(j,i-m,k,1) + wv_trans_gpu(j,i-m+l,k,1) 

                        uui   = wv_trans_gpu(j,i-m,k,2)
                        vvi   = wv_trans_gpu(j,i-m,k,3)
                        wwi   = wv_trans_gpu(j,i-m,k,4)
                        ppi   = wv_trans_gpu(j,i-m,k,1)*temperature_trans_gpu(j,i-m,k)
                        enti  = wv_trans_gpu(j,i-m,k,5)

                        uuip   = wv_trans_gpu(j,i-m+l,k,2)
                        vvip   = wv_trans_gpu(j,i-m+l,k,3)
                        wwip   = wv_trans_gpu(j,i-m+l,k,4)
                        ppip   = wv_trans_gpu(j,i-m+l,k,1)*temperature_trans_gpu(j,i-m+l,k)
                        entip  = wv_trans_gpu(j,i-m+l,k,5)

                        uvs1 = uvs1 + (uui+uuip)*(2._mykind)*rhom
                        uvs2 = uvs2 + (uui+uuip)*(uui+uuip)*rhom
                        uvs3 = uvs3 + (uui+uuip)*(vvi+vvip)*rhom
                        uvs4 = uvs4 + (uui+uuip)*(wwi+wwip)*rhom
                        uvs5 = uvs5 + (uui+uuip)*(enti+entip)*rhom
                        uvs6 = uvs6 + (2._mykind)*(ppi+ppip)
                    enddo
                    ft1 = ft1 + dcoe_gpu(l,lmax)*uvs1
                    ft2 = ft2 + dcoe_gpu(l,lmax)*uvs2
                    ft3 = ft3 + dcoe_gpu(l,lmax)*uvs3
                    ft4 = ft4 + dcoe_gpu(l,lmax)*uvs4
                    ft5 = ft5 + dcoe_gpu(l,lmax)*uvs5
                    ft6 = ft6 + dcoe_gpu(l,lmax)*uvs6
                enddo
!
                fhat_trans_gpu(j,i,k,1) = 0.25_mykind*ft1
                fhat_trans_gpu(j,i,k,2) = 0.25_mykind*ft2 + 0.5_mykind*ft6
                fhat_trans_gpu(j,i,k,3) = 0.25_mykind*ft3
                fhat_trans_gpu(j,i,k,4) = 0.25_mykind*ft4
                fhat_trans_gpu(j,i,k,5) = 0.25_mykind*ft5
            else
                ip = i + 1
!
!               Compute Roe average
!
!               Left state (node i)
                ri        =  1._mykind/wv_trans_gpu(j,i,k,1)
                uu        =  wv_trans_gpu(j,i,k,2)
                vv        =  wv_trans_gpu(j,i,k,3)
                ww        =  wv_trans_gpu(j,i,k,4)
                qq        =  0.5_mykind * (uu*uu  +vv*vv + ww*ww)
                h         =  wv_trans_gpu(j,i,k,5)
!               Right state (node i+1)
                rip       =  1._mykind/wv_trans_gpu(j,ip,k,1)
                up        =  wv_trans_gpu(j,ip,k,2)
                vp        =  wv_trans_gpu(j,ip,k,3)
                wp        =  wv_trans_gpu(j,ip,k,4)
                qqp       =  0.5_mykind * (up*up  +vp*vp +wp*wp)
                hp        =  wv_trans_gpu(j,ip,k,5)
!               average state
                r         =  wv_trans_gpu(j,ip,k,1)*ri
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

!               left eigenvectors matrix (at Roe state)
!
                b2        =   gm1/cc
                b1        =   b2 * qq
                el(1,1)   =   0.5_mykind * (b1     + uu * ci)
                el(2,1)   =  -0.5_mykind * (b2 * uu +     ci)
                el(3,1)   =  -0.5_mykind * (b2 * vv         )
                el(4,1)   =  -0.5_mykind * (b2 * ww         )
                el(5,1)   =   0.5_mykind * b2
                el(1,2)   =   1._mykind - b1
                el(2,2)   =   b2*uu
                el(3,2)   =   b2*vv
                el(4,2)   =   b2*ww
                el(5,2)   =  -b2
                el(1,3)   =   0.5_mykind * (b1     - uu * ci)
                el(2,3)   =  -0.5_mykind * (b2 * uu -     ci)
                el(3,3)   =  -0.5_mykind * (b2 * vv         )
                el(4,3)   =  -0.5_mykind * (b2 * ww         )
                el(5,3)   =   0.5_mykind * b2
                el(1,4)   =  -vv
                el(2,4)   =   0._mykind
                el(3,4)   =   1._mykind
                el(4,4)   =   0._mykind
                el(5,4)   =   0._mykind
                el(1,5)   =  -ww
                el(2,5)   =   0._mykind
                el(3,5)   =   0._mykind
                el(4,5)   =   1._mykind
                el(5,5)   =   0._mykind
!
!               right eigenvectors matrix (at Roe state)
!
                er(1,1)   =  1._mykind
                er(2,1)   =  1._mykind
                er(3,1)   =  1._mykind
                er(4,1)   =  0._mykind
                er(5,1)   =  0._mykind
                er(1,2)   =  uu -  c
                er(2,2)   =  uu
                er(3,2)   =  uu +  c
                er(4,2)   =  0._mykind
                er(5,2)   =  0._mykind
                er(1,3)   =  vv
                er(2,3)   =  vv
                er(3,3)   =  vv
                er(4,3)   =  1._mykind
                er(5,3)   =  0._mykind
                er(1,4)   =  ww
                er(2,4)   =  ww
                er(3,4)   =  ww
                er(4,4)   =  0._mykind
                er(5,4)   =  1._mykind
                er(1,5)   =  h  - uu * c
                er(2,5)   =  qq
                er(3,5)   =  h  + uu * c
                er(4,5)   =  vv
                er(5,5)   =  ww

                do m=1,5 ! loop on characteristic fields
                    evmax(m) = -1._mykind
                enddo
                do l=1,ng2 ! LLF
                    ll = i + l - iweno
!
                    rho  = wv_trans_gpu(j,ll,k,1)
                    rhou = wv_trans_gpu(j,ll,k,2)*rho
                    ri   = 1._mykind/rho
                    uu   = rhou*ri
                    pp   = rho*temperature_trans_gpu(j,ll,k)
                    gpr  = gamma * pp * ri
                    c    = sqrt (gpr)
                    ev(1) = abs(uu-c)
                    ev(2) = abs(uu)
                    ev(3) = abs(uu+c)
                    ev(4) = ev(2)
                    ev(5) = ev(2)
                    do m=1,5
                        evmax(m) = max(ev(m),evmax(m))
                    enddo
                enddo
                do l=1,ng2 ! loop over the stencil centered at face i
                    ll = i + l - iweno

                    rho  = wv_trans_gpu(j,ll,k,1)
                    rhou = wv_trans_gpu(j,ll,k,2)*rho
                    ri   = 1._mykind/rho
                    uu   = rhou*ri
                    pp   = rho*temperature_trans_gpu(j,ll,k)
                    fi(1)  =       rho * wv_trans_gpu(j,ll,k,2)
                    fi(2)  = uu *  rho * wv_trans_gpu(j,ll,k,2)  + pp
                    fi(3)  = uu *  rho * wv_trans_gpu(j,ll,k,3)
                    fi(4)  = uu *  rho * wv_trans_gpu(j,ll,k,4)
                    fi(5)  = uu *  rho * wv_trans_gpu(j,ll,k,5) 
                    do m=1,5
                        wc = 0._mykind
                        gc = 0._mykind

                        do mm=1,1
                            wc = wc + el(mm,m) * rho
                            gc = gc + el(mm,m) * fi(mm)
                        enddo
                        do mm=2,4
                            wc = wc + el(mm,m) * wv_trans_gpu(j,ll,k,mm)*rho
                            gc = gc + el(mm,m) * fi(mm)
                        enddo
                        do mm=5,5
                            wc = wc + el(mm,m) * (wv_trans_gpu(j,ll,k,mm) * rho - pp)
                            gc = gc + el(mm,m) * fi(mm)
                        enddo
                        gplus (m,l,j,k) = 0.5_mykind * (gc + evmax(m) * wc)
                        gminus(m,l,j,k) = gc - gplus(m,l,j,k)
                    enddo
                enddo
!
!               Reconstruction of the '+' and '-' fluxes
!
                call wenorec(5,gplus(:,:,j,k),gminus(:,:,j,k),gl,gr,iweno)
!
                do m=1,5
                    ghat(m) = gl(m) + gr(m) ! char. flux
                enddo

!               !Return to conservative fluxes
                do m=1,5
                    fhat_trans_gpu(j,i,k,m) = 0._mykind
                    do mm=1,5
                       fhat_trans_gpu(j,i,k,m) = fhat_trans_gpu(j,i,k,m) + er(mm,m) * ghat(mm)
                    enddo
                enddo

            endif
        enddo
!
!       Update net flux 
        if (num_call == last_call) then
            do i=i1,i2 ! loop on the inner nodes
                do m=1,5
                    fl_trans_gpu(j,i,k,m) = (fhat_trans_gpu(j,i,k,m)-fhat_trans_gpu(j,i-1,k,m))*dcsidx_gpu(i)
                enddo
            enddo
        endif
!
    endsubroutine euler_i_kernel

    attributes(global) subroutine euler_j_kernel_central(nv, nx, ny, nz, ng, &
        jstart, jend, lmax, iflow, &
        temperature_gpu, fhat_gpu, fl_gpu, dcoe_gpu, detady_gpu, wv_gpu)
        implicit none
        integer, parameter :: mykind = MYKIND
        ! Passed arguments
        integer, value :: nv, nx, ny, nz, ng, iflow
        integer, value :: jstart, jend, lmax
        real(mykind) :: wv_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,nv)
        real(mykind) :: temperature_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng)
        real(mykind) :: fhat_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,6)
        real(mykind) :: fl_gpu(nx,ny,nz,nv)
        real(mykind) :: dcoe_gpu(4,4)
        real(mykind) :: detady_gpu(ny)
        ! Local variables
        integer :: i, j, k, m, l
        real(mykind) :: fh1, fh2, fh3, fh4, fh5
        real(mykind) :: rhom, uui, vvi, wwi, ppi, enti
        real(mykind) :: uuip, vvip, wwip, ppip, entip
        real(mykind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(mykind) :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uv_part
 
        i = blockDim%x * (blockIdx%x - 1) + threadIdx%x 
        k = blockDim%y * (blockIdx%y - 1) + threadIdx%y
        if (i > nx .or. k > nz) return

        do j=jstart-1,jend
            ft1 = 0._mykind
            ft2 = 0._mykind
            ft3 = 0._mykind
            ft4 = 0._mykind
            ft5 = 0._mykind
            ft6 = 0._mykind
            do l=1,lmax
                uvs1 = 0._mykind
                uvs2 = 0._mykind
                uvs3 = 0._mykind
                uvs4 = 0._mykind
                uvs5 = 0._mykind
                uvs6 = 0._mykind
                do m=0,l-1
                    rhom  = wv_gpu(i,j-m,k,1) + wv_gpu(i,j-m+l,k,1) 

                    uui   = wv_gpu(i,j-m,k,2)
                    vvi   = wv_gpu(i,j-m,k,3)
                    wwi   = wv_gpu(i,j-m,k,4)
                    ppi   = wv_gpu(i,j-m,k,1)*temperature_gpu(i,j-m,k)
                    enti  = wv_gpu(i,j-m,k,5)

                    uuip  = wv_gpu(i,j-m+l,k,2)
                    vvip  = wv_gpu(i,j-m+l,k,3)
                    wwip  = wv_gpu(i,j-m+l,k,4)
                    ppip  = wv_gpu(i,j-m+l,k,1)*temperature_gpu(i,j-m+l,k)
                    entip = wv_gpu(i,j-m+l,k,5)

                    uv_part = (vvi+vvip) * rhom
                    uvs1 = uvs1 + uv_part * (2._mykind)
                    uvs2 = uvs2 + uv_part * (uui+uuip)
                    uvs3 = uvs3 + uv_part * (vvi+vvip)
                    uvs4 = uvs4 + uv_part * (wwi+wwip)
                    uvs5 = uvs5 + uv_part * (enti+entip)
                    uvs6 = uvs6 + (2._mykind)*(ppi+ppip)
                enddo
                ft1 = ft1 + dcoe_gpu(l,lmax)*uvs1
                ft2 = ft2 + dcoe_gpu(l,lmax)*uvs2
                ft3 = ft3 + dcoe_gpu(l,lmax)*uvs3
                ft4 = ft4 + dcoe_gpu(l,lmax)*uvs4
                ft5 = ft5 + dcoe_gpu(l,lmax)*uvs5
                ft6 = ft6 + dcoe_gpu(l,lmax)*uvs6
            enddo
            fh1 = 0.25_mykind*ft1
            fh2 = 0.25_mykind*ft2
            fh3 = 0.25_mykind*ft3
            fh4 = 0.25_mykind*ft4
            fh5 = 0.25_mykind*ft5
            if (iflow==0) then
                if (j==0.or.j==ny) then
                    fh1 = 0._mykind
                    fh2 = 0._mykind
                    fh3 = 0._mykind
                    fh4 = 0._mykind
                    fh5 = 0._mykind
                endif
            endif
            fh3 = fh3 + 0.5_mykind*ft6

            fhat_gpu(i,j,k,1) = fh1
            fhat_gpu(i,j,k,2) = fh2
            fhat_gpu(i,j,k,3) = fh3
            fhat_gpu(i,j,k,4) = fh4
            fhat_gpu(i,j,k,5) = fh5
        enddo

!       Update net flux 
        do j=jstart,jend ! loop on the inner nodes
            do m=1,5
                fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + (fhat_gpu(i,j,k,m)-fhat_gpu(i,j-1,k,m))*detady_gpu(j)
            enddo
        enddo

    endsubroutine euler_j_kernel_central

    attributes(global) subroutine euler_j_kernel(nv, nx, ny, nz, ng, ng2, gamma, gm1, &
        iweno, jstart, jend, lmax, iorder, iflow, &
        w_gpu, temperature_gpu, ducros_gpu, &
        fhat_gpu, fl_gpu, dcoe_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, gplus, gminus, wv_gpu )
        implicit none
        integer, parameter :: mykind = MYKIND
        ! Passed arguments
        integer, value :: nv, nx, ny, nz, ng, ng2, iorder, iflow
        real(mykind), value :: gamma, gm1
        integer, value :: iweno, jstart, jend, lmax
        real(mykind) :: w_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,nv)
        real(mykind) :: wv_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,nv)
        real(mykind) :: temperature_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng)
        logical :: ducros_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng)
        real(mykind) :: fhat_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,6)
        real(mykind) :: fl_gpu(nx,ny,nz,nv)
        real(mykind) :: dcoe_gpu(4,4)
        real(mykind) :: dcsidx_gpu(nx), detady_gpu(ny), dzitdz_gpu(nz)
        real(mykind) :: gplus(5,2*iweno,nx,nz)
        real(mykind) :: gminus(5,2*iweno,nx,nz)
        ! Local variables
        integer :: i, j, k, jm, m, jj, l, ip, jp, ll, mm
        logical :: ishk
        real(mykind) :: fh1, fh2, fh3, fh4, fh5
        real(mykind) :: rhom, uui, vvi, wwi, ppi, enti
        real(mykind) :: uuip, vvip, wwip, ppip, entip
        real(mykind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(mykind) :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6
        real(mykind) :: ri, uu, vv, ww, qq, h
        real(mykind) :: rip, up, vp, wp, qqp, hp
        real(mykind) :: r, rp1, cc, c, ci, b1, b2
        real(mykind) :: gpr, pp
        real(mykind) :: df
        ! for weno
        real(mykind), dimension(5) :: ev, evmax, fi, ghat, gl, gr
        real(mykind), dimension(5,5) :: el, er
        real(mykind) :: rho, rhov, wc, gc
 
        i = blockDim%x * (blockIdx%x - 1) + threadIdx%x 
        k = blockDim%y * (blockIdx%y - 1) + threadIdx%y
        if(i > nx .or. k > nz) return

        do j=jstart-1,jend

            ishk = .false.
            do jj=j-iweno+1,j+iweno
                if (ducros_gpu(i,jj,k)) ishk = .true.
            enddo
!
            if (.not.ishk) then
                ft1 = 0._mykind
                ft2 = 0._mykind
                ft3 = 0._mykind
                ft4 = 0._mykind
                ft5 = 0._mykind
                ft6 = 0._mykind
                do l=1,lmax
                    uvs1 = 0._mykind
                    uvs2 = 0._mykind
                    uvs3 = 0._mykind
                    uvs4 = 0._mykind
                    uvs5 = 0._mykind
                    uvs6 = 0._mykind
                    do m=0,l-1
                        rhom  = wv_gpu(i,j-m,k,1) + wv_gpu(i,j-m+l,k,1) 

                        uui   = wv_gpu(i,j-m,k,2)
                        vvi   = wv_gpu(i,j-m,k,3)
                        wwi   = wv_gpu(i,j-m,k,4)
                        ppi   = wv_gpu(i,j-m,k,1)*temperature_gpu(i,j-m,k)
                        enti  = wv_gpu(i,j-m,k,5)

                        uuip  = wv_gpu(i,j-m+l,k,2)
                        vvip  = wv_gpu(i,j-m+l,k,3)
                        wwip  = wv_gpu(i,j-m+l,k,4)
                        ppip  = wv_gpu(i,j-m+l,k,1)*temperature_gpu(i,j-m+l,k)
                        entip = wv_gpu(i,j-m+l,k,5)

                        uvs1 = uvs1 + (vvi+vvip)*(2._mykind)*rhom
                        uvs2 = uvs2 + (vvi+vvip)*(uui+uuip)*rhom
                        uvs3 = uvs3 + (vvi+vvip)*(vvi+vvip)*rhom
                        uvs4 = uvs4 + (vvi+vvip)*(wwi+wwip)*rhom
                        uvs5 = uvs5 + (vvi+vvip)*(enti+entip)*rhom
                        uvs6 = uvs6 + (2._mykind)*(ppi+ppip)
                    enddo
                    ft1 = ft1 + dcoe_gpu(l,lmax)*uvs1
                    ft2 = ft2 + dcoe_gpu(l,lmax)*uvs2
                    ft3 = ft3 + dcoe_gpu(l,lmax)*uvs3
                    ft4 = ft4 + dcoe_gpu(l,lmax)*uvs4
                    ft5 = ft5 + dcoe_gpu(l,lmax)*uvs5
                    ft6 = ft6 + dcoe_gpu(l,lmax)*uvs6
                enddo
                fh1 = 0.25_mykind*ft1
                fh2 = 0.25_mykind*ft2
                fh3 = 0.25_mykind*ft3
                fh4 = 0.25_mykind*ft4
                fh5 = 0.25_mykind*ft5
                if (iflow==0) then
                    if (j==0.or.j==ny) then
                        fh1 = 0._mykind
                        fh2 = 0._mykind
                        fh3 = 0._mykind
                        fh4 = 0._mykind
                        fh5 = 0._mykind
                    endif
                endif
                fh3 = fh3 + 0.5_mykind*ft6

                fhat_gpu(i,j,k,1) = fh1
                fhat_gpu(i,j,k,2) = fh2
                fhat_gpu(i,j,k,3) = fh3
                fhat_gpu(i,j,k,4) = fh4
                fhat_gpu(i,j,k,5) = fh5
            else
!
                jp = j + 1
!
!               Compute Roe average
!
!               Left state (node i)
                ri        =  1._mykind/w_gpu(i,j,k,1)
                uu        =  w_gpu(i,j,k,2) * ri
                vv        =  w_gpu(i,j,k,3) * ri
                ww        =  w_gpu(i,j,k,4) * ri
                qq        =  0.5_mykind * (uu*uu  +vv*vv + ww*ww)
!               pp        =  gm1 * (w_gpu(i,j,k,5) - w_gpu(i,j,k,1) * qq)
!               h         =  (w_gpu(i,j,k,5)  + pp) * ri
                h         =  gamma*w_gpu(i,j,k,5)*ri-gm1*qq
!               Right state (node i+1)
                rip       =  1._mykind/w_gpu(i,jp,k,1)
                up        =  w_gpu(i,jp,k,2) * rip
                vp        =  w_gpu(i,jp,k,3) * rip
                wp        =  w_gpu(i,jp,k,4) * rip
                qqp       =  0.5_mykind * (up*up  +vp*vp +wp*wp)
!               ppp       =  gm1 * (w_gpu(i,jp,k,5) - w_gpu(i,jp,k,1) * qqp)
!               hp        =  (w_gpu(i,jp,k,5)  + ppp) * rip
                hp        =  gamma*w_gpu(i,jp,k,5)*rip-gm1*qqp
!               average state
                r         =  w_gpu(i,jp,k,1)*ri
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
!               left eigenvectors matrix (at Roe state)
!
                b2     =   gm1/cc
                b1     =   b2 * qq
                el(1,1) =   0.5_mykind * (b1     + vv * ci)
                el(2,1) =  -0.5_mykind * (b2 * uu         )
                el(3,1) =  -0.5_mykind * (b2 * vv +     ci)
                el(4,1) =  -0.5_mykind * (b2 * ww         )
                el(5,1) =   0.5_mykind * b2
                el(1,2) =   1._mykind - b1
                el(2,2) =   b2 * uu
                el(3,2) =   b2 * vv
                el(4,2) =   b2 * ww
                el(5,2) =  -b2
                el(1,3) =   0.5_mykind * (b1     - vv * ci)
                el(2,3) =  -0.5_mykind * (b2 * uu         )
                el(3,3) =  -0.5_mykind * (b2 * vv -     ci)
                el(4,3) =  -0.5_mykind * (b2 * ww         )
                el(5,3) =   0.5_mykind * b2
                el(1,4) =  -ww
                el(2,4) =   0._mykind
                el(3,4) =   0._mykind
                el(4,4) =   1._mykind
                el(5,4) =   0._mykind
                el(1,5) =   uu
                el(2,5) =  -1._mykind
                el(3,5) =   0._mykind
                el(4,5) =   0._mykind
                el(5,5) =   0._mykind
!
!               right eigenvectors matrix (at Roe state)
!
                er(1,1) =  1._mykind
                er(2,1) =  1._mykind
                er(3,1) =  1._mykind
                er(4,1) =  0._mykind
                er(5,1) =  0._mykind
                er(1,2) =  uu
                er(2,2) =  uu
                er(3,2) =  uu
                er(4,2) =  0._mykind
                er(5,2) = -1._mykind
                er(1,3) =  vv - c
                er(2,3) =  vv
                er(3,3) =  vv + c
                er(4,3) =  0._mykind
                er(5,3) =  0._mykind
                er(1,4) =  ww
                er(2,4) =  ww
                er(3,4) =  ww
                er(4,4) =  1._mykind
                er(5,4) =  0._mykind
                er(1,5) =  h  - vv * c
                er(2,5) =  qq
                er(3,5) =  h  + vv * c
                er(4,5) =  ww
                er(5,5) = -uu

                do m=1,5 ! loop on characteristic fields
                    evmax(m) = -1._mykind
                enddo
                do l=1,ng2 ! LLF
                    ll = j + l - iweno
!
                    rho  = w_gpu(i,ll,k,1)
                    rhov = w_gpu(i,ll,k,3)
                    ri   = 1._mykind/rho
                    vv   = rhov*ri
                    pp   = rho*temperature_gpu(i,ll,k)
                    gpr  = gamma * pp * ri
                    c    = sqrt (gpr)
                    ev(1) = abs(vv-c)
                    ev(2) = abs(vv)
                    ev(3) = abs(vv+c)
                    ev(4) = ev(2)
                    ev(5) = ev(2)
                    do m=1,5
                        evmax(m) = max(ev(m),evmax(m))
                    enddo
                enddo
                do l=1,ng2 ! loop over the stencil centered at face i
                    ll = j + l - iweno

                    rho  = w_gpu(i,ll,k,1)
                    rhov = w_gpu(i,ll,k,3)
                    ri   = 1._mykind/rho
                    vv   = rhov*ri
                    pp   = rho*temperature_gpu(i,ll,k)
                    fi(1)  =       w_gpu(i,ll,k,3)
                    fi(2)  = vv *  w_gpu(i,ll,k,2)
                    fi(3)  = vv *  w_gpu(i,ll,k,3)  + pp
                    fi(4)  = vv *  w_gpu(i,ll,k,4)
                    fi(5)  = vv * (w_gpu(i,ll,k,5)  + pp)
                    do m=1,5
                        wc = 0._mykind
                        gc = 0._mykind

                        do mm=1,5
                            wc = wc + el(mm,m) * w_gpu(i,ll,k,mm)
                            gc = gc + el(mm,m) * fi(mm)
                        enddo
                        gplus (m,l,i,k) = 0.5_mykind * (gc + evmax(m) * wc)
                        gminus(m,l,i,k) = gc - gplus(m,l,i,k)
                    enddo
                enddo
!
!               Reconstruction of the '+' and '-' fluxes
!
                call wenorec(5,gplus(:,:,i,k),gminus(:,:,i,k),gl,gr,iweno)
!
                do m=1,5
                    ghat(m) = gl(m) + gr(m) ! char. flux
                enddo

!               !Return to conservative fluxes
                do m=1,5
                    fhat_gpu(i,j,k,m) = 0._mykind
                    do mm=1,5
                       fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * ghat(mm)
                    enddo
                enddo

            endif
        enddo

!       Update net flux 
        do j=jstart,jend ! loop on the inner nodes
            do m=1,5
                df = (fhat_gpu(i,j,k,m)-fhat_gpu(i,j-1,k,m))*detady_gpu(j)
                fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + df
            enddo
        enddo

    endsubroutine euler_j_kernel

    attributes(global) subroutine euler_k_kernel_central(nv, nx, ny, nz, ng, &
        kstart, kend, lmax, &
        temperature_gpu, fhat_gpu, fl_gpu, dcoe_gpu, dzitdz_gpu, wv_gpu)
        implicit none
        integer, parameter :: mykind = MYKIND
        ! Passed arguments
        integer, value :: nv, nx, ny, nz, ng
        integer, value :: kstart, kend, lmax
        real(mykind) :: wv_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,nv)
        real(mykind) :: temperature_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng)
        real(mykind) :: fhat_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,6)
        real(mykind) :: fl_gpu(nx,ny,nz,nv)
        real(mykind) :: dcoe_gpu(4,4)
        real(mykind) :: dzitdz_gpu(nz)
        ! Local variables
        integer :: i, j, k, m, l
        real(mykind) :: rhom, uui, vvi, wwi, ppi, enti
        real(mykind) :: uuip, vvip, wwip, ppip, entip
        real(mykind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(mykind) :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uv_part
 
        i = blockDim%x * (blockIdx%x - 1) + threadIdx%x 
        j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
        if(i > nx .or. j > ny) return

        do k=kstart-1,kend
            ft1 = 0._mykind
            ft2 = 0._mykind
            ft3 = 0._mykind
            ft4 = 0._mykind
            ft5 = 0._mykind
            ft6 = 0._mykind
            do l=1,lmax
                uvs1 = 0._mykind
                uvs2 = 0._mykind
                uvs3 = 0._mykind
                uvs4 = 0._mykind
                uvs5 = 0._mykind
                uvs6 = 0._mykind
                do m=0,l-1
                    rhom  = wv_gpu(i,j,k-m,1) + wv_gpu(i,j,k-m+l,1) 

                    uui   = wv_gpu(i,j,k-m,2)
                    vvi   = wv_gpu(i,j,k-m,3)
                    wwi   = wv_gpu(i,j,k-m,4)
                    ppi   = wv_gpu(i,j,k-m,1)*temperature_gpu(i,j,k-m)
                    enti  = wv_gpu(i,j,k-m,5)

                    uuip   = wv_gpu(i,j,k-m+l,2)
                    vvip   = wv_gpu(i,j,k-m+l,3)
                    wwip   = wv_gpu(i,j,k-m+l,4)
                    ppip   = wv_gpu(i,j,k-m+l,1)*temperature_gpu(i,j,k-m+l)
                    entip  = wv_gpu(i,j,k-m+l,5)

                    uv_part = (wwi+wwip) * rhom
                    uvs1 = uvs1 + uv_part * (2._mykind)
                    uvs2 = uvs2 + uv_part * (uui+uuip)
                    uvs3 = uvs3 + uv_part * (vvi+vvip)
                    uvs4 = uvs4 + uv_part * (wwi+wwip)
                    uvs5 = uvs5 + uv_part * (enti+entip)
                    uvs6 = uvs6 + (2._mykind)*(ppi+ppip)
                enddo
                ft1 = ft1 + dcoe_gpu(l,lmax)*uvs1
                ft2 = ft2 + dcoe_gpu(l,lmax)*uvs2
                ft3 = ft3 + dcoe_gpu(l,lmax)*uvs3
                ft4 = ft4 + dcoe_gpu(l,lmax)*uvs4
                ft5 = ft5 + dcoe_gpu(l,lmax)*uvs5
                ft6 = ft6 + dcoe_gpu(l,lmax)*uvs6
            enddo
            fhat_gpu(i,j,k,1) = 0.25_mykind*ft1
            fhat_gpu(i,j,k,2) = 0.25_mykind*ft2
            fhat_gpu(i,j,k,3) = 0.25_mykind*ft3
            fhat_gpu(i,j,k,4) = 0.25_mykind*ft4 + 0.5_mykind*ft6
            fhat_gpu(i,j,k,5) = 0.25_mykind*ft5
        enddo

!       Update net flux 
        do k=kstart,kend ! loop on the inner nodes
            do m=1,5
                fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + (fhat_gpu(i,j,k,m)-fhat_gpu(i,j,k-1,m))*dzitdz_gpu(k)
            enddo
        enddo

    endsubroutine euler_k_kernel_central

    attributes(global) subroutine euler_k_kernel(nv, nx, ny, nz, ng, ng2, gamma, gm1, &
        iweno, kstart, kend, lmax, iorder, iflow, &
        w_gpu, temperature_gpu, ducros_gpu, &
        fhat_gpu, fl_gpu, dcoe_gpu, dcsidx_gpu, detady_gpu, dzitdz_gpu, gplus, gminus, wv_gpu )
        implicit none
        integer, parameter :: mykind = MYKIND
        ! Passed arguments
        integer, value :: nv, nx, ny, nz, ng, ng2, iorder, iflow
        real(mykind), value :: gamma, gm1
        integer, value :: iweno, kstart, kend, lmax
        real(mykind) :: w_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,nv)
        real(mykind) :: wv_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,nv)
        real(mykind) :: temperature_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng)
        logical :: ducros_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng)
        real(mykind) :: fhat_gpu(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng,6)
        real(mykind) :: fl_gpu(nx,ny,nz,nv)
        real(mykind) :: dcoe_gpu(4,4)
        real(mykind) :: dcsidx_gpu(nx), detady_gpu(ny), dzitdz_gpu(nz)
        real(mykind) :: gplus(5,2*iweno,nx,ny)
        real(mykind) :: gminus(5,2*iweno,nx,ny)
        ! Local variables
        integer :: i, j, k, km, m, kk, l, ip, jp, kp, ll, mm
        logical :: ishk
        real(mykind) :: fh1, fh2, fh3, fh4, fh5
        real(mykind) :: rhom, uui, vvi, wwi, ppi, enti
        real(mykind) :: uuip, vvip, wwip, ppip, entip
        real(mykind) :: ft1, ft2, ft3, ft4, ft5, ft6
        real(mykind) :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6
        real(mykind) :: ri, uu, vv, ww, qq, h
        real(mykind) :: rip, up, vp, wp, qqp, hp
        real(mykind) :: r, rp1, cc, c, ci, b1, b2
        real(mykind) :: gpr, pp
        real(mykind) :: df
        ! for weno
        real(mykind), dimension(5) :: ev, evmax, fi, ghat, gl, gr
        real(mykind), dimension(5,5) :: el, er
        real(mykind) :: rho, rhow, wc, gc
 
        i = blockDim%x * (blockIdx%x - 1) + threadIdx%x 
        j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
        if(i > nx .or. j > ny) return

        do k=kstart-1,kend

            ishk = .false.
            do kk=k-iweno+1,k+iweno
                if (ducros_gpu(i,j,kk)) ishk = .true.
            enddo
!
            if (.not.ishk) then
                ft1 = 0._mykind
                ft2 = 0._mykind
                ft3 = 0._mykind
                ft4 = 0._mykind
                ft5 = 0._mykind
                ft6 = 0._mykind
                do l=1,lmax
                    uvs1 = 0._mykind
                    uvs2 = 0._mykind
                    uvs3 = 0._mykind
                    uvs4 = 0._mykind
                    uvs5 = 0._mykind
                    uvs6 = 0._mykind
                    do m=0,l-1
                        rhom  = wv_gpu(i,j,k-m,1) + wv_gpu(i,j,k-m+l,1) 

                        uui   = wv_gpu(i,j,k-m,2)
                        vvi   = wv_gpu(i,j,k-m,3)
                        wwi   = wv_gpu(i,j,k-m,4)
                        ppi   = wv_gpu(i,j,k-m,1)*temperature_gpu(i,j,k-m)
                        enti  = wv_gpu(i,j,k-m,5)

                        uuip   = wv_gpu(i,j,k-m+l,2)
                        vvip   = wv_gpu(i,j,k-m+l,3)
                        wwip   = wv_gpu(i,j,k-m+l,4)
                        ppip   = wv_gpu(i,j,k-m+l,1)*temperature_gpu(i,j,k-m+l)
                        entip  = wv_gpu(i,j,k-m+l,5)

                        uvs1 = uvs1 + (wwi+wwip)*(2._mykind)*rhom
                        uvs2 = uvs2 + (wwi+wwip)*(uui+uuip)*rhom
                        uvs3 = uvs3 + (wwi+wwip)*(vvi+vvip)*rhom
                        uvs4 = uvs4 + (wwi+wwip)*(wwi+wwip)*rhom
                        uvs5 = uvs5 + (wwi+wwip)*(enti+entip)*rhom
                        uvs6 = uvs6 + (2._mykind)*(ppi+ppip)
                    enddo
                    ft1 = ft1 + dcoe_gpu(l,lmax)*uvs1
                    ft2 = ft2 + dcoe_gpu(l,lmax)*uvs2
                    ft3 = ft3 + dcoe_gpu(l,lmax)*uvs3
                    ft4 = ft4 + dcoe_gpu(l,lmax)*uvs4
                    ft5 = ft5 + dcoe_gpu(l,lmax)*uvs5
                    ft6 = ft6 + dcoe_gpu(l,lmax)*uvs6
                enddo
                fhat_gpu(i,j,k,1) = 0.25_mykind*ft1
                fhat_gpu(i,j,k,2) = 0.25_mykind*ft2
                fhat_gpu(i,j,k,3) = 0.25_mykind*ft3
                fhat_gpu(i,j,k,4) = 0.25_mykind*ft4 + 0.5_mykind*ft6
                fhat_gpu(i,j,k,5) = 0.25_mykind*ft5
            else
!
                kp = k + 1
!
!               Compute Roe average
!
!               Left state (node i)
                ri        =  1._mykind/w_gpu(i,j,k,1)
                uu        =  w_gpu(i,j,k,2) * ri
                vv        =  w_gpu(i,j,k,3) * ri
                ww        =  w_gpu(i,j,k,4) * ri
                qq        =  0.5_mykind * (uu*uu  +vv*vv + ww*ww)
!               pp        =  gm1 * (w_gpu(i,j,k,5) - w_gpu(i,j,k,1) * qq)
!               h         =  (w_gpu(i,j,k,5)  + pp) * ri
                h         =  gamma*w_gpu(i,j,k,5)*ri-gm1*qq
!               Right state (node i+1)
                rip       =  1._mykind/w_gpu(i,j,kp,1)
                up        =  w_gpu(i,j,kp,2) * rip
                vp        =  w_gpu(i,j,kp,3) * rip
                wp        =  w_gpu(i,j,kp,4) * rip
                qqp       =  0.5_mykind * (up*up  +vp*vp +wp*wp)
!               ppp       =  gm1 * (w_gpu(i,jp,k,5) - w_gpu(i,jp,k,1) * qqp)
!               hp        =  (w_gpu(i,jp,k,5)  + ppp) * rip
                hp        =  gamma*w_gpu(i,j,kp,5)*rip-gm1*qqp
!               average state
                r         =  w_gpu(i,j,kp,1)*ri
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
!               left eigenvectors matrix (at Roe state)
!
                b2     =   gm1/cc
                b1     =   b2 * qq
                el(1,1) =   0.5_mykind * (b1      + ww* ci)
                el(2,1) =  -0.5_mykind * (b2 * uu         )
                el(3,1) =  -0.5_mykind * (b2 * vv         )
                el(4,1) =  -0.5_mykind * (b2 * ww +     ci)
                el(5,1) =   0.5_mykind * b2
                el(1,2) =   1._mykind - b1
                el(2,2) =   b2 * uu
                el(3,2) =   b2 * vv
                el(4,2) =   b2 * ww
                el(5,2) =  -b2
                el(1,3) =   0.5_mykind * (b1     - ww * ci)
                el(2,3) =  -0.5_mykind * (b2 * uu         )
                el(3,3) =  -0.5_mykind * (b2 * vv         )
                el(4,3) =  -0.5_mykind * (b2 * ww -     ci)
                el(5,3) =   0.5_mykind * b2
                el(1,4) =  -uu
                el(2,4) =   1._mykind
                el(3,4) =   0._mykind
                el(4,4) =   0._mykind
                el(5,4) =   0._mykind
                el(1,5) =  -vv
                el(2,5) =   0._mykind
                el(3,5) =   1._mykind
                el(4,5) =   0._mykind
                el(5,5) =   0._mykind
!
!               right eigenvectors matrix (at Roe state)
!
                er(1,1) =  1._mykind
                er(2,1) =  1._mykind
                er(3,1) =  1._mykind
                er(4,1) =  0._mykind
                er(5,1) =  0._mykind
                er(1,2) =  uu
                er(2,2) =  uu
                er(3,2) =  uu
                er(4,2) =  1._mykind
                er(5,2) =  0._mykind
                er(1,3) =  vv
                er(2,3) =  vv
                er(3,3) =  vv
                er(4,3) =  0._mykind
                er(5,3) =  1._mykind
                er(1,4) =  ww - c
                er(2,4) =  ww
                er(3,4) =  ww + c
                er(4,4) =  0._mykind
                er(5,4) =  0._mykind
                er(1,5) =  h  - ww * c
                er(2,5) =  qq
                er(3,5) =  h  + ww * c
                er(4,5) =  uu
                er(5,5) =  vv

                do m=1,5 ! loop on characteristic fields
                    evmax(m) = -1._mykind
                enddo
                do l=1,ng2 ! LLF
                    ll = k + l - iweno
!
                    rho  = w_gpu(i,j,ll,1)
                    rhow = w_gpu(i,j,ll,4)
                    ri   = 1._mykind/rho
                    ww   = rhow*ri
                    pp   = rho*temperature_gpu(i,j,ll)
                    gpr  = gamma * pp * ri
                    c    = sqrt (gpr)
                    ev(1) = abs(ww-c)
                    ev(2) = abs(ww)
                    ev(3) = abs(ww+c)
                    ev(4) = ev(2)
                    ev(5) = ev(2)
                    do m=1,5
                        evmax(m) = max(ev(m),evmax(m))
                    enddo
                enddo
                do l=1,ng2 ! loop over the stencil centered at face i
                    ll = k + l - iweno

                    rho  = w_gpu(i,j,ll,1)
                    rhow = w_gpu(i,j,ll,4)
                    ri   = 1._mykind/rho
                    ww   = rhow*ri
                    pp   = rho*temperature_gpu(i,j,ll)
                    fi(1)  =       w_gpu(i,j,ll,4)
                    fi(2)  = ww *  w_gpu(i,j,ll,2)
                    fi(3)  = ww *  w_gpu(i,j,ll,3)
                    fi(4)  = ww *  w_gpu(i,j,ll,4)  + pp
                    fi(5)  = ww * (w_gpu(i,j,ll,5)  + pp)
                    do m=1,5
                        wc = 0._mykind
                        gc = 0._mykind

                        do mm=1,5
                            wc = wc + el(mm,m) * w_gpu(i,j,ll,mm)
                            gc = gc + el(mm,m) * fi(mm)
                        enddo
                        gplus (m,l,i,j) = 0.5_mykind * (gc + evmax(m) * wc)
                        gminus(m,l,i,j) = gc - gplus(m,l,i,j)
                    enddo
                enddo
!
!               Reconstruction of the '+' and '-' fluxes
!
                call wenorec(5,gplus(:,:,i,j),gminus(:,:,i,j),gl,gr,iweno)
!
                do m=1,5
                    ghat(m) = gl(m) + gr(m) ! char. flux
                enddo

!               !Return to conservative fluxes
                do m=1,5
                    fhat_gpu(i,j,k,m) = 0._mykind
                    do mm=1,5
                       fhat_gpu(i,j,k,m) = fhat_gpu(i,j,k,m) + er(mm,m) * ghat(mm)
                    enddo
                enddo
            endif
        enddo

!       Update net flux 
        do k=kstart,kend ! loop on the inner nodes
            do m=1,5
                df = (fhat_gpu(i,j,k,m)-fhat_gpu(i,j,k-1,m))*dzitdz_gpu(k)
                fl_gpu(i,j,k,m) = fl_gpu(i,j,k,m) + df
            enddo
        enddo

    endsubroutine euler_k_kernel
#endif

#ifdef USE_CUDA
attributes(device) &
#endif
subroutine wenorec(nvar,vp,vm,vminus,vplus,iweno)
!    
     implicit none
     integer, parameter :: mykind = MYKIND
!
!    Passed arguments
     integer :: nvar, iweno
     real(mykind),dimension(nvar,2*iweno) :: vm,vp
     real(mykind),dimension(nvar) :: vminus,vplus
!    
!    Local variables
     real(mykind),dimension(-1:4) :: dwe           ! linear weights
     real(mykind),dimension(-1:4) :: alfp,alfm     ! alpha_l
     real(mykind),dimension(-1:4) :: alfp_map,alfm_map ! alpha_l
     real(mykind),dimension(-1:4) :: betap,betam   ! beta_l
     real(mykind),dimension(-1:4) :: omp,omm       ! WENO weights
!    
     integer :: r,i,j,k,l,m
     real(mykind) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4,summ,sump
     real(mykind) :: x,y,y2
!    
     if (iweno==1) then ! Godunov
!    
         i = iweno ! index of intermediate node to perform reconstruction
!    
         vminus(1:nvar) = vp(1:nvar,i)
         vplus (1:nvar) = vm(1:nvar,i+1)
!    
     elseif (iweno==2) then ! WENO-3
!    
         i = iweno ! index of intermediate node to perform reconstruction
!    
         dwe(1)   = 2._mykind/3._mykind
         dwe(0)   = 1._mykind/3._mykind
!    
         do m=1,nvar
!    
             betap(0)  = (vp(m,i  )-vp(m,i-1))**2
             betap(1)  = (vp(m,i+1)-vp(m,i  ))**2
             betam(0)  = (vm(m,i+2)-vm(m,i+1))**2
             betam(1)  = (vm(m,i+1)-vm(m,i  ))**2
!    
             sump = 0._mykind
             summ = 0._mykind
             do l=0,1
                 alfp(l) = dwe(l)/(0.000001_mykind+betap(l))**2
                 alfm(l) = dwe(l)/(0.000001_mykind+betam(l))**2
                 sump = sump + alfp(l)
                 summ = summ + alfm(l)
             enddo
             do l=0,1
                 omp(l) = alfp(l)/sump
                 omm(l) = alfm(l)/summ
             enddo
!    
             vminus(m) = omp(0) *(-vp(m,i-1)+3*vp(m,i  )) + omp(1) *( vp(m,i  )+ vp(m,i+1))
             vplus(m)  = omm(0) *(-vm(m,i+2)+3*vm(m,i+1)) + omm(1) *( vm(m,i  )+ vm(m,i+1))
!    
         enddo ! end of m-loop
!    
         do m=1,nvar
             vminus(m) = 0.5_mykind*vminus(m)
             vplus(m)  = 0.5_mykind*vplus(m)
         enddo
!    
     elseif (iweno==3) then ! WENO-5
!    
      i = iweno ! index of intermediate node to perform reconstruction
!    
      dwe( 0) = 1._mykind/10._mykind
      dwe( 1) = 6._mykind/10._mykind
      dwe( 2) = 3._mykind/10._mykind
!     JS
      d0 = 13._mykind/12._mykind
      d1 = 1._mykind/4._mykind
!     Weights for polynomial reconstructions
      c0 = 1._mykind/3._mykind
      c1 = 5._mykind/6._mykind
      c2 =-1._mykind/6._mykind
      c3 =-7._mykind/6._mykind
      c4 =11._mykind/6._mykind
!    
      do m=1,nvar
!    
       betap(2) = d0*(     vp(m,i)-2._mykind*vp(m,i+1)+vp(m,i+2))**2+d1*(3._mykind*vp(m,i)-4._mykind*vp(m,i+1)+vp(m,i+2))**2
       betap(1) = d0*(     vp(m,i-1)-2._mykind*vp(m,i)+vp(m,i+1))**2+d1*(     vp(m,i-1)-vp(m,i+1) )**2
       betap(0) = d0*(     vp(m,i)-2._mykind*vp(m,i-1)+vp(m,i-2))**2+d1*(3._mykind*vp(m,i)-4._mykind*vp(m,i-1)+vp(m,i-2))**2
!    
       betam(2) = d0*(     vm(m,i+1)-2._mykind*vm(m,i)+vm(m,i-1))**2+d1*(3._mykind*vm(m,i+1)-4._mykind*vm(m,i)+vm(m,i-1))**2
       betam(1) = d0*(     vm(m,i+2)-2._mykind*vm(m,i+1)+vm(m,i))**2+d1*(     vm(m,i+2)-vm(m,i) )**2
       betam(0) = d0*(     vm(m,i+1)-2._mykind*vm(m,i+2)+vm(m,i+3))**2+d1*(3._mykind*vm(m,i+1)-4._mykind*vm(m,i+2)+vm(m,i+3))**2
!    
       sump = 0._mykind
       summ = 0._mykind
       do l=0,2
        alfp(l) = dwe(  l)/(0.000001_mykind+betap(l))**2
        alfm(l) = dwe(  l)/(0.000001_mykind+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,2
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m)   = omp(2)*(c0*vp(m,i  )+c1*vp(m,i+1)+c2*vp(m,i+2)) + &
         & omp(1)*(c2*vp(m,i-1)+c1*vp(m,i  )+c0*vp(m,i+1)) + omp(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i  ))
       vplus(m)   = omm(2)*(c0*vm(m,i+1)+c1*vm(m,i  )+c2*vm(m,i-1)) +  &
         & omm(1)*(c2*vm(m,i+2)+c1*vm(m,i+1)+c0*vm(m,i  )) + omm(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
!    
      enddo ! end of m-loop 
!    
     elseif (iweno==4) then ! WENO-7
!    
      i = iweno ! index of intermediate node to perform reconstruction
!    
      dwe( 0) = 1._mykind/35._mykind
      dwe( 1) = 12._mykind/35._mykind
      dwe( 2) = 18._mykind/35._mykind
      dwe( 3) = 4._mykind/35._mykind
!    
!     JS weights
      d1 = 1._mykind/36._mykind
      d2 = 13._mykind/12._mykind
      d3 = 781._mykind/720._mykind
!    
      do m=1,nvar
!    
       betap(3)= d1*(-11*vp(m,  i)+18*vp(m,i+1)- 9*vp(m,i+2)+ 2*vp(m,i+3))**2+&
       &  d2*(  2*vp(m,  i)- 5*vp(m,i+1)+ 4*vp(m,i+2)-   vp(m,i+3))**2+ &
       & d3*(   -vp(m,  i)+ 3*vp(m,i+1)- 3*vp(m,i+2)+   vp(m,i+3))**2
       betap(2)= d1*(- 2*vp(m,i-1)- 3*vp(m,i  )+ 6*vp(m,i+1)-   vp(m,i+2))**2+&
       &  d2*(    vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1)             )**2+&
       &  d3*(   -vp(m,i-1)+ 3*vp(m,i  )- 3*vp(m,i+1)+   vp(m,i+2))**2
       betap(1)= d1*(    vp(m,i-2)- 6*vp(m,i-1)+ 3*vp(m,i  )+ 2*vp(m,i+1))**2+&
       &  d2*( vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1))**2+ &
       &  d3*(   -vp(m,i-2)+ 3*vp(m,i-1)- 3*vp(m,i  )+   vp(m,i+1))**2
       betap(0)= d1*(- 2*vp(m,i-3)+ 9*vp(m,i-2)-18*vp(m,i-1)+11*vp(m,i  ))**2+&
       &  d2*(-   vp(m,i-3)+ 4*vp(m,i-2)- 5*vp(m,i-1)+ 2*vp(m,i  ))**2+&
       &  d3*(   -vp(m,i-3)+ 3*vp(m,i-2)- 3*vp(m,i-1)+   vp(m,i  ))**2
!    
       betam(3)= d1*(-11*vm(m,i+1)+18*vm(m,i  )- 9*vm(m,i-1)+ 2*vm(m,i-2))**2+&
       &  d2*(  2*vm(m,i+1)- 5*vm(m,i  )+ 4*vm(m,i-1)-   vm(m,i-2))**2+&
       &  d3*(   -vm(m,i+1)+ 3*vm(m,i  )- 3*vm(m,i-1)+   vm(m,i-2))**2
       betam(2)= d1*(- 2*vm(m,i+2)- 3*vm(m,i+1)+ 6*vm(m,i  )-   vm(m,i-1))**2+&
       &  d2*(    vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  )             )**2+&
       &  d3*(   -vm(m,i+2)+ 3*vm(m,i+1)- 3*vm(m,i  )+   vm(m,i-1))**2
       betam(1)= d1*(    vm(m,i+3)- 6*vm(m,i+2)+ 3*vm(m,i+1)+ 2*vm(m,i  ))**2+&
       &  d2*(                 vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  ))**2+&
       &  d3*(   -vm(m,i+3)+ 3*vm(m,i+2)- 3*vm(m,i+1)+   vm(m,i  ))**2
       betam(0)= d1*(- 2*vm(m,i+4)+ 9*vm(m,i+3)-18*vm(m,i+2)+11*vm(m,i+1))**2+&
       &  d2*(-   vm(m,i+4)+ 4*vm(m,i+3)- 5*vm(m,i+2)+ 2*vm(m,i+1))**2+&
       &  d3*(   -vm(m,i+4)+ 3*vm(m,i+3)- 3*vm(m,i+2)+   vm(m,i+1))**2 
!    
       sump = 0._mykind
       summ = 0._mykind
       do l=0,3
        alfp(l) = dwe(  l)/(0.000001_mykind+betap(l))**2
        alfm(l) = dwe(  l)/(0.000001_mykind+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,3
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!    
       vminus(m)   = omp(3)*( 6*vp(m,i  )+26*vp(m,i+1)-10*vp(m,i+2)+ 2*vp(m,i+3))+&
        omp(2)*(-2*vp(m,i-1)+14*vp(m,i  )+14*vp(m,i+1)- 2*vp(m,i+2))+&
        omp(1)*( 2*vp(m,i-2)-10*vp(m,i-1)+26*vp(m,i  )+ 6*vp(m,i+1))+&
        omp(0)*(-6*vp(m,i-3)+26*vp(m,i-2)-46*vp(m,i-1)+50*vp(m,i  ))
       vplus(m)   =  omm(3)*( 6*vm(m,i+1)+26*vm(m,i  )-10*vm(m,i-1)+ 2*vm(m,i-2))+&
        omm(2)*(-2*vm(m,i+2)+14*vm(m,i+1)+14*vm(m,i  )- 2*vm(m,i-1))+&
        omm(1)*( 2*vm(m,i+3)-10*vm(m,i+2)+26*vm(m,i+1)+ 6*vm(m,i  ))+&
        omm(0)*(-6*vm(m,i+4)+26*vm(m,i+3)-46*vm(m,i+2)+50*vm(m,i+1))
!    
      enddo ! end of m-loop 
!    
      vminus = vminus/24._mykind
      vplus  = vplus /24._mykind
!    
     else
      write(*,*) 'Error! WENO scheme not implemented'
      stop
     endif

endsubroutine wenorec
!
#ifdef USE_CUDA
attributes(device) &
#endif
subroutine wenorecsymbo(nvar,vp,vm,vminus,vplus,iweno)
!    
     implicit none
     integer, parameter :: mykind = MYKIND
!
!    Passed arguments
     integer :: nvar, iweno
     real(mykind),dimension(nvar,2*iweno) :: vm,vp
     real(mykind),dimension(nvar) :: vminus,vplus
!    
!    Local variables
     real(mykind),dimension(-1:4) :: dwe           ! linear weights
     real(mykind),dimension(-1:4) :: alfp,alfm     ! alpha_l
     real(mykind),dimension(-1:4) :: betap,betam   ! beta_l
     real(mykind),dimension(-1:4) :: omp,omm       ! WENO weights
!    
     integer :: r,i,j,k,l,m
     real(mykind) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4,summ,sump
!    
     if (iweno==3) then
!    
      i = iweno ! index of intermediate node to perform reconstruction
!    
      dwe( 0) = 0.094647545896
      dwe( 1) = 0.428074212384
      dwe( 2) = 0.408289331408
      dwe( 3) = 0.068988910311
!     JS
      d0 = 13._mykind/12._mykind
      d1 = 1._mykind/4._mykind
!     Weights for polynomial reconstructions
      c0 = 1._mykind/3._mykind
      c1 = 5._mykind/6._mykind
      c2 =-1._mykind/6._mykind
      c3 =-7._mykind/6._mykind
      c4 =11._mykind/6._mykind
!    
      do m=1,nvar
!    
       betap(3) = d0*(vp(m,i+1)-2._mykind*vp(m,i+2)+vp(m,i+3))**2+d1*(5._mykind*vp(m,i+1)-8._mykind*vp(m,i+2)+&
       3._mykind*vp(m,i+3))**2
       betap(2) = d0*(vp(m,i)  -2._mykind*vp(m,i+1)+vp(m,i+2))**2+d1*(3._mykind*vp(m,i)  -4._mykind*vp(m,i+1)+vp(m,i+2))**2
       betap(1) = d0*(vp(m,i-1)-2._mykind*vp(m,i)  +vp(m,i+1))**2+d1*(          vp(m,i-1)-vp(m,i+1) )**2
       betap(0) = d0*(vp(m,i)  -2._mykind*vp(m,i-1)+vp(m,i-2))**2+d1*(3._mykind*vp(m,i)-4._mykind*vp(m,i-1)+vp(m,i-2))**2
       betap(3) = max(betap(0),betap(1),betap(2),betap(3))
! 
       betam(3) = d0*(vm(m,i  )-2._mykind*vm(m,i-1)+vm(m,i-2))**2+d1*(5._mykind*vm(m,i  )-8._mykind*vm(m,i-1)+&
       3._mykind*vm(m,i-2))**2
       betam(2) = d0*(vm(m,i+1)-2._mykind*vm(m,i)  +vm(m,i-1))**2+d1*(3._mykind*vm(m,i+1)-4._mykind*vm(m,i)+vm(m,i-1))**2
       betam(1) = d0*(vm(m,i+2)-2._mykind*vm(m,i+1)+vm(m,i  ))**2+d1*(          vm(m,i+2)-          vm(m,i) )**2
       betam(0) = d0*(vm(m,i+1)-2._mykind*vm(m,i+2)+vm(m,i+3))**2+d1*(3._mykind*vm(m,i+1)-4._mykind*vm(m,i+2)+vm(m,i+3))**2
       betam(3) = max(betam(0),betam(1),betam(2),betam(3))
!
       sump  = 0._mykind
       summ  = 0._mykind
       do l=0,3
        alfp(l) = dwe(l)/(0.000000000001_mykind+betap(l))**2
        alfm(l) = dwe(l)/(0.000000000001_mykind+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,3
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m) = omp(3)*(c4*vp(m,i+1)+c3*vp(m,i+2)+c0*vp(m,i+3)) + &
                   omp(2)*(c0*vp(m,i  )+c1*vp(m,i+1)+c2*vp(m,i+2)) + &
                   omp(1)*(c2*vp(m,i-1)+c1*vp(m,i  )+c0*vp(m,i+1)) + &
                   omp(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i  ))
       vplus(m)  = omm(3)*(c4*vm(m,i  )+c3*vm(m,i-1)+c0*vm(m,i-2)) + &
                   omm(2)*(c0*vm(m,i+1)+c1*vm(m,i  )+c2*vm(m,i-1)) + &
                   omm(1)*(c2*vm(m,i+2)+c1*vm(m,i+1)+c0*vm(m,i  )) + &
                   omm(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
!    
      enddo ! end of m-loop 
!    
     endif
!
endsubroutine wenorecsymbo
!
#ifdef USE_CUDA
attributes(device) &
#endif
subroutine upwrec(nvar,vp,vm,vminus,vplus,iweno)
!    
     implicit none
     integer, parameter :: mykind = MYKIND
!
!    Passed arguments
     integer :: nvar, iweno
     real(mykind),dimension(nvar,2*iweno) :: vm,vp
     real(mykind),dimension(nvar) :: vminus,vplus
!    
!    Local variables
     real(mykind),dimension(-1:4) :: dwe           ! linear weights
     real(mykind),dimension(-1:4) :: alfp,alfm     ! alpha_l
     real(mykind),dimension(-1:4) :: betap,betam   ! beta_l
     real(mykind),dimension(-1:4) :: omp,omm       ! WENO weights
!    
     integer :: r,i,j,k,l,m
     real(mykind) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4,summ,sump
!    
     if (iweno==1) then ! Godunov
!    
      i = iweno ! index of intermediate node to perform reconstruction
!    
      vminus(1:nvar) = vp(1:nvar,i)
      vplus (1:nvar) = vm(1:nvar,i+1)
!    
     elseif (iweno==2) then ! WENO-3
!    
      i = iweno ! index of intermediate node to perform reconstruction
!
      vminus(1:nvar) =  -1*vp(1:nvar,i-1)+5*vp(1:nvar,i)  +2*vp(1:nvar,i+1)
      vplus (1:nvar) =   2*vm(1:nvar,i  )+5*vm(1:nvar,i+1)-1*vm(1:nvar,i+2)
      vminus(1:nvar) = vminus(1:nvar)/6._mykind
      vplus (1:nvar) = vplus (1:nvar)/6._mykind
!
     elseif (iweno==3) then ! WENO-5
!    
      i = iweno ! index of intermediate node to perform reconstruction
!
      vminus(1:nvar) =  2*vp(1:nvar,i-2)-13*vp(1:nvar,i-1)+47*vp(1:nvar,i)  +27*vp(1:nvar,i+1)-3*vp(1:nvar,i+2)
      vplus (1:nvar) = -3*vm(1:nvar,i-1)+27*vm(1:nvar,i  )+47*vm(1:nvar,i+1)-13*vm(1:nvar,i+2)+2*vm(1:nvar,i+3)
      vminus(1:nvar) = vminus(1:nvar)/60._mykind
      vplus (1:nvar) = vplus (1:nvar)/60._mykind
!    
     endif
!
end subroutine upwrec
!
#ifdef USE_CUDA
attributes(device) &
#endif
subroutine weno5map(nvar,vp,vm,vminus,vplus,iweno)
!    
     implicit none
     integer, parameter :: mykind = MYKIND
!
!    Passed arguments
     integer :: nvar, iweno
     real(mykind),dimension(nvar,2*iweno) :: vm,vp
     real(mykind),dimension(nvar) :: vminus,vplus
!    
!    Local variables
     real(mykind),dimension(-1:4) :: dwe           ! linear weights
     real(mykind),dimension(-1:4) :: alfp,alfm     ! alpha_l
     real(mykind),dimension(-1:4) :: alfp_map,alfm_map ! alpha_l
     real(mykind),dimension(-1:4) :: betap,betam   ! beta_l
     real(mykind),dimension(-1:4) :: omp,omm       ! WENO weights
!    
     integer :: r,i,j,k,l,m
     real(mykind) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4,summ,sump
     real(mykind) :: x,y,y2
!    
     if (iweno==3) then ! WENO-5 mapped
!    
      i = iweno ! index of intermediate node to perform reconstruction
!    
      dwe( 0) = 1._mykind/10._mykind
      dwe( 1) = 6._mykind/10._mykind
      dwe( 2) = 3._mykind/10._mykind
!     JS
      d0 = 13._mykind/12._mykind
      d1 = 1._mykind/4._mykind
!     Weights for polynomial reconstructions
      c0 = 1._mykind/3._mykind
      c1 = 5._mykind/6._mykind
      c2 =-1._mykind/6._mykind
      c3 =-7._mykind/6._mykind
      c4 =11._mykind/6._mykind
!    
      do m=1,nvar
!    
       betap(2) = d0*(     vp(m,i)-2._mykind*vp(m,i+1)+vp(m,i+2))**2+d1*(3._mykind*vp(m,i)-4._mykind*vp(m,i+1)+vp(m,i+2))**2
       betap(1) = d0*(     vp(m,i-1)-2._mykind*vp(m,i)+vp(m,i+1))**2+d1*(     vp(m,i-1)-vp(m,i+1) )**2
       betap(0) = d0*(     vp(m,i)-2._mykind*vp(m,i-1)+vp(m,i-2))**2+d1*(3._mykind*vp(m,i)-4._mykind*vp(m,i-1)+vp(m,i-2))**2
!    
       betam(2) = d0*(     vm(m,i+1)-2._mykind*vm(m,i)+vm(m,i-1))**2+d1*(3._mykind*vm(m,i+1)-4._mykind*vm(m,i)+vm(m,i-1))**2
       betam(1) = d0*(     vm(m,i+2)-2._mykind*vm(m,i+1)+vm(m,i))**2+d1*(     vm(m,i+2)-vm(m,i) )**2
       betam(0) = d0*(     vm(m,i+1)-2._mykind*vm(m,i+2)+vm(m,i+3))**2+d1*(3._mykind*vm(m,i+1)-4._mykind*vm(m,i+2)+vm(m,i+3))**2
!    
       sump = 0._mykind
       summ = 0._mykind
       do l=0,2
        alfp(l) = dwe(  l)/(0.000001_mykind+betap(l))**2
        alfm(l) = dwe(  l)/(0.000001_mykind+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,2
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
!********************************************
!      Mapping procedure
       do l=0,2
        x  = omp(l)
        y  = dwe(l)
        y2 = y*y
        alfp_map(l) = x*(y+y2-3._mykind*y*x+x*x)/(y2+x*(1._mykind-2._mykind*y))
        x = omm(l)
        alfm_map(l) = x*(y+y2-3._mykind*y*x+x*x)/(y2+x*(1._mykind-2._mykind*y))
       enddo
!
       sump = 0._mykind
       summ = 0._mykind
       do l=0,2
        sump = sump + alfp_map(l)
        summ = summ + alfm_map(l)
       enddo
       do l=0,2
        omp(l) = alfp_map(l)/sump
        omm(l) = alfm_map(l)/summ
       enddo
!*********************************************
       vminus(m)   = omp(2)*(c0*vp(m,i  )+c1*vp(m,i+1)+c2*vp(m,i+2)) + &
         & omp(1)*(c2*vp(m,i-1)+c1*vp(m,i  )+c0*vp(m,i+1)) + omp(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i  ))
       vplus(m)   = omm(2)*(c0*vm(m,i+1)+c1*vm(m,i  )+c2*vm(m,i-1)) +  &
         & omm(1)*(c2*vm(m,i+2)+c1*vm(m,i+1)+c0*vm(m,i  )) + omm(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
!    
      enddo ! end of m-loop 
!    
     endif
!
endsubroutine weno5map
!
#ifdef USE_CUDA
attributes(device) &
#endif
subroutine wenocu6(nvar,vp,vm,vminus,vplus,iweno)
!    
     implicit none
     integer, parameter :: mykind = MYKIND
!
!    Passed arguments
     integer :: nvar, iweno
     real(mykind),dimension(nvar,2*iweno) :: vm,vp
     real(mykind),dimension(nvar) :: vminus,vplus
!    
!    Local variables
     real(mykind),dimension(-1:4) :: dwe           ! linear weights
     real(mykind),dimension(-1:4) :: alfp,alfm     ! alpha_l
     real(mykind),dimension(-1:4) :: betap,betam   ! beta_l
     real(mykind),dimension(-1:4) :: omp,omm       ! WENO weights
!    
     integer :: r,i,j,k,l,m
     real(mykind) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4,summ,sump
     real(mykind) :: bigc,eps40,tau6p,tau6m
!    
     if (iweno==3) then ! WENO-CU6
!    
      i = iweno ! index of intermediate node to perform reconstruction
!    
      dwe( 0) = 1._mykind/20._mykind
      dwe( 1) = 9._mykind/20._mykind
      dwe( 2) = 9._mykind/20._mykind
      dwe( 3) = 1._mykind/20._mykind
!     JS
      d0 = 13._mykind/12._mykind
      d1 = 1._mykind/4._mykind
!     Weights for polynomial reconstructions
      c0 = 1._mykind/3._mykind
      c1 = 5._mykind/6._mykind
      c2 =-1._mykind/6._mykind
      c3 =-7._mykind/6._mykind
      c4 =11._mykind/6._mykind
!    
      do m=1,nvar
!    
       betap(3) = 271779*vp(m,i-2)**2+&
                         vp(m,i-2)*(2380800*vp(m,i-1)+4086352 *vp(m,i)-3462252 *vp(m,i+1)+1458762 *vp(m,i+2)-245620 *vp(m,i+3))+&
                         vp(m,i-1)*(5653317*vp(m,i-1)-20427884*vp(m,i)+17905032*vp(m,i+1)-7727988 *vp(m,i+2)+1325006*vp(m,i+3))+&
                         vp(m,i  )*(                 +19510972*vp(m,i)-35817664*vp(m,i+1)+15929912*vp(m,i+2)-2792660*vp(m,i+3))+&
                         vp(m,i+1)*(                                 +17195652 *vp(m,i+1)-15880404*vp(m,i+2)+2863984*vp(m,i+3))+&
                         vp(m,i+2)*(                                                      +3824847*vp(m,i+2)-1429976*vp(m,i+3))+&
                  139633*vp(m,i+3)**2
       betap(3) = betap(3)/120960._mykind
       betap(2) = d0*(     vp(m,i)  -2._mykind*vp(m,i+1)+vp(m,i+2))**2+d1*(3._mykind*vp(m,i)  -4._mykind*vp(m,i+1)+vp(m,i+2))**2
       betap(1) = d0*(     vp(m,i-1)-2._mykind*vp(m,i)+vp(m,i+1))**2+d1*(     vp(m,i-1)-vp(m,i+1) )**2
       betap(0) = d0*(     vp(m,i)-2._mykind*vp(m,i-1)+vp(m,i-2))**2+d1*(3._mykind*vp(m,i)-4._mykind*vp(m,i-1)+vp(m,i-2))**2
 
       betam(3) = 271779*vm(m,i+3)**2+&
                         vm(m,i+3)*(2380800*vm(m,i+2)+4086352 *vm(m,i+1)-3462252 *vm(m,i)+1458762 *vm(m,i-1)-245620 *vm(m,i-2))+&
                         vm(m,i+2)*(5653317*vm(m,i+2)-20427884*vm(m,i+1)+17905032*vm(m,i)-7727988 *vm(m,i-1)+1325006*vm(m,i-2))+&
                         vm(m,i+1)*(                 +19510972*vm(m,i+1)-35817664*vm(m,i)+15929912*vm(m,i-1)-2792660*vm(m,i-2))+&
                         vm(m,i  )*(                                   +17195652 *vm(m,i)-15880404*vm(m,i-1)+2863984*vm(m,i-2))+&
                         vm(m,i-1)*(                                                      +3824847*vm(m,i-1)-1429976*vm(m,i-2))+&
                  139633*vm(m,i-2)**2
       betam(3) = betam(3)/120960._mykind
       betam(2) = d0*(     vm(m,i+1)-2._mykind*vm(m,i)  +vm(m,i-1))**2+d1*(3._mykind*vm(m,i+1)-4._mykind*vm(m,i)+vm(m,i-1))**2
       betam(1) = d0*(     vm(m,i+2)-2._mykind*vm(m,i+1)+vm(m,i  ))**2+d1*(          vm(m,i+2)-          vm(m,i) )**2
       betam(0) = d0*(     vm(m,i+1)-2._mykind*vm(m,i+2)+vm(m,i+3))**2+d1*(3._mykind*vm(m,i+1)-4._mykind*vm(m,i+2)+vm(m,i+3))**2
!
       tau6p = betap(3)-(betap(0)+betap(2)+4*betap(1))/6._mykind
       tau6m = betam(3)-(betam(0)+betam(2)+4*betam(1))/6._mykind
!    
       bigc  = 20._mykind
       eps40 = 1.D-40
       sump  = 0._mykind
       summ  = 0._mykind
       do l=0,3
        alfp(l) = dwe(l)*(bigc+tau6p/(betap(l)+eps40))
        alfm(l) = dwe(l)*(bigc+tau6m/(betam(l)+eps40))
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,3
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m) = omp(3)*(c4*vp(m,i+1)+c3*vp(m,i+2)+c0*vp(m,i+3)) + &
                   omp(2)*(c0*vp(m,i  )+c1*vp(m,i+1)+c2*vp(m,i+2)) + &
                   omp(1)*(c2*vp(m,i-1)+c1*vp(m,i  )+c0*vp(m,i+1)) + &
                   omp(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i  ))
       vplus(m)  = omm(3)*(c4*vm(m,i  )+c3*vm(m,i-1)+c0*vm(m,i-2)) + &
                   omm(2)*(c0*vm(m,i+1)+c1*vm(m,i  )+c2*vm(m,i-1)) + &
                   omm(1)*(c2*vm(m,i+2)+c1*vm(m,i+1)+c0*vm(m,i  )) + &
                   omm(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
!    
      enddo ! end of m-loop 
!    
     endif
!
endsubroutine wenocu6
!
#ifdef USE_CUDA
attributes(device) &
#endif
subroutine weno5z(nvar,vp,vm,vminus,vplus,iweno)
!    
     implicit none
     integer, parameter :: mykind = MYKIND
!
!    Passed arguments
     integer :: nvar, iweno
     real(mykind),dimension(nvar,2*iweno) :: vm,vp
     real(mykind),dimension(nvar) :: vminus,vplus
!    
!    Local variables
     real(mykind),dimension(-1:4) :: dwe           ! linear weights
     real(mykind),dimension(-1:4) :: alfp,alfm     ! alpha_l
     real(mykind),dimension(-1:4) :: alfp_map,alfm_map ! alpha_l
     real(mykind),dimension(-1:4) :: betap,betam   ! beta_l
     real(mykind),dimension(-1:4) :: betazp,betazm ! betaz_l
     real(mykind),dimension(-1:4) :: omp,omm       ! WENO weights
!    
     integer :: r,i,j,k,l,m
     real(mykind) :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4,summ,sump,eps40,tau5p,tau5m
!    
     if (iweno==3) then ! WENO-5 mapped
!    
      i = iweno ! index of intermediate node to perform reconstruction
!    
      dwe( 0) = 1._mykind/10._mykind
      dwe( 1) = 6._mykind/10._mykind
      dwe( 2) = 3._mykind/10._mykind
!     JS
      d0 = 13._mykind/12._mykind
      d1 = 1._mykind/4._mykind
!     Weights for polynomial reconstructions
      c0 = 1._mykind/3._mykind
      c1 = 5._mykind/6._mykind
      c2 =-1._mykind/6._mykind
      c3 =-7._mykind/6._mykind
      c4 =11._mykind/6._mykind
!    
      do m=1,nvar
!    
       betap(2) = d0*(     vp(m,i)-2._mykind*vp(m,i+1)+vp(m,i+2))**2+d1*(3._mykind*vp(m,i)-4._mykind*vp(m,i+1)+vp(m,i+2))**2
       betap(1) = d0*(     vp(m,i-1)-2._mykind*vp(m,i)+vp(m,i+1))**2+d1*(     vp(m,i-1)-vp(m,i+1) )**2
       betap(0) = d0*(     vp(m,i)-2._mykind*vp(m,i-1)+vp(m,i-2))**2+d1*(3._mykind*vp(m,i)-4._mykind*vp(m,i-1)+vp(m,i-2))**2
!    
       betam(2) = d0*(     vm(m,i+1)-2._mykind*vm(m,i)+vm(m,i-1))**2+d1*(3._mykind*vm(m,i+1)-4._mykind*vm(m,i)+vm(m,i-1))**2
       betam(1) = d0*(     vm(m,i+2)-2._mykind*vm(m,i+1)+vm(m,i))**2+d1*(     vm(m,i+2)-vm(m,i) )**2
       betam(0) = d0*(     vm(m,i+1)-2._mykind*vm(m,i+2)+vm(m,i+3))**2+d1*(3._mykind*vm(m,i+1)-4._mykind*vm(m,i+2)+vm(m,i+3))**2
!
       tau5p = abs(betap(0)-betap(2))
       tau5m = abs(betam(0)-betam(2))
       eps40 = 1.D-40
!
       do l=0,2
        betazp(l) = (betap(l)+eps40)/(betap(l)+eps40+tau5p)
        betazm(l) = (betam(l)+eps40)/(betam(l)+eps40+tau5m)
       enddo
!
       sump = 0._mykind
       summ = 0._mykind
       do l=0,2
        alfp(l) = dwe(l)/betazp(l)
        alfm(l) = dwe(l)/betazm(l)
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,2
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m)   = omp(2)*(c0*vp(m,i  )+c1*vp(m,i+1)+c2*vp(m,i+2)) + &
         & omp(1)*(c2*vp(m,i-1)+c1*vp(m,i  )+c0*vp(m,i+1)) + omp(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i  ))
       vplus(m)   = omm(2)*(c0*vm(m,i+1)+c1*vm(m,i  )+c2*vm(m,i-1)) +  &
         & omm(1)*(c2*vm(m,i+2)+c1*vm(m,i+1)+c0*vm(m,i  )) + omm(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
!    
      enddo ! end of m-loop 
!    
     endif
!
endsubroutine weno5z
!
end module mod_euler
