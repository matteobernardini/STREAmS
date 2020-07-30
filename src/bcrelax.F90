subroutine bcrelax(ilat)
!
! Apply nonreflecting boundary conditions
!
 use mod_streams
 implicit none
!
 real(mykind) :: dwdxi1,dwdxi2,dwdxi3,dwdxi4,dwdxi5
 real(mykind) :: dwdxo1,dwdxo2,dwdxo3,dwdxo4,dwdxo5
 real(mykind) :: dwcdxi1,dwcdxi2,dwcdxi3,dwcdxi4,dwcdxi5
 real(mykind) :: dwcdxo1,dwcdxo2,dwcdxo3,dwcdxo4,dwcdxo5
 real(mykind) :: dwcdx1,dwcdx2,dwcdx3,dwcdx4,dwcdx5
 real(mykind) :: dwdyi1,dwdyi2,dwdyi3,dwdyi4,dwdyi5
 real(mykind) :: dwdyo1,dwdyo2,dwdyo3,dwdyo4,dwdyo5
 real(mykind) :: dwcdyi1,dwcdyi2,dwcdyi3,dwcdyi4,dwcdyi5
 real(mykind) :: dwcdyo1,dwcdyo2,dwcdyo3,dwcdyo4,dwcdyo5
 real(mykind) :: dwcdy1,dwcdy2,dwcdy3,dwcdy4,dwcdy5
 real(mykind) :: wloc0_1,wloc0_2,wloc0_3,wloc0_4,wloc0_5
 real(mykind) :: wloc1_1,wloc1_2,wloc1_3,wloc1_4,wloc1_5
 real(mykind) :: wloc2_1,wloc2_2,wloc2_3,wloc2_4,wloc2_5
 real(mykind) :: wloc3_1,wloc3_2,wloc3_3,wloc3_4,wloc3_5
!
 real(mykind) :: jac_1_1,jac_1_2,jac_1_3,jac_1_4,jac_1_5
 real(mykind) :: jac_2_1,jac_2_2,jac_2_3,jac_2_4,jac_2_5
 real(mykind) :: jac_3_1,jac_3_2,jac_3_3,jac_3_4,jac_3_5
 real(mykind) :: jac_4_1,jac_4_2,jac_4_3,jac_4_4,jac_4_5
 real(mykind) :: jac_5_1,jac_5_2,jac_5_3,jac_5_4,jac_5_5
 real(mykind) :: jacinv_1_1,jacinv_1_2,jacinv_1_3,jacinv_1_4,jacinv_1_5
 real(mykind) :: jacinv_2_1,jacinv_2_2,jacinv_2_3,jacinv_2_4,jacinv_2_5
 real(mykind) :: jacinv_3_1,jacinv_3_2,jacinv_3_3,jacinv_3_4,jacinv_3_5
 real(mykind) :: jacinv_4_1,jacinv_4_2,jacinv_4_3,jacinv_4_4,jacinv_4_5
 real(mykind) :: jacinv_5_1,jacinv_5_2,jacinv_5_3,jacinv_5_4,jacinv_5_5
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
 real(mykind) :: ev1,ev2,ev3,ev4,ev5
 real(mykind) :: cnr1,cnr2,cnr3
!
 integer :: ilat,i,j,k,l,m,ll,mm,ibcord
 real(mykind) :: c,cc,ci,df,p2,pp,qq,rho,ri
 real(mykind) :: uu,vv,ww,xx
!
!Check the bak version for mask
 ibcord = 2
 cnr1   = -1.5_mykind
 cnr2   =  2._mykind
 cnr3   = -0.5_mykind
!
 if (ilat==1) then ! left side
 elseif (ilat==2) then ! right side
!
! Redefining f(u)_x on the boundary
!
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do j=1,ny
!
    i   = nx
!
    rho =  w_gpu(i,j,k,1) 
    ri  =  1._mykind/w_gpu(i,j,k,1)
    uu  =  w_gpu(i,j,k,2) * ri
    vv  =  w_gpu(i,j,k,3) * ri
    ww  =  w_gpu(i,j,k,4) * ri
    qq  =  0.5_mykind * (uu*uu  + vv*vv + ww*ww)
!   pp  =  gm1 * (w_gpu(i,j,k,5) - w_gpu(i,j,k,1) * qq)
    pp  =  rho*temperature_gpu(i,j,k)
    cc  =  gamma * pp * ri
    c   =  sqrt(cc)
    ci  =  1._mykind/c
!   Jacobian of conservative/primitive transformation (Eqn. (A.5) of Lodato et al, JCP 2008)
    jac_1_1  =  1._mykind
    jac_1_2  =  0._mykind
    jac_1_3  =  0._mykind
    jac_1_4  =  0._mykind
    jac_1_5  =  0._mykind
    jac_2_1  =  uu
    jac_2_2  =  rho
    jac_2_3  =  0._mykind
    jac_2_4  =  0._mykind
    jac_2_5  =  0._mykind
    jac_3_1  =  vv
    jac_3_2  =  0._mykind
    jac_3_3  =  rho
    jac_3_4  =  0._mykind
    jac_3_5  =  0._mykind
    jac_4_1  =  ww
    jac_4_2  =  0._mykind
    jac_4_3  =  0._mykind
    jac_4_4  =  rho
    jac_4_5  =  0._mykind
    jac_5_1  =  qq
    jac_5_2  =  w_gpu(i,j,k,2)
    jac_5_3  =  w_gpu(i,j,k,3)
    jac_5_4  =  w_gpu(i,j,k,4)
    jac_5_5  =  gm
!   Jacobian of inverse conservative/primitive transformation (Eqn. (A.5) of Lodato et al, JCP 2008)
    jacinv_1_1 =  1._mykind
    jacinv_1_2 =  0._mykind
    jacinv_1_3 =  0._mykind
    jacinv_1_4 =  0._mykind
    jacinv_1_5 =  0._mykind
    jacinv_2_1 = -uu*ri
    jacinv_2_2 =  ri
    jacinv_2_3 =  0._mykind
    jacinv_2_4 =  0._mykind
    jacinv_2_5 =  0._mykind
    jacinv_3_1 = -vv*ri
    jacinv_3_2 =  0._mykind
    jacinv_3_3 =  ri
    jacinv_3_4 =  0._mykind
    jacinv_3_5 =  0._mykind
    jacinv_4_1 = -ww*ri
    jacinv_4_2 =  0._mykind
    jacinv_4_3 =  0._mykind
    jacinv_4_4 =  ri
    jacinv_4_5 =  0._mykind
    jacinv_5_1 =  gm1*qq
    jacinv_5_2 = -gm1*uu
    jacinv_5_3 = -gm1*vv
    jacinv_5_4 = -gm1*ww
    jacinv_5_5 =  gm1
!   left eigenvectors matrix (Eqn._mykind (A.12_mykind) of Lodato et al, JCP 2008)
    el_1_1 =  0._mykind
    el_1_2 = -rho*c
    el_1_3 =  0._mykind
    el_1_4 =  0._mykind
    el_1_5 =  1._mykind
    el_2_1 =  cc
    el_2_2 =  0._mykind
    el_2_3 =  0._mykind
    el_2_4 =  0._mykind
    el_2_5 = -1._mykind
    el_3_1 =  0._mykind
    el_3_2 =  0._mykind
    el_3_3 =  1._mykind
    el_3_4 =  0._mykind
    el_3_5 =  0._mykind
    el_4_1 =  0._mykind
    el_4_2 =  0._mykind
    el_4_3 =  0._mykind
    el_4_4 =  1._mykind
    el_4_5 =  0._mykind
    el_5_1 =  0._mykind
    el_5_2 =  rho*c
    el_5_3 =  0._mykind
    el_5_4 =  0._mykind
    el_5_5 =  1._mykind
!   left eigenvectors matrix (Eqn._mykind (A.11_mykind) of Lodato et al, JCP 2008)
    er_1_1 =  0.5_mykind/cc
    er_2_1 = -0.5_mykind*ri*ci
    er_3_1 =  0._mykind
    er_4_1 =  0._mykind
    er_5_1 =  0.5_mykind
    er_1_2 =  1._mykind/cc
    er_2_2 =  0._mykind
    er_3_2 =  0._mykind
    er_4_2 =  0._mykind
    er_5_2 =  0._mykind
    er_1_3 =  0._mykind
    er_2_3 =  0._mykind
    er_3_3 =  1._mykind
    er_4_3 =  0._mykind
    er_5_3 =  0._mykind
    er_1_4 =  0._mykind
    er_2_4 =  0._mykind
    er_3_4 =  0._mykind
    er_4_4 =  1._mykind
    er_5_4 =  0._mykind
    er_1_5 =  0.5_mykind/cc
    er_2_5 =  0.5_mykind*ri*ci
    er_3_5 =  0._mykind
    er_4_5 =  0._mykind
    er_5_5 =  0.5_mykind
!
!   Eigenvalues 
    ev1    =  uu-c
    ev2    =  uu
    ev3    =  uu
    ev4    =  uu
    ev5    =  uu+c
!
!   Local array of conservative variables
    wloc0_1 = winf_gpu(1)
    wloc0_2 = winf_gpu(2)
    wloc0_3 = winf_gpu(3)
    wloc0_4 = winf_gpu(4)
    wloc0_5 = winf_gpu(5)
    l=1
    ll = i + 1 - l
    wloc1_1 = w_gpu(ll,j,k,1)
    wloc1_2 = w_gpu(ll,j,k,2)
    wloc1_3 = w_gpu(ll,j,k,3)
    wloc1_4 = w_gpu(ll,j,k,4)
    wloc1_5 = w_gpu(ll,j,k,5)
    l=2
    ll = i + 1 - l
    wloc2_1 = w_gpu(ll,j,k,1)
    wloc2_2 = w_gpu(ll,j,k,2)
    wloc2_3 = w_gpu(ll,j,k,3)
    wloc2_4 = w_gpu(ll,j,k,4)
    wloc2_5 = w_gpu(ll,j,k,5)
    l=3
    ll = i + 1 - l
    wloc3_1 = w_gpu(ll,j,k,1)
    wloc3_2 = w_gpu(ll,j,k,2)
    wloc3_3 = w_gpu(ll,j,k,3)
    wloc3_4 = w_gpu(ll,j,k,4)
    wloc3_5 = w_gpu(ll,j,k,5)
!
!   Derivatives of conservative variables
!   Inner derivatives
    dwdxi1 = -(cnr1*wloc1_1+cnr2*wloc2_1+cnr3*wloc3_1)
    dwdxi2 = -(cnr1*wloc1_2+cnr2*wloc2_2+cnr3*wloc3_2)
    dwdxi3 = -(cnr1*wloc1_3+cnr2*wloc2_3+cnr3*wloc3_3)
    dwdxi4 = -(cnr1*wloc1_4+cnr2*wloc2_4+cnr3*wloc3_4)
    dwdxi5 = -(cnr1*wloc1_5+cnr2*wloc2_5+cnr3*wloc3_5)
!   Outer derivatives
    dwdxo1 = wloc0_1-wloc1_1
    dwdxo2 = wloc0_2-wloc1_2
    dwdxo3 = wloc0_3-wloc1_3
    dwdxo4 = wloc0_4-wloc1_4
    dwdxo5 = wloc0_5-wloc1_5
!
!   Derivatives of characteristic variables
    !m=1
    dwcdxi1 = 0._mykind
    dwcdxo1 = 0._mykind
    !mm=1
    dwcdxi1 =dwcdxi1 +el_1_1*jacinv_1_1*dwdxi1 !  inner
    dwcdxi1 =dwcdxi1 +el_1_1*jacinv_1_2*dwdxi2 !  inner
    dwcdxi1 =dwcdxi1 +el_1_1*jacinv_1_3*dwdxi3 !  inner
    dwcdxi1 =dwcdxi1 +el_1_1*jacinv_1_4*dwdxi4 !  inner
    dwcdxi1 =dwcdxi1 +el_1_1*jacinv_1_5*dwdxi5 !  inner
    !mm=2
    dwcdxi1 =dwcdxi1 +el_1_2*jacinv_2_1*dwdxi1 !  inner
    dwcdxi1 =dwcdxi1 +el_1_2*jacinv_2_2*dwdxi2 !  inner
    dwcdxi1 =dwcdxi1 +el_1_2*jacinv_2_3*dwdxi3 !  inner
    dwcdxi1 =dwcdxi1 +el_1_2*jacinv_2_4*dwdxi4 !  inner
    dwcdxi1 =dwcdxi1 +el_1_2*jacinv_2_5*dwdxi5 !  inner
    !mm=3
    dwcdxi1 =dwcdxi1 +el_1_3*jacinv_3_1*dwdxi1 !  inner
    dwcdxi1 =dwcdxi1 +el_1_3*jacinv_3_2*dwdxi2 !  inner
    dwcdxi1 =dwcdxi1 +el_1_3*jacinv_3_3*dwdxi3 !  inner
    dwcdxi1 =dwcdxi1 +el_1_3*jacinv_3_4*dwdxi4 !  inner
    dwcdxi1 =dwcdxi1 +el_1_3*jacinv_3_5*dwdxi5 !  inner
    !mm=4
    dwcdxi1 =dwcdxi1 +el_1_4*jacinv_4_1*dwdxi1 !  inner
    dwcdxi1 =dwcdxi1 +el_1_4*jacinv_4_2*dwdxi2 !  inner
    dwcdxi1 =dwcdxi1 +el_1_4*jacinv_4_3*dwdxi3 !  inner
    dwcdxi1 =dwcdxi1 +el_1_4*jacinv_4_4*dwdxi4 !  inner
    dwcdxi1 =dwcdxi1 +el_1_4*jacinv_4_5*dwdxi5 !  inner
    !mm=5
    dwcdxi1 =dwcdxi1 +el_1_5*jacinv_5_1*dwdxi1 !  inner
    dwcdxo1 =dwcdxo1 +el_1_5*jacinv_5_1*dwdxo1 !  outer
    dwcdxi1 =dwcdxi1 +el_1_5*jacinv_5_2*dwdxi2 !  inner
    dwcdxo1 =dwcdxo1 +el_1_5*jacinv_5_2*dwdxo2 !  outer
    dwcdxi1 =dwcdxi1 +el_1_5*jacinv_5_3*dwdxi3 !  inner
    dwcdxo1 =dwcdxo1 +el_1_5*jacinv_5_3*dwdxo3 !  outer
    dwcdxi1 =dwcdxi1 +el_1_5*jacinv_5_4*dwdxi4 !  inner
    dwcdxo1 =dwcdxo1 +el_1_5*jacinv_5_4*dwdxo4 !  outer
    dwcdxi1 =dwcdxi1 +el_1_5*jacinv_5_5*dwdxi5 !  inner
    dwcdxo1 =dwcdxo1 +el_1_5*jacinv_5_5*dwdxo5 !  outer
    !m=2
    dwcdxi2 = 0._mykind
    dwcdxo2 = 0._mykind
    !mm=1
    dwcdxi2 =dwcdxi2 +el_2_1*jacinv_1_1*dwdxi1 !  inner
    dwcdxi2 =dwcdxi2 +el_2_1*jacinv_1_2*dwdxi2 !  inner
    dwcdxi2 =dwcdxi2 +el_2_1*jacinv_1_3*dwdxi3 !  inner
    dwcdxi2 =dwcdxi2 +el_2_1*jacinv_1_4*dwdxi4 !  inner
    dwcdxi2 =dwcdxi2 +el_2_1*jacinv_1_5*dwdxi5 !  inner
    !mm=2
    dwcdxi2 =dwcdxi2 +el_2_2*jacinv_2_1*dwdxi1 !  inner
    dwcdxi2 =dwcdxi2 +el_2_2*jacinv_2_2*dwdxi2 !  inner
    dwcdxi2 =dwcdxi2 +el_2_2*jacinv_2_3*dwdxi3 !  inner
    dwcdxi2 =dwcdxi2 +el_2_2*jacinv_2_4*dwdxi4 !  inner
    dwcdxi2 =dwcdxi2 +el_2_2*jacinv_2_5*dwdxi5 !  inner
    !mm=3
    dwcdxi2 =dwcdxi2 +el_2_3*jacinv_3_1*dwdxi1 !  inner
    dwcdxi2 =dwcdxi2 +el_2_3*jacinv_3_2*dwdxi2 !  inner
    dwcdxi2 =dwcdxi2 +el_2_3*jacinv_3_3*dwdxi3 !  inner
    dwcdxi2 =dwcdxi2 +el_2_3*jacinv_3_4*dwdxi4 !  inner
    dwcdxi2 =dwcdxi2 +el_2_3*jacinv_3_5*dwdxi5 !  inner
    !mm=4
    dwcdxi2 =dwcdxi2 +el_2_4*jacinv_4_1*dwdxi1 !  inner
    dwcdxi2 =dwcdxi2 +el_2_4*jacinv_4_2*dwdxi2 !  inner
    dwcdxi2 =dwcdxi2 +el_2_4*jacinv_4_3*dwdxi3 !  inner
    dwcdxi2 =dwcdxi2 +el_2_4*jacinv_4_4*dwdxi4 !  inner
    dwcdxi2 =dwcdxi2 +el_2_4*jacinv_4_5*dwdxi5 !  inner
    !mm=5
    dwcdxi2 =dwcdxi2 +el_2_5*jacinv_5_1*dwdxi1 !  inner
    dwcdxo2 =dwcdxo2 +el_2_5*jacinv_5_1*dwdxo1 !  outer
    dwcdxi2 =dwcdxi2 +el_2_5*jacinv_5_2*dwdxi2 !  inner
    dwcdxo2 =dwcdxo2 +el_2_5*jacinv_5_2*dwdxo2 !  outer
    dwcdxi2 =dwcdxi2 +el_2_5*jacinv_5_3*dwdxi3 !  inner
    dwcdxo2 =dwcdxo2 +el_2_5*jacinv_5_3*dwdxo3 !  outer
    dwcdxi2 =dwcdxi2 +el_2_5*jacinv_5_4*dwdxi4 !  inner
    dwcdxo2 =dwcdxo2 +el_2_5*jacinv_5_4*dwdxo4 !  outer
    dwcdxi2 =dwcdxi2 +el_2_5*jacinv_5_5*dwdxi5 !  inner
    dwcdxo2 =dwcdxo2 +el_2_5*jacinv_5_5*dwdxo5 !  outer
    !m=3
    dwcdxi3 = 0._mykind
    dwcdxo3 = 0._mykind
    !mm=1
    dwcdxi3 =dwcdxi3 +el_3_1*jacinv_1_1*dwdxi1 !  inner
    dwcdxi3 =dwcdxi3 +el_3_1*jacinv_1_2*dwdxi2 !  inner
    dwcdxi3 =dwcdxi3 +el_3_1*jacinv_1_3*dwdxi3 !  inner
    dwcdxi3 =dwcdxi3 +el_3_1*jacinv_1_4*dwdxi4 !  inner
    dwcdxi3 =dwcdxi3 +el_3_1*jacinv_1_5*dwdxi5 !  inner
    !mm=2
    dwcdxi3 =dwcdxi3 +el_3_2*jacinv_2_1*dwdxi1 !  inner
    dwcdxi3 =dwcdxi3 +el_3_2*jacinv_2_2*dwdxi2 !  inner
    dwcdxi3 =dwcdxi3 +el_3_2*jacinv_2_3*dwdxi3 !  inner
    dwcdxi3 =dwcdxi3 +el_3_2*jacinv_2_4*dwdxi4 !  inner
    dwcdxi3 =dwcdxi3 +el_3_2*jacinv_2_5*dwdxi5 !  inner
    !mm=3
    dwcdxi3 =dwcdxi3 +el_3_3*jacinv_3_1*dwdxi1 !  inner
    dwcdxi3 =dwcdxi3 +el_3_3*jacinv_3_2*dwdxi2 !  inner
    dwcdxi3 =dwcdxi3 +el_3_3*jacinv_3_3*dwdxi3 !  inner
    dwcdxi3 =dwcdxi3 +el_3_3*jacinv_3_4*dwdxi4 !  inner
    dwcdxi3 =dwcdxi3 +el_3_3*jacinv_3_5*dwdxi5 !  inner
    !mm=4
    dwcdxi3 =dwcdxi3 +el_3_4*jacinv_4_1*dwdxi1 !  inner
    dwcdxi3 =dwcdxi3 +el_3_4*jacinv_4_2*dwdxi2 !  inner
    dwcdxi3 =dwcdxi3 +el_3_4*jacinv_4_3*dwdxi3 !  inner
    dwcdxi3 =dwcdxi3 +el_3_4*jacinv_4_4*dwdxi4 !  inner
    dwcdxi3 =dwcdxi3 +el_3_4*jacinv_4_5*dwdxi5 !  inner
    !mm=5
    dwcdxi3 =dwcdxi3 +el_3_5*jacinv_5_1*dwdxi1 !  inner
    dwcdxo3 =dwcdxo3 +el_3_5*jacinv_5_1*dwdxo1 !  outer
    dwcdxi3 =dwcdxi3 +el_3_5*jacinv_5_2*dwdxi2 !  inner
    dwcdxo3 =dwcdxo3 +el_3_5*jacinv_5_2*dwdxo2 !  outer
    dwcdxi3 =dwcdxi3 +el_3_5*jacinv_5_3*dwdxi3 !  inner
    dwcdxo3 =dwcdxo3 +el_3_5*jacinv_5_3*dwdxo3 !  outer
    dwcdxi3 =dwcdxi3 +el_3_5*jacinv_5_4*dwdxi4 !  inner
    dwcdxo3 =dwcdxo3 +el_3_5*jacinv_5_4*dwdxo4 !  outer
    dwcdxi3 =dwcdxi3 +el_3_5*jacinv_5_5*dwdxi5 !  inner
    dwcdxo3 =dwcdxo3 +el_3_5*jacinv_5_5*dwdxo5 !  outer
    !m=4
    dwcdxi4 = 0._mykind
    dwcdxo4 = 0._mykind
    !mm=1
    dwcdxi4 =dwcdxi4 +el_4_1*jacinv_1_1*dwdxi1 !  inner
    dwcdxi4 =dwcdxi4 +el_4_1*jacinv_1_2*dwdxi2 !  inner
    dwcdxi4 =dwcdxi4 +el_4_1*jacinv_1_3*dwdxi3 !  inner
    dwcdxi4 =dwcdxi4 +el_4_1*jacinv_1_4*dwdxi4 !  inner
    dwcdxi4 =dwcdxi4 +el_4_1*jacinv_1_5*dwdxi5 !  inner
    !mm=2
    dwcdxi4 =dwcdxi4 +el_4_2*jacinv_2_1*dwdxi1 !  inner
    dwcdxi4 =dwcdxi4 +el_4_2*jacinv_2_2*dwdxi2 !  inner
    dwcdxi4 =dwcdxi4 +el_4_2*jacinv_2_3*dwdxi3 !  inner
    dwcdxi4 =dwcdxi4 +el_4_2*jacinv_2_4*dwdxi4 !  inner
    dwcdxi4 =dwcdxi4 +el_4_2*jacinv_2_5*dwdxi5 !  inner
    !mm=3
    dwcdxi4 =dwcdxi4 +el_4_3*jacinv_3_1*dwdxi1 !  inner
    dwcdxi4 =dwcdxi4 +el_4_3*jacinv_3_2*dwdxi2 !  inner
    dwcdxi4 =dwcdxi4 +el_4_3*jacinv_3_3*dwdxi3 !  inner
    dwcdxi4 =dwcdxi4 +el_4_3*jacinv_3_4*dwdxi4 !  inner
    dwcdxi4 =dwcdxi4 +el_4_3*jacinv_3_5*dwdxi5 !  inner
    !mm=4
    dwcdxi4 =dwcdxi4 +el_4_4*jacinv_4_1*dwdxi1 !  inner
    dwcdxi4 =dwcdxi4 +el_4_4*jacinv_4_2*dwdxi2 !  inner
    dwcdxi4 =dwcdxi4 +el_4_4*jacinv_4_3*dwdxi3 !  inner
    dwcdxi4 =dwcdxi4 +el_4_4*jacinv_4_4*dwdxi4 !  inner
    dwcdxi4 =dwcdxi4 +el_4_4*jacinv_4_5*dwdxi5 !  inner
    !mm=5
    dwcdxi4 =dwcdxi4 +el_4_5*jacinv_5_1*dwdxi1 !  inner
    dwcdxo4 =dwcdxo4 +el_4_5*jacinv_5_1*dwdxo1 !  outer
    dwcdxi4 =dwcdxi4 +el_4_5*jacinv_5_2*dwdxi2 !  inner
    dwcdxo4 =dwcdxo4 +el_4_5*jacinv_5_2*dwdxo2 !  outer
    dwcdxi4 =dwcdxi4 +el_4_5*jacinv_5_3*dwdxi3 !  inner
    dwcdxo4 =dwcdxo4 +el_4_5*jacinv_5_3*dwdxo3 !  outer
    dwcdxi4 =dwcdxi4 +el_4_5*jacinv_5_4*dwdxi4 !  inner
    dwcdxo4 =dwcdxo4 +el_4_5*jacinv_5_4*dwdxo4 !  outer
    dwcdxi4 =dwcdxi4 +el_4_5*jacinv_5_5*dwdxi5 !  inner
    dwcdxo4 =dwcdxo4 +el_4_5*jacinv_5_5*dwdxo5 !  outer
    !m=5
    dwcdxi5 = 0._mykind
    dwcdxo5 = 0._mykind
    !mm=1
    dwcdxi5 =dwcdxi5 +el_5_1*jacinv_1_1*dwdxi1 !  inner
    dwcdxi5 =dwcdxi5 +el_5_1*jacinv_1_2*dwdxi2 !  inner
    dwcdxi5 =dwcdxi5 +el_5_1*jacinv_1_3*dwdxi3 !  inner
    dwcdxi5 =dwcdxi5 +el_5_1*jacinv_1_4*dwdxi4 !  inner
    dwcdxi5 =dwcdxi5 +el_5_1*jacinv_1_5*dwdxi5 !  inner
    !mm=2
    dwcdxi5 =dwcdxi5 +el_5_2*jacinv_2_1*dwdxi1 !  inner
    dwcdxi5 =dwcdxi5 +el_5_2*jacinv_2_2*dwdxi2 !  inner
    dwcdxi5 =dwcdxi5 +el_5_2*jacinv_2_3*dwdxi3 !  inner
    dwcdxi5 =dwcdxi5 +el_5_2*jacinv_2_4*dwdxi4 !  inner
    dwcdxi5 =dwcdxi5 +el_5_2*jacinv_2_5*dwdxi5 !  inner
    !mm=3
    dwcdxi5 =dwcdxi5 +el_5_3*jacinv_3_1*dwdxi1 !  inner
    dwcdxi5 =dwcdxi5 +el_5_3*jacinv_3_2*dwdxi2 !  inner
    dwcdxi5 =dwcdxi5 +el_5_3*jacinv_3_3*dwdxi3 !  inner
    dwcdxi5 =dwcdxi5 +el_5_3*jacinv_3_4*dwdxi4 !  inner
    dwcdxi5 =dwcdxi5 +el_5_3*jacinv_3_5*dwdxi5 !  inner
    !mm=4
    dwcdxi5 =dwcdxi5 +el_5_4*jacinv_4_1*dwdxi1 !  inner
    dwcdxi5 =dwcdxi5 +el_5_4*jacinv_4_2*dwdxi2 !  inner
    dwcdxi5 =dwcdxi5 +el_5_4*jacinv_4_3*dwdxi3 !  inner
    dwcdxi5 =dwcdxi5 +el_5_4*jacinv_4_4*dwdxi4 !  inner
    dwcdxi5 =dwcdxi5 +el_5_4*jacinv_4_5*dwdxi5 !  inner
    !mm=5
    dwcdxi5 =dwcdxi5 +el_5_5*jacinv_5_1*dwdxi1 !  inner
    dwcdxo5 =dwcdxo5 +el_5_5*jacinv_5_1*dwdxo1 !  outer
    dwcdxi5 =dwcdxi5 +el_5_5*jacinv_5_2*dwdxi2 !  inner
    dwcdxo5 =dwcdxo5 +el_5_5*jacinv_5_2*dwdxo2 !  outer
    dwcdxi5 =dwcdxi5 +el_5_5*jacinv_5_3*dwdxi3 !  inner
    dwcdxo5 =dwcdxo5 +el_5_5*jacinv_5_3*dwdxo3 !  outer
    dwcdxi5 =dwcdxi5 +el_5_5*jacinv_5_4*dwdxi4 !  inner
    dwcdxo5 =dwcdxo5 +el_5_5*jacinv_5_4*dwdxo4 !  outer
    dwcdxi5 =dwcdxi5 +el_5_5*jacinv_5_5*dwdxi5 !  inner
    dwcdxo5 =dwcdxo5 +el_5_5*jacinv_5_5*dwdxo5 !  outer
!
!   Enforce LODI relations
    !m=1
    if (ev1>0._mykind) then ! Waves pointing out of the domain
     dwcdx1 = dwcdxi1
    else ! Waves entering the domain
     dwcdx1 = ibcnr_gpu(2)*dwcdxo1 ! n.r with or without relaxation
    endif
    !m=2
    if (ev2>0._mykind) then ! Waves pointing out of the domain
     dwcdx2 = dwcdxi2
    else ! Waves entering the domain
     dwcdx2 = ibcnr_gpu(2)*dwcdxo2 ! n.r. with or without relaxation
    endif
    !m=3
    if (ev3>0._mykind) then ! Waves pointing out of the domain
     dwcdx3 = dwcdxi3
    else ! Waves entering the domain
     dwcdx3 = ibcnr_gpu(2)*dwcdxo3 ! n.r with or without relaxation
    endif
    !m=4
    if (ev4>0._mykind) then ! Waves pointing out of the domain
     dwcdx4 = dwcdxi4
    else ! Waves entering the domain
     dwcdx4 = ibcnr_gpu(2)*dwcdxo4 ! n.r with or without relaxation
    endif
    !m=5
    if (ev5>0._mykind) then ! Waves pointing out of the domain
     dwcdx5 = dwcdxi5
    else ! Waves entering the domain
     dwcdx5 = ibcnr_gpu(2)*dwcdxo5 ! n.r with or without relaxation
    endif

!   Amplitude of characteristic waves
    dwcdx1 = dwcdx1 * ev1
    dwcdx2 = dwcdx2 * ev2
    dwcdx3 = dwcdx3 * ev3
    dwcdx4 = dwcdx4 * ev4
    dwcdx5 = dwcdx5 * ev5
!   Return to conservative variables 
    !m=1
    df = 0._mykind
    !mm=1
    df = df + jac_1_1*er_1_1*dwcdx1
    df = df + jac_1_1*er_1_2*dwcdx2
    df = df + jac_1_1*er_1_3*dwcdx3
    df = df + jac_1_1*er_1_4*dwcdx4
    df = df + jac_1_1*er_1_5*dwcdx5
    !mm=2
    df = df + jac_1_2*er_2_1*dwcdx1
    df = df + jac_1_2*er_2_2*dwcdx2
    df = df + jac_1_2*er_2_3*dwcdx3
    df = df + jac_1_2*er_2_4*dwcdx4
    df = df + jac_1_2*er_2_5*dwcdx5
    !mm=3
    df = df + jac_1_3*er_3_1*dwcdx1
    df = df + jac_1_3*er_3_2*dwcdx2
    df = df + jac_1_3*er_3_3*dwcdx3
    df = df + jac_1_3*er_3_4*dwcdx4
    df = df + jac_1_3*er_3_5*dwcdx5
    !mm=4
    df = df + jac_1_4*er_4_1*dwcdx1
    df = df + jac_1_4*er_4_2*dwcdx2
    df = df + jac_1_4*er_4_3*dwcdx3
    df = df + jac_1_4*er_4_4*dwcdx4
    df = df + jac_1_4*er_4_5*dwcdx5
    !mm=5
    df = df + jac_1_5*er_5_1*dwcdx1
    df = df + jac_1_5*er_5_2*dwcdx2
    df = df + jac_1_5*er_5_3*dwcdx3
    df = df + jac_1_5*er_5_4*dwcdx4
    df = df + jac_1_5*er_5_5*dwcdx5
    fl_gpu(i,j,k,1) = fl_gpu(i,j,k,1) + df * dcsidx_gpu(i)
    !m=2
    df = 0._mykind
    !mm=1
    df = df + jac_2_1*er_1_1*dwcdx1
    df = df + jac_2_1*er_1_2*dwcdx2
    df = df + jac_2_1*er_1_3*dwcdx3
    df = df + jac_2_1*er_1_4*dwcdx4
    df = df + jac_2_1*er_1_5*dwcdx5
    !mm=2
    df = df + jac_2_2*er_2_1*dwcdx1
    df = df + jac_2_2*er_2_2*dwcdx2
    df = df + jac_2_2*er_2_3*dwcdx3
    df = df + jac_2_2*er_2_4*dwcdx4
    df = df + jac_2_2*er_2_5*dwcdx5
    !mm=3
    df = df + jac_2_3*er_3_1*dwcdx1
    df = df + jac_2_3*er_3_2*dwcdx2
    df = df + jac_2_3*er_3_3*dwcdx3
    df = df + jac_2_3*er_3_4*dwcdx4
    df = df + jac_2_3*er_3_5*dwcdx5
    !mm=4
    df = df + jac_2_4*er_4_1*dwcdx1
    df = df + jac_2_4*er_4_2*dwcdx2
    df = df + jac_2_4*er_4_3*dwcdx3
    df = df + jac_2_4*er_4_4*dwcdx4
    df = df + jac_2_4*er_4_5*dwcdx5
    !mm=5
    df = df + jac_2_5*er_5_1*dwcdx1
    df = df + jac_2_5*er_5_2*dwcdx2
    df = df + jac_2_5*er_5_3*dwcdx3
    df = df + jac_2_5*er_5_4*dwcdx4
    df = df + jac_2_5*er_5_5*dwcdx5
    fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) + df * dcsidx_gpu(i)
    !m=3
    df = 0._mykind
    !mm=1
    df = df + jac_3_1*er_1_1*dwcdx1
    df = df + jac_3_1*er_1_2*dwcdx2
    df = df + jac_3_1*er_1_3*dwcdx3
    df = df + jac_3_1*er_1_4*dwcdx4
    df = df + jac_3_1*er_1_5*dwcdx5
    !mm=2
    df = df + jac_3_2*er_2_1*dwcdx1
    df = df + jac_3_2*er_2_2*dwcdx2
    df = df + jac_3_2*er_2_3*dwcdx3
    df = df + jac_3_2*er_2_4*dwcdx4
    df = df + jac_3_2*er_2_5*dwcdx5
    !mm=3
    df = df + jac_3_3*er_3_1*dwcdx1
    df = df + jac_3_3*er_3_2*dwcdx2
    df = df + jac_3_3*er_3_3*dwcdx3
    df = df + jac_3_3*er_3_4*dwcdx4
    df = df + jac_3_3*er_3_5*dwcdx5
    !mm=4
    df = df + jac_3_4*er_4_1*dwcdx1
    df = df + jac_3_4*er_4_2*dwcdx2
    df = df + jac_3_4*er_4_3*dwcdx3
    df = df + jac_3_4*er_4_4*dwcdx4
    df = df + jac_3_4*er_4_5*dwcdx5
    !mm=5
    df = df + jac_3_5*er_5_1*dwcdx1
    df = df + jac_3_5*er_5_2*dwcdx2
    df = df + jac_3_5*er_5_3*dwcdx3
    df = df + jac_3_5*er_5_4*dwcdx4
    df = df + jac_3_5*er_5_5*dwcdx5
    fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) + df * dcsidx_gpu(i)
    !m=4
    df = 0._mykind
    !mm=1
    df = df + jac_4_1*er_1_1*dwcdx1
    df = df + jac_4_1*er_1_2*dwcdx2
    df = df + jac_4_1*er_1_3*dwcdx3
    df = df + jac_4_1*er_1_4*dwcdx4
    df = df + jac_4_1*er_1_5*dwcdx5
    !mm=2
    df = df + jac_4_2*er_2_1*dwcdx1
    df = df + jac_4_2*er_2_2*dwcdx2
    df = df + jac_4_2*er_2_3*dwcdx3
    df = df + jac_4_2*er_2_4*dwcdx4
    df = df + jac_4_2*er_2_5*dwcdx5
    !mm=3
    df = df + jac_4_3*er_3_1*dwcdx1
    df = df + jac_4_3*er_3_2*dwcdx2
    df = df + jac_4_3*er_3_3*dwcdx3
    df = df + jac_4_3*er_3_4*dwcdx4
    df = df + jac_4_3*er_3_5*dwcdx5
    !mm=4
    df = df + jac_4_4*er_4_1*dwcdx1
    df = df + jac_4_4*er_4_2*dwcdx2
    df = df + jac_4_4*er_4_3*dwcdx3
    df = df + jac_4_4*er_4_4*dwcdx4
    df = df + jac_4_4*er_4_5*dwcdx5
    !mm=5
    df = df + jac_4_5*er_5_1*dwcdx1
    df = df + jac_4_5*er_5_2*dwcdx2
    df = df + jac_4_5*er_5_3*dwcdx3
    df = df + jac_4_5*er_5_4*dwcdx4
    df = df + jac_4_5*er_5_5*dwcdx5
    fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) + df * dcsidx_gpu(i)
    !m=5
    df = 0._mykind
    !mm=1
    df = df + jac_5_1*er_1_1*dwcdx1
    df = df + jac_5_1*er_1_2*dwcdx2
    df = df + jac_5_1*er_1_3*dwcdx3
    df = df + jac_5_1*er_1_4*dwcdx4
    df = df + jac_5_1*er_1_5*dwcdx5
    !mm=2
    df = df + jac_5_2*er_2_1*dwcdx1
    df = df + jac_5_2*er_2_2*dwcdx2
    df = df + jac_5_2*er_2_3*dwcdx3
    df = df + jac_5_2*er_2_4*dwcdx4
    df = df + jac_5_2*er_2_5*dwcdx5
    !mm=3
    df = df + jac_5_3*er_3_1*dwcdx1
    df = df + jac_5_3*er_3_2*dwcdx2
    df = df + jac_5_3*er_3_3*dwcdx3
    df = df + jac_5_3*er_3_4*dwcdx4
    df = df + jac_5_3*er_3_5*dwcdx5
    !mm=4
    df = df + jac_5_4*er_4_1*dwcdx1
    df = df + jac_5_4*er_4_2*dwcdx2
    df = df + jac_5_4*er_4_3*dwcdx3
    df = df + jac_5_4*er_4_4*dwcdx4
    df = df + jac_5_4*er_4_5*dwcdx5
    !mm=5
    df = df + jac_5_5*er_5_1*dwcdx1
    df = df + jac_5_5*er_5_2*dwcdx2
    df = df + jac_5_5*er_5_3*dwcdx3
    df = df + jac_5_5*er_5_4*dwcdx4
    df = df + jac_5_5*er_5_5*dwcdx5
    fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) + df * dcsidx_gpu(i)
!
   enddo  ! end of j-loop 
  enddo  ! end of k-loop
  !@cuf iercuda=cudaDeviceSynchronize()
!
 elseif (ilat==3) then  ! lower side
 elseif (ilat==4) then  ! upper side
!
! Redefining g(u)_y on the boundary
!
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
!
    j = ny
!
    rho =  w_gpu(i,j,k,1)
    ri  =  1._mykind/rho
    uu  =  w_gpu(i,j,k,2) * ri
    vv  =  w_gpu(i,j,k,3) * ri
    ww  =  w_gpu(i,j,k,4) * ri
    qq  =  0.5_mykind * (uu*uu + vv*vv + ww*ww)
!   pp  =  gm1 * (w_gpu(i,j,k,5) - w_gpu(i,j,k,1) * qq)
    pp  =  rho*temperature_gpu(i,j,k)
    cc  =  gamma * pp * ri
    c   =  sqrt(cc)
    ci  =  1._mykind/c 
!   Jacobian of conservative/primitive transformation (Eqn. (A.5) of Lodato et al, JCP 2008)
    jac_1_1  =  1._mykind
    jac_1_2  =  0._mykind
    jac_1_3  =  0._mykind
    jac_1_4  =  0._mykind
    jac_1_5  =  0._mykind
    jac_2_1  =  uu
    jac_2_2  =  rho
    jac_2_3  =  0._mykind
    jac_2_4  =  0._mykind
    jac_2_5  =  0._mykind
    jac_3_1  =  vv
    jac_3_2  =  0._mykind
    jac_3_3  =  rho
    jac_3_4  =  0._mykind
    jac_3_5  =  0._mykind
    jac_4_1  =  ww
    jac_4_2  =  0._mykind
    jac_4_3  =  0._mykind
    jac_4_4  =  rho
    jac_4_5  =  0._mykind
    jac_5_1  =  qq
    jac_5_2  =  w_gpu(i,j,k,2)
    jac_5_3  =  w_gpu(i,j,k,3)
    jac_5_4  =  w_gpu(i,j,k,4)
    jac_5_5  =  gm
!   Jacobian of inverse conservative/primitive transformation (Eqn. (A.6) of Lodato et al, JCP 2008)
    jacinv_1_1 =  1._mykind
    jacinv_1_2 =  0._mykind
    jacinv_1_3 =  0._mykind
    jacinv_1_4 =  0._mykind
    jacinv_1_5 =  0._mykind
    jacinv_2_1 = -uu*ri
    jacinv_2_2 =  ri
    jacinv_2_3 =  0._mykind
    jacinv_2_4 =  0._mykind
    jacinv_2_5 =  0._mykind
    jacinv_3_1 = -vv*ri
    jacinv_3_2 =  0._mykind
    jacinv_3_3 =  ri
    jacinv_3_4 =  0._mykind
    jacinv_3_5 =  0._mykind
    jacinv_4_1 = -ww*ri
    jacinv_4_2 =  0._mykind
    jacinv_4_3 =  0._mykind
    jacinv_4_4 =  ri
    jacinv_4_5 =  0._mykind
    jacinv_5_1 =  gm1*qq
    jacinv_5_2 = -gm1*uu
    jacinv_5_3 = -gm1*vv
    jacinv_5_4 = -gm1*ww
    jacinv_5_5 =  gm1
!   left eigenvectors matrix (Eqn._mykind (A.12_mykind) of Lodato et al, JCP 2008)
    el_1_1  =  0._mykind
    el_1_2  =  0._mykind
    el_1_3  = -rho*c
    el_1_4  =  0._mykind
    el_1_5  =  1._mykind
    el_2_1  =  0._mykind
    el_2_2  =  1._mykind
    el_2_3  =  0._mykind
    el_2_4  =  0._mykind
    el_2_5  =  0._mykind
    el_3_1  =  cc
    el_3_2  =  0._mykind
    el_3_3  =  0._mykind
    el_3_4  =  0._mykind
    el_3_5  = -1._mykind
    el_4_1  =  0._mykind
    el_4_2  =  0._mykind
    el_4_3  =  0._mykind
    el_4_4  =  1._mykind
    el_4_5  =  0._mykind
    el_5_1  =  0._mykind
    el_5_2  =  0._mykind
    el_5_3  =  rho*c
    el_5_4  =  0._mykind
    el_5_5  =  1._mykind
!   left eigenvectors matrix (Eqn._mykind (A.11_mykind) of Lodato et al, JCP 2008)
    er_1_1  =  0.5_mykind/cc
    er_1_2  =  0._mykind
    er_1_3  =  1._mykind/cc
    er_1_4  =  0._mykind
    er_1_5  =  0.5_mykind/cc
    er_2_1  =  0._mykind
    er_2_2  =  1._mykind
    er_2_3  =  0._mykind
    er_2_4  =  0._mykind
    er_2_5  =  0._mykind
    er_3_1  = -0.5_mykind*ri*ci
    er_3_2  =  0._mykind
    er_3_3  =  0._mykind
    er_3_4  =  0._mykind
    er_3_5  =  0.5_mykind*ri*ci
    er_4_1  =  0._mykind
    er_4_2  =  0._mykind
    er_4_3  =  0._mykind 
    er_4_4  =  1._mykind
    er_4_5  =  0._mykind
    er_5_1  =  0.5_mykind
    er_5_2  =  0._mykind
    er_5_3  =  0._mykind
    er_5_4  =  0._mykind
    er_5_5  =  0.5_mykind
!
!   Eigenvalues 
    ev1     =  vv-c
    ev2     =  vv
    ev3     =  vv
    ev4     =  vv
    ev5     =  vv+c
!
!   Local array of conservative variables
    wloc0_1 = winf_gpu(1)
    wloc0_2 = winf_gpu(2)
    wloc0_3 = winf_gpu(3)
    wloc0_4 = winf_gpu(4)
    wloc0_5 = winf_gpu(5)
    l=1
    ll = j + 1 - l
    wloc1_1 = w_gpu(i,ll,k,1)
    wloc1_2 = w_gpu(i,ll,k,2)
    wloc1_3 = w_gpu(i,ll,k,3)
    wloc1_4 = w_gpu(i,ll,k,4)
    wloc1_5 = w_gpu(i,ll,k,5)
    l=2
    ll = j + 1 - l
    wloc2_1 = w_gpu(i,ll,k,1)
    wloc2_2 = w_gpu(i,ll,k,2)
    wloc2_3 = w_gpu(i,ll,k,3)
    wloc2_4 = w_gpu(i,ll,k,4)
    wloc2_5 = w_gpu(i,ll,k,5)
    l=3
    ll = j + 1 - l
    wloc3_1 = w_gpu(i,ll,k,1)
    wloc3_2 = w_gpu(i,ll,k,2)
    wloc3_3 = w_gpu(i,ll,k,3)
    wloc3_4 = w_gpu(i,ll,k,4)
    wloc3_5 = w_gpu(i,ll,k,5)
!
!   Derivatives of conservative variables
!   Inner derivatives
    dwdyi1 = -(cnr1*wloc1_1+cnr2*wloc2_1+cnr3*wloc3_1)
    dwdyi2 = -(cnr1*wloc1_2+cnr2*wloc2_2+cnr3*wloc3_2)
    dwdyi3 = -(cnr1*wloc1_3+cnr2*wloc2_3+cnr3*wloc3_3)
    dwdyi4 = -(cnr1*wloc1_4+cnr2*wloc2_4+cnr3*wloc3_4)
    dwdyi5 = -(cnr1*wloc1_5+cnr2*wloc2_5+cnr3*wloc3_5)
!   Outer derivatives
    dwdyo1 = wloc0_1-wloc1_1
    dwdyo2 = wloc0_2-wloc1_2
    dwdyo3 = wloc0_3-wloc1_3
    dwdyo4 = wloc0_4-wloc1_4
    dwdyo5 = wloc0_5-wloc1_5
!
!   Derivatives of characteristic variables
    !m=1
    dwcdyi1 = 0._mykind
    dwcdyo1 = 0._mykind
    !mm=1
    dwcdyi1 = dwcdyi1+el_1_1*jacinv_1_1*dwdyi1 !  inner
    dwcdyo1 = dwcdyo1+el_1_1*jacinv_1_1*dwdyo1 !  outer
    dwcdyi1 = dwcdyi1+el_1_1*jacinv_1_2*dwdyi2 !  inner
    dwcdyo1 = dwcdyo1+el_1_1*jacinv_1_2*dwdyo2 !  outer
    dwcdyi1 = dwcdyi1+el_1_1*jacinv_1_3*dwdyi3 !  inner
    dwcdyo1 = dwcdyo1+el_1_1*jacinv_1_3*dwdyo3 !  outer
    dwcdyi1 = dwcdyi1+el_1_1*jacinv_1_4*dwdyi4 !  inner
    dwcdyo1 = dwcdyo1+el_1_1*jacinv_1_4*dwdyo4 !  outer
    dwcdyi1 = dwcdyi1+el_1_1*jacinv_1_5*dwdyi5 !  inner
    dwcdyo1 = dwcdyo1+el_1_1*jacinv_1_5*dwdyo5 !  outer
    !mm=2
    dwcdyi1 = dwcdyi1+el_1_2*jacinv_2_1*dwdyi1 !  inner
    dwcdyo1 = dwcdyo1+el_1_2*jacinv_2_1*dwdyo1 !  outer
    dwcdyi1 = dwcdyi1+el_1_2*jacinv_2_2*dwdyi2 !  inner
    dwcdyo1 = dwcdyo1+el_1_2*jacinv_2_2*dwdyo2 !  outer
    dwcdyi1 = dwcdyi1+el_1_2*jacinv_2_3*dwdyi3 !  inner
    dwcdyo1 = dwcdyo1+el_1_2*jacinv_2_3*dwdyo3 !  outer
    dwcdyi1 = dwcdyi1+el_1_2*jacinv_2_4*dwdyi4 !  inner
    dwcdyo1 = dwcdyo1+el_1_2*jacinv_2_4*dwdyo4 !  outer
    dwcdyi1 = dwcdyi1+el_1_2*jacinv_2_5*dwdyi5 !  inner
    dwcdyo1 = dwcdyo1+el_1_2*jacinv_2_5*dwdyo5 !  outer
    !mm=3
    dwcdyi1 = dwcdyi1+el_1_3*jacinv_3_1*dwdyi1 !  inner
    dwcdyo1 = dwcdyo1+el_1_3*jacinv_3_1*dwdyo1 !  outer
    dwcdyi1 = dwcdyi1+el_1_3*jacinv_3_2*dwdyi2 !  inner
    dwcdyo1 = dwcdyo1+el_1_3*jacinv_3_2*dwdyo2 !  outer
    dwcdyi1 = dwcdyi1+el_1_3*jacinv_3_3*dwdyi3 !  inner
    dwcdyo1 = dwcdyo1+el_1_3*jacinv_3_3*dwdyo3 !  outer
    dwcdyi1 = dwcdyi1+el_1_3*jacinv_3_4*dwdyi4 !  inner
    dwcdyo1 = dwcdyo1+el_1_3*jacinv_3_4*dwdyo4 !  outer
    dwcdyi1 = dwcdyi1+el_1_3*jacinv_3_5*dwdyi5 !  inner
    dwcdyo1 = dwcdyo1+el_1_3*jacinv_3_5*dwdyo5 !  outer
    !mm=4
    dwcdyi1 = dwcdyi1+el_1_4*jacinv_4_1*dwdyi1 !  inner
    dwcdyo1 = dwcdyo1+el_1_4*jacinv_4_1*dwdyo1 !  outer
    dwcdyi1 = dwcdyi1+el_1_4*jacinv_4_2*dwdyi2 !  inner
    dwcdyo1 = dwcdyo1+el_1_4*jacinv_4_2*dwdyo2 !  outer
    dwcdyi1 = dwcdyi1+el_1_4*jacinv_4_3*dwdyi3 !  inner
    dwcdyo1 = dwcdyo1+el_1_4*jacinv_4_3*dwdyo3 !  outer
    dwcdyi1 = dwcdyi1+el_1_4*jacinv_4_4*dwdyi4 !  inner
    dwcdyo1 = dwcdyo1+el_1_4*jacinv_4_4*dwdyo4 !  outer
    dwcdyi1 = dwcdyi1+el_1_4*jacinv_4_5*dwdyi5 !  inner
    dwcdyo1 = dwcdyo1+el_1_4*jacinv_4_5*dwdyo5 !  outer
    !mm=5
    dwcdyi1 = dwcdyi1+el_1_5*jacinv_5_1*dwdyi1 !  inner
    dwcdyo1 = dwcdyo1+el_1_5*jacinv_5_1*dwdyo1 !  outer
    dwcdyi1 = dwcdyi1+el_1_5*jacinv_5_2*dwdyi2 !  inner
    dwcdyo1 = dwcdyo1+el_1_5*jacinv_5_2*dwdyo2 !  outer
    dwcdyi1 = dwcdyi1+el_1_5*jacinv_5_3*dwdyi3 !  inner
    dwcdyo1 = dwcdyo1+el_1_5*jacinv_5_3*dwdyo3 !  outer
    dwcdyi1 = dwcdyi1+el_1_5*jacinv_5_4*dwdyi4 !  inner
    dwcdyo1 = dwcdyo1+el_1_5*jacinv_5_4*dwdyo4 !  outer
    dwcdyi1 = dwcdyi1+el_1_5*jacinv_5_5*dwdyi5 !  inner
    dwcdyo1 = dwcdyo1+el_1_5*jacinv_5_5*dwdyo5 !  outer
    !m=2
    dwcdyi2 = 0._mykind
    dwcdyo2 = 0._mykind
    !mm=1
    dwcdyi2 = dwcdyi2+el_2_1*jacinv_1_1*dwdyi1 !  inner
    dwcdyo2 = dwcdyo2+el_2_1*jacinv_1_1*dwdyo1 !  outer
    dwcdyi2 = dwcdyi2+el_2_1*jacinv_1_2*dwdyi2 !  inner
    dwcdyo2 = dwcdyo2+el_2_1*jacinv_1_2*dwdyo2 !  outer
    dwcdyi2 = dwcdyi2+el_2_1*jacinv_1_3*dwdyi3 !  inner
    dwcdyo2 = dwcdyo2+el_2_1*jacinv_1_3*dwdyo3 !  outer
    dwcdyi2 = dwcdyi2+el_2_1*jacinv_1_4*dwdyi4 !  inner
    dwcdyo2 = dwcdyo2+el_2_1*jacinv_1_4*dwdyo4 !  outer
    dwcdyi2 = dwcdyi2+el_2_1*jacinv_1_5*dwdyi5 !  inner
    dwcdyo2 = dwcdyo2+el_2_1*jacinv_1_5*dwdyo5 !  outer
    !mm=2
    dwcdyi2 = dwcdyi2+el_2_2*jacinv_2_1*dwdyi1 !  inner
    dwcdyo2 = dwcdyo2+el_2_2*jacinv_2_1*dwdyo1 !  outer
    dwcdyi2 = dwcdyi2+el_2_2*jacinv_2_2*dwdyi2 !  inner
    dwcdyo2 = dwcdyo2+el_2_2*jacinv_2_2*dwdyo2 !  outer
    dwcdyi2 = dwcdyi2+el_2_2*jacinv_2_3*dwdyi3 !  inner
    dwcdyo2 = dwcdyo2+el_2_2*jacinv_2_3*dwdyo3 !  outer
    dwcdyi2 = dwcdyi2+el_2_2*jacinv_2_4*dwdyi4 !  inner
    dwcdyo2 = dwcdyo2+el_2_2*jacinv_2_4*dwdyo4 !  outer
    dwcdyi2 = dwcdyi2+el_2_2*jacinv_2_5*dwdyi5 !  inner
    dwcdyo2 = dwcdyo2+el_2_2*jacinv_2_5*dwdyo5 !  outer
    !mm=3
    dwcdyi2 = dwcdyi2+el_2_3*jacinv_3_1*dwdyi1 !  inner
    dwcdyo2 = dwcdyo2+el_2_3*jacinv_3_1*dwdyo1 !  outer
    dwcdyi2 = dwcdyi2+el_2_3*jacinv_3_2*dwdyi2 !  inner
    dwcdyo2 = dwcdyo2+el_2_3*jacinv_3_2*dwdyo2 !  outer
    dwcdyi2 = dwcdyi2+el_2_3*jacinv_3_3*dwdyi3 !  inner
    dwcdyo2 = dwcdyo2+el_2_3*jacinv_3_3*dwdyo3 !  outer
    dwcdyi2 = dwcdyi2+el_2_3*jacinv_3_4*dwdyi4 !  inner
    dwcdyo2 = dwcdyo2+el_2_3*jacinv_3_4*dwdyo4 !  outer
    dwcdyi2 = dwcdyi2+el_2_3*jacinv_3_5*dwdyi5 !  inner
    dwcdyo2 = dwcdyo2+el_2_3*jacinv_3_5*dwdyo5 !  outer
    !mm=4
    dwcdyi2 = dwcdyi2+el_2_4*jacinv_4_1*dwdyi1 !  inner
    dwcdyo2 = dwcdyo2+el_2_4*jacinv_4_1*dwdyo1 !  outer
    dwcdyi2 = dwcdyi2+el_2_4*jacinv_4_2*dwdyi2 !  inner
    dwcdyo2 = dwcdyo2+el_2_4*jacinv_4_2*dwdyo2 !  outer
    dwcdyi2 = dwcdyi2+el_2_4*jacinv_4_3*dwdyi3 !  inner
    dwcdyo2 = dwcdyo2+el_2_4*jacinv_4_3*dwdyo3 !  outer
    dwcdyi2 = dwcdyi2+el_2_4*jacinv_4_4*dwdyi4 !  inner
    dwcdyo2 = dwcdyo2+el_2_4*jacinv_4_4*dwdyo4 !  outer
    dwcdyi2 = dwcdyi2+el_2_4*jacinv_4_5*dwdyi5 !  inner
    dwcdyo2 = dwcdyo2+el_2_4*jacinv_4_5*dwdyo5 !  outer
    !mm=5
    dwcdyi2 = dwcdyi2+el_2_5*jacinv_5_1*dwdyi1 !  inner
    dwcdyo2 = dwcdyo2+el_2_5*jacinv_5_1*dwdyo1 !  outer
    dwcdyi2 = dwcdyi2+el_2_5*jacinv_5_2*dwdyi2 !  inner
    dwcdyo2 = dwcdyo2+el_2_5*jacinv_5_2*dwdyo2 !  outer
    dwcdyi2 = dwcdyi2+el_2_5*jacinv_5_3*dwdyi3 !  inner
    dwcdyo2 = dwcdyo2+el_2_5*jacinv_5_3*dwdyo3 !  outer
    dwcdyi2 = dwcdyi2+el_2_5*jacinv_5_4*dwdyi4 !  inner
    dwcdyo2 = dwcdyo2+el_2_5*jacinv_5_4*dwdyo4 !  outer
    dwcdyi2 = dwcdyi2+el_2_5*jacinv_5_5*dwdyi5 !  inner
    dwcdyo2 = dwcdyo2+el_2_5*jacinv_5_5*dwdyo5 !  outer
    !m=3
    dwcdyi3 = 0._mykind
    dwcdyo3 = 0._mykind
    !mm=1
    dwcdyi3 = dwcdyi3+el_3_1*jacinv_1_1*dwdyi1 !  inner
    dwcdyo3 = dwcdyo3+el_3_1*jacinv_1_1*dwdyo1 !  outer
    dwcdyi3 = dwcdyi3+el_3_1*jacinv_1_2*dwdyi2 !  inner
    dwcdyo3 = dwcdyo3+el_3_1*jacinv_1_2*dwdyo2 !  outer
    dwcdyi3 = dwcdyi3+el_3_1*jacinv_1_3*dwdyi3 !  inner
    dwcdyo3 = dwcdyo3+el_3_1*jacinv_1_3*dwdyo3 !  outer
    dwcdyi3 = dwcdyi3+el_3_1*jacinv_1_4*dwdyi4 !  inner
    dwcdyo3 = dwcdyo3+el_3_1*jacinv_1_4*dwdyo4 !  outer
    dwcdyi3 = dwcdyi3+el_3_1*jacinv_1_5*dwdyi5 !  inner
    dwcdyo3 = dwcdyo3+el_3_1*jacinv_1_5*dwdyo5 !  outer
    !mm=2
    dwcdyi3 = dwcdyi3+el_3_2*jacinv_2_1*dwdyi1 !  inner
    dwcdyo3 = dwcdyo3+el_3_2*jacinv_2_1*dwdyo1 !  outer
    dwcdyi3 = dwcdyi3+el_3_2*jacinv_2_2*dwdyi2 !  inner
    dwcdyo3 = dwcdyo3+el_3_2*jacinv_2_2*dwdyo2 !  outer
    dwcdyi3 = dwcdyi3+el_3_2*jacinv_2_3*dwdyi3 !  inner
    dwcdyo3 = dwcdyo3+el_3_2*jacinv_2_3*dwdyo3 !  outer
    dwcdyi3 = dwcdyi3+el_3_2*jacinv_2_4*dwdyi4 !  inner
    dwcdyo3 = dwcdyo3+el_3_2*jacinv_2_4*dwdyo4 !  outer
    dwcdyi3 = dwcdyi3+el_3_2*jacinv_2_5*dwdyi5 !  inner
    dwcdyo3 = dwcdyo3+el_3_2*jacinv_2_5*dwdyo5 !  outer
    !mm=3
    dwcdyi3 = dwcdyi3+el_3_3*jacinv_3_1*dwdyi1 !  inner
    dwcdyo3 = dwcdyo3+el_3_3*jacinv_3_1*dwdyo1 !  outer
    dwcdyi3 = dwcdyi3+el_3_3*jacinv_3_2*dwdyi2 !  inner
    dwcdyo3 = dwcdyo3+el_3_3*jacinv_3_2*dwdyo2 !  outer
    dwcdyi3 = dwcdyi3+el_3_3*jacinv_3_3*dwdyi3 !  inner
    dwcdyo3 = dwcdyo3+el_3_3*jacinv_3_3*dwdyo3 !  outer
    dwcdyi3 = dwcdyi3+el_3_3*jacinv_3_4*dwdyi4 !  inner
    dwcdyo3 = dwcdyo3+el_3_3*jacinv_3_4*dwdyo4 !  outer
    dwcdyi3 = dwcdyi3+el_3_3*jacinv_3_5*dwdyi5 !  inner
    dwcdyo3 = dwcdyo3+el_3_3*jacinv_3_5*dwdyo5 !  outer
    !mm=4
    dwcdyi3 = dwcdyi3+el_3_4*jacinv_4_1*dwdyi1 !  inner
    dwcdyo3 = dwcdyo3+el_3_4*jacinv_4_1*dwdyo1 !  outer
    dwcdyi3 = dwcdyi3+el_3_4*jacinv_4_2*dwdyi2 !  inner
    dwcdyo3 = dwcdyo3+el_3_4*jacinv_4_2*dwdyo2 !  outer
    dwcdyi3 = dwcdyi3+el_3_4*jacinv_4_3*dwdyi3 !  inner
    dwcdyo3 = dwcdyo3+el_3_4*jacinv_4_3*dwdyo3 !  outer
    dwcdyi3 = dwcdyi3+el_3_4*jacinv_4_4*dwdyi4 !  inner
    dwcdyo3 = dwcdyo3+el_3_4*jacinv_4_4*dwdyo4 !  outer
    dwcdyi3 = dwcdyi3+el_3_4*jacinv_4_5*dwdyi5 !  inner
    dwcdyo3 = dwcdyo3+el_3_4*jacinv_4_5*dwdyo5 !  outer
    !mm=5
    dwcdyi3 = dwcdyi3+el_3_5*jacinv_5_1*dwdyi1 !  inner
    dwcdyo3 = dwcdyo3+el_3_5*jacinv_5_1*dwdyo1 !  outer
    dwcdyi3 = dwcdyi3+el_3_5*jacinv_5_2*dwdyi2 !  inner
    dwcdyo3 = dwcdyo3+el_3_5*jacinv_5_2*dwdyo2 !  outer
    dwcdyi3 = dwcdyi3+el_3_5*jacinv_5_3*dwdyi3 !  inner
    dwcdyo3 = dwcdyo3+el_3_5*jacinv_5_3*dwdyo3 !  outer
    dwcdyi3 = dwcdyi3+el_3_5*jacinv_5_4*dwdyi4 !  inner
    dwcdyo3 = dwcdyo3+el_3_5*jacinv_5_4*dwdyo4 !  outer
    dwcdyi3 = dwcdyi3+el_3_5*jacinv_5_5*dwdyi5 !  inner
    dwcdyo3 = dwcdyo3+el_3_5*jacinv_5_5*dwdyo5 !  outer
    !m=4
    dwcdyi4 = 0._mykind
    dwcdyo4 = 0._mykind
    !mm=1
    dwcdyi4 = dwcdyi4+el_4_1*jacinv_1_1*dwdyi1 !  inner
    dwcdyo4 = dwcdyo4+el_4_1*jacinv_1_1*dwdyo1 !  outer
    dwcdyi4 = dwcdyi4+el_4_1*jacinv_1_2*dwdyi2 !  inner
    dwcdyo4 = dwcdyo4+el_4_1*jacinv_1_2*dwdyo2 !  outer
    dwcdyi4 = dwcdyi4+el_4_1*jacinv_1_3*dwdyi3 !  inner
    dwcdyo4 = dwcdyo4+el_4_1*jacinv_1_3*dwdyo3 !  outer
    dwcdyi4 = dwcdyi4+el_4_1*jacinv_1_4*dwdyi4 !  inner
    dwcdyo4 = dwcdyo4+el_4_1*jacinv_1_4*dwdyo4 !  outer
    dwcdyi4 = dwcdyi4+el_4_1*jacinv_1_5*dwdyi5 !  inner
    dwcdyo4 = dwcdyo4+el_4_1*jacinv_1_5*dwdyo5 !  outer
    !mm=2
    dwcdyi4 = dwcdyi4+el_4_2*jacinv_2_1*dwdyi1 !  inner
    dwcdyo4 = dwcdyo4+el_4_2*jacinv_2_1*dwdyo1 !  outer
    dwcdyi4 = dwcdyi4+el_4_2*jacinv_2_2*dwdyi2 !  inner
    dwcdyo4 = dwcdyo4+el_4_2*jacinv_2_2*dwdyo2 !  outer
    dwcdyi4 = dwcdyi4+el_4_2*jacinv_2_3*dwdyi3 !  inner
    dwcdyo4 = dwcdyo4+el_4_2*jacinv_2_3*dwdyo3 !  outer
    dwcdyi4 = dwcdyi4+el_4_2*jacinv_2_4*dwdyi4 !  inner
    dwcdyo4 = dwcdyo4+el_4_2*jacinv_2_4*dwdyo4 !  outer
    dwcdyi4 = dwcdyi4+el_4_2*jacinv_2_5*dwdyi5 !  inner
    dwcdyo4 = dwcdyo4+el_4_2*jacinv_2_5*dwdyo5 !  outer
    !mm=3
    dwcdyi4 = dwcdyi4+el_4_3*jacinv_3_1*dwdyi1 !  inner
    dwcdyo4 = dwcdyo4+el_4_3*jacinv_3_1*dwdyo1 !  outer
    dwcdyi4 = dwcdyi4+el_4_3*jacinv_3_2*dwdyi2 !  inner
    dwcdyo4 = dwcdyo4+el_4_3*jacinv_3_2*dwdyo2 !  outer
    dwcdyi4 = dwcdyi4+el_4_3*jacinv_3_3*dwdyi3 !  inner
    dwcdyo4 = dwcdyo4+el_4_3*jacinv_3_3*dwdyo3 !  outer
    dwcdyi4 = dwcdyi4+el_4_3*jacinv_3_4*dwdyi4 !  inner
    dwcdyo4 = dwcdyo4+el_4_3*jacinv_3_4*dwdyo4 !  outer
    dwcdyi4 = dwcdyi4+el_4_3*jacinv_3_5*dwdyi5 !  inner
    dwcdyo4 = dwcdyo4+el_4_3*jacinv_3_5*dwdyo5 !  outer
    !mm=4
    dwcdyi4 = dwcdyi4+el_4_4*jacinv_4_1*dwdyi1 !  inner
    dwcdyo4 = dwcdyo4+el_4_4*jacinv_4_1*dwdyo1 !  outer
    dwcdyi4 = dwcdyi4+el_4_4*jacinv_4_2*dwdyi2 !  inner
    dwcdyo4 = dwcdyo4+el_4_4*jacinv_4_2*dwdyo2 !  outer
    dwcdyi4 = dwcdyi4+el_4_4*jacinv_4_3*dwdyi3 !  inner
    dwcdyo4 = dwcdyo4+el_4_4*jacinv_4_3*dwdyo3 !  outer
    dwcdyi4 = dwcdyi4+el_4_4*jacinv_4_4*dwdyi4 !  inner
    dwcdyo4 = dwcdyo4+el_4_4*jacinv_4_4*dwdyo4 !  outer
    dwcdyi4 = dwcdyi4+el_4_4*jacinv_4_5*dwdyi5 !  inner
    dwcdyo4 = dwcdyo4+el_4_4*jacinv_4_5*dwdyo5 !  outer
    !mm=5
    dwcdyi4 = dwcdyi4+el_4_5*jacinv_5_1*dwdyi1 !  inner
    dwcdyo4 = dwcdyo4+el_4_5*jacinv_5_1*dwdyo1 !  outer
    dwcdyi4 = dwcdyi4+el_4_5*jacinv_5_2*dwdyi2 !  inner
    dwcdyo4 = dwcdyo4+el_4_5*jacinv_5_2*dwdyo2 !  outer
    dwcdyi4 = dwcdyi4+el_4_5*jacinv_5_3*dwdyi3 !  inner
    dwcdyo4 = dwcdyo4+el_4_5*jacinv_5_3*dwdyo3 !  outer
    dwcdyi4 = dwcdyi4+el_4_5*jacinv_5_4*dwdyi4 !  inner
    dwcdyo4 = dwcdyo4+el_4_5*jacinv_5_4*dwdyo4 !  outer
    dwcdyi4 = dwcdyi4+el_4_5*jacinv_5_5*dwdyi5 !  inner
    dwcdyo4 = dwcdyo4+el_4_5*jacinv_5_5*dwdyo5 !  outer
    !m=5
    dwcdyi5 = 0._mykind
    dwcdyo5 = 0._mykind
    !mm=1
    dwcdyi5 = dwcdyi5+el_5_1*jacinv_1_1*dwdyi1 !  inner
    dwcdyo5 = dwcdyo5+el_5_1*jacinv_1_1*dwdyo1 !  outer
    dwcdyi5 = dwcdyi5+el_5_1*jacinv_1_2*dwdyi2 !  inner
    dwcdyo5 = dwcdyo5+el_5_1*jacinv_1_2*dwdyo2 !  outer
    dwcdyi5 = dwcdyi5+el_5_1*jacinv_1_3*dwdyi3 !  inner
    dwcdyo5 = dwcdyo5+el_5_1*jacinv_1_3*dwdyo3 !  outer
    dwcdyi5 = dwcdyi5+el_5_1*jacinv_1_4*dwdyi4 !  inner
    dwcdyo5 = dwcdyo5+el_5_1*jacinv_1_4*dwdyo4 !  outer
    dwcdyi5 = dwcdyi5+el_5_1*jacinv_1_5*dwdyi5 !  inner
    dwcdyo5 = dwcdyo5+el_5_1*jacinv_1_5*dwdyo5 !  outer
    !mm=2
    dwcdyi5 = dwcdyi5+el_5_2*jacinv_2_1*dwdyi1 !  inner
    dwcdyo5 = dwcdyo5+el_5_2*jacinv_2_1*dwdyo1 !  outer
    dwcdyi5 = dwcdyi5+el_5_2*jacinv_2_2*dwdyi2 !  inner
    dwcdyo5 = dwcdyo5+el_5_2*jacinv_2_2*dwdyo2 !  outer
    dwcdyi5 = dwcdyi5+el_5_2*jacinv_2_3*dwdyi3 !  inner
    dwcdyo5 = dwcdyo5+el_5_2*jacinv_2_3*dwdyo3 !  outer
    dwcdyi5 = dwcdyi5+el_5_2*jacinv_2_4*dwdyi4 !  inner
    dwcdyo5 = dwcdyo5+el_5_2*jacinv_2_4*dwdyo4 !  outer
    dwcdyi5 = dwcdyi5+el_5_2*jacinv_2_5*dwdyi5 !  inner
    dwcdyo5 = dwcdyo5+el_5_2*jacinv_2_5*dwdyo5 !  outer
    !mm=3
    dwcdyi5 = dwcdyi5+el_5_3*jacinv_3_1*dwdyi1 !  inner
    dwcdyo5 = dwcdyo5+el_5_3*jacinv_3_1*dwdyo1 !  outer
    dwcdyi5 = dwcdyi5+el_5_3*jacinv_3_2*dwdyi2 !  inner
    dwcdyo5 = dwcdyo5+el_5_3*jacinv_3_2*dwdyo2 !  outer
    dwcdyi5 = dwcdyi5+el_5_3*jacinv_3_3*dwdyi3 !  inner
    dwcdyo5 = dwcdyo5+el_5_3*jacinv_3_3*dwdyo3 !  outer
    dwcdyi5 = dwcdyi5+el_5_3*jacinv_3_4*dwdyi4 !  inner
    dwcdyo5 = dwcdyo5+el_5_3*jacinv_3_4*dwdyo4 !  outer
    dwcdyi5 = dwcdyi5+el_5_3*jacinv_3_5*dwdyi5 !  inner
    dwcdyo5 = dwcdyo5+el_5_3*jacinv_3_5*dwdyo5 !  outer
    !mm=4
    dwcdyi5 = dwcdyi5+el_5_4*jacinv_4_1*dwdyi1 !  inner
    dwcdyo5 = dwcdyo5+el_5_4*jacinv_4_1*dwdyo1 !  outer
    dwcdyi5 = dwcdyi5+el_5_4*jacinv_4_2*dwdyi2 !  inner
    dwcdyo5 = dwcdyo5+el_5_4*jacinv_4_2*dwdyo2 !  outer
    dwcdyi5 = dwcdyi5+el_5_4*jacinv_4_3*dwdyi3 !  inner
    dwcdyo5 = dwcdyo5+el_5_4*jacinv_4_3*dwdyo3 !  outer
    dwcdyi5 = dwcdyi5+el_5_4*jacinv_4_4*dwdyi4 !  inner
    dwcdyo5 = dwcdyo5+el_5_4*jacinv_4_4*dwdyo4 !  outer
    dwcdyi5 = dwcdyi5+el_5_4*jacinv_4_5*dwdyi5 !  inner
    dwcdyo5 = dwcdyo5+el_5_4*jacinv_4_5*dwdyo5 !  outer
    !mm=5
    dwcdyi5 = dwcdyi5+el_5_5*jacinv_5_1*dwdyi1 !  inner
    dwcdyo5 = dwcdyo5+el_5_5*jacinv_5_1*dwdyo1 !  outer
    dwcdyi5 = dwcdyi5+el_5_5*jacinv_5_2*dwdyi2 !  inner
    dwcdyo5 = dwcdyo5+el_5_5*jacinv_5_2*dwdyo2 !  outer
    dwcdyi5 = dwcdyi5+el_5_5*jacinv_5_3*dwdyi3 !  inner
    dwcdyo5 = dwcdyo5+el_5_5*jacinv_5_3*dwdyo3 !  outer
    dwcdyi5 = dwcdyi5+el_5_5*jacinv_5_4*dwdyi4 !  inner
    dwcdyo5 = dwcdyo5+el_5_5*jacinv_5_4*dwdyo4 !  outer
    dwcdyi5 = dwcdyi5+el_5_5*jacinv_5_5*dwdyi5 !  inner
    dwcdyo5 = dwcdyo5+el_5_5*jacinv_5_5*dwdyo5 !  outer
!
!   Enforce LODI relations
    if (ev1>0._mykind) then ! Waves pointing out of the domain
     dwcdy1 = dwcdyi1
    else ! Waves entering the domain
     dwcdy1 = ibcnr_gpu(4)*dwcdyo1 ! n.r with or without relaxation
    endif
    if (ev2>0._mykind) then ! Waves pointing out of the domain
     dwcdy2 = dwcdyi2
    else ! Waves entering the domain
     dwcdy2 = ibcnr_gpu(4)*dwcdyo2 ! n.r with or without relaxation
    endif
    if (ev3>0._mykind) then ! Waves pointing out of the domain
     dwcdy3 = dwcdyi3
    else ! Waves entering the domain
     dwcdy3 = ibcnr_gpu(4)*dwcdyo3 ! n.r with or without relaxation
    endif
    if (ev4>0._mykind) then ! Waves pointing out of the domain
     dwcdy4 = dwcdyi4
    else ! Waves entering the domain
     dwcdy4 = ibcnr_gpu(4)*dwcdyo4 ! n.r with or without relaxation
    endif
    if (ev5>0._mykind) then ! Waves pointing out of the domain
     dwcdy5 = dwcdyi5
    else ! Waves entering the domain
     dwcdy5 = ibcnr_gpu(4)*dwcdyo5 ! n.r with or without relaxation
    endif
!   Amplitude of characteristic waves
    dwcdy1 = dwcdy1 * ev1
    dwcdy2 = dwcdy2 * ev2
    dwcdy3 = dwcdy3 * ev3
    dwcdy4 = dwcdy4 * ev4
    dwcdy5 = dwcdy5 * ev5
!   Exploit LODI relations (just in case)
!   if (ibcnr_gpu(4)==2) then ! Poinsot & Lele subsonic inflow treatment 
!   endif
!
!   Return to conservative variables 
    !m=1
    df = 0._mykind
    !mm=1
    df = df + jac_1_1*er_1_1*dwcdy1
    df = df + jac_1_1*er_1_2*dwcdy2
    df = df + jac_1_1*er_1_3*dwcdy3
    df = df + jac_1_1*er_1_4*dwcdy4
    df = df + jac_1_1*er_1_5*dwcdy5
    !mm=2
    df = df + jac_1_2*er_2_1*dwcdy1
    df = df + jac_1_2*er_2_2*dwcdy2
    df = df + jac_1_2*er_2_3*dwcdy3
    df = df + jac_1_2*er_2_4*dwcdy4
    df = df + jac_1_2*er_2_5*dwcdy5
    !mm=3
    df = df + jac_1_3*er_3_1*dwcdy1
    df = df + jac_1_3*er_3_2*dwcdy2
    df = df + jac_1_3*er_3_3*dwcdy3
    df = df + jac_1_3*er_3_4*dwcdy4
    df = df + jac_1_3*er_3_5*dwcdy5
    !mm=4
    df = df + jac_1_4*er_4_1*dwcdy1
    df = df + jac_1_4*er_4_2*dwcdy2
    df = df + jac_1_4*er_4_3*dwcdy3
    df = df + jac_1_4*er_4_4*dwcdy4
    df = df + jac_1_4*er_4_5*dwcdy5
    !mm=5
    df = df + jac_1_5*er_5_1*dwcdy1
    df = df + jac_1_5*er_5_2*dwcdy2
    df = df + jac_1_5*er_5_3*dwcdy3
    df = df + jac_1_5*er_5_4*dwcdy4
    df = df + jac_1_5*er_5_5*dwcdy5
    fl_gpu(i,j,k,1) = fl_gpu(i,j,k,1) + df * detady_gpu(j)
    !m=2
    df = 0._mykind
    !mm=1
    df = df + jac_2_1*er_1_1*dwcdy1
    df = df + jac_2_1*er_1_2*dwcdy2
    df = df + jac_2_1*er_1_3*dwcdy3
    df = df + jac_2_1*er_1_4*dwcdy4
    df = df + jac_2_1*er_1_5*dwcdy5
    !mm=2
    df = df + jac_2_2*er_2_1*dwcdy1
    df = df + jac_2_2*er_2_2*dwcdy2
    df = df + jac_2_2*er_2_3*dwcdy3
    df = df + jac_2_2*er_2_4*dwcdy4
    df = df + jac_2_2*er_2_5*dwcdy5
    !mm=3
    df = df + jac_2_3*er_3_1*dwcdy1
    df = df + jac_2_3*er_3_2*dwcdy2
    df = df + jac_2_3*er_3_3*dwcdy3
    df = df + jac_2_3*er_3_4*dwcdy4
    df = df + jac_2_3*er_3_5*dwcdy5
    !mm=4
    df = df + jac_2_4*er_4_1*dwcdy1
    df = df + jac_2_4*er_4_2*dwcdy2
    df = df + jac_2_4*er_4_3*dwcdy3
    df = df + jac_2_4*er_4_4*dwcdy4
    df = df + jac_2_4*er_4_5*dwcdy5
    !mm=5
    df = df + jac_2_5*er_5_1*dwcdy1
    df = df + jac_2_5*er_5_2*dwcdy2
    df = df + jac_2_5*er_5_3*dwcdy3
    df = df + jac_2_5*er_5_4*dwcdy4
    df = df + jac_2_5*er_5_5*dwcdy5
    fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) + df * detady_gpu(j)
    !m=3
    df = 0._mykind
    !mm=1
    df = df + jac_3_1*er_1_1*dwcdy1
    df = df + jac_3_1*er_1_2*dwcdy2
    df = df + jac_3_1*er_1_3*dwcdy3
    df = df + jac_3_1*er_1_4*dwcdy4
    df = df + jac_3_1*er_1_5*dwcdy5
    !mm=2
    df = df + jac_3_2*er_2_1*dwcdy1
    df = df + jac_3_2*er_2_2*dwcdy2
    df = df + jac_3_2*er_2_3*dwcdy3
    df = df + jac_3_2*er_2_4*dwcdy4
    df = df + jac_3_2*er_2_5*dwcdy5
    !mm=3
    df = df + jac_3_3*er_3_1*dwcdy1
    df = df + jac_3_3*er_3_2*dwcdy2
    df = df + jac_3_3*er_3_3*dwcdy3
    df = df + jac_3_3*er_3_4*dwcdy4
    df = df + jac_3_3*er_3_5*dwcdy5
    !mm=4
    df = df + jac_3_4*er_4_1*dwcdy1
    df = df + jac_3_4*er_4_2*dwcdy2
    df = df + jac_3_4*er_4_3*dwcdy3
    df = df + jac_3_4*er_4_4*dwcdy4
    df = df + jac_3_4*er_4_5*dwcdy5
    !mm=5
    df = df + jac_3_5*er_5_1*dwcdy1
    df = df + jac_3_5*er_5_2*dwcdy2
    df = df + jac_3_5*er_5_3*dwcdy3
    df = df + jac_3_5*er_5_4*dwcdy4
    df = df + jac_3_5*er_5_5*dwcdy5
    fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) + df * detady_gpu(j)
    !m=4
    df = 0._mykind
    !mm=1
    df = df + jac_4_1*er_1_1*dwcdy1
    df = df + jac_4_1*er_1_2*dwcdy2
    df = df + jac_4_1*er_1_3*dwcdy3
    df = df + jac_4_1*er_1_4*dwcdy4
    df = df + jac_4_1*er_1_5*dwcdy5
    !mm=2
    df = df + jac_4_2*er_2_1*dwcdy1
    df = df + jac_4_2*er_2_2*dwcdy2
    df = df + jac_4_2*er_2_3*dwcdy3
    df = df + jac_4_2*er_2_4*dwcdy4
    df = df + jac_4_2*er_2_5*dwcdy5
    !mm=3
    df = df + jac_4_3*er_3_1*dwcdy1
    df = df + jac_4_3*er_3_2*dwcdy2
    df = df + jac_4_3*er_3_3*dwcdy3
    df = df + jac_4_3*er_3_4*dwcdy4
    df = df + jac_4_3*er_3_5*dwcdy5
    !mm=4
    df = df + jac_4_4*er_4_1*dwcdy1
    df = df + jac_4_4*er_4_2*dwcdy2
    df = df + jac_4_4*er_4_3*dwcdy3
    df = df + jac_4_4*er_4_4*dwcdy4
    df = df + jac_4_4*er_4_5*dwcdy5
    !mm=5
    df = df + jac_4_5*er_5_1*dwcdy1
    df = df + jac_4_5*er_5_2*dwcdy2
    df = df + jac_4_5*er_5_3*dwcdy3
    df = df + jac_4_5*er_5_4*dwcdy4
    df = df + jac_4_5*er_5_5*dwcdy5
    fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) + df * detady_gpu(j)
    !m=5
    df = 0._mykind
    !mm=1
    df = df + jac_5_1*er_1_1*dwcdy1
    df = df + jac_5_1*er_1_2*dwcdy2
    df = df + jac_5_1*er_1_3*dwcdy3
    df = df + jac_5_1*er_1_4*dwcdy4
    df = df + jac_5_1*er_1_5*dwcdy5
    !mm=2
    df = df + jac_5_2*er_2_1*dwcdy1
    df = df + jac_5_2*er_2_2*dwcdy2
    df = df + jac_5_2*er_2_3*dwcdy3
    df = df + jac_5_2*er_2_4*dwcdy4
    df = df + jac_5_2*er_2_5*dwcdy5
    !mm=3
    df = df + jac_5_3*er_3_1*dwcdy1
    df = df + jac_5_3*er_3_2*dwcdy2
    df = df + jac_5_3*er_3_3*dwcdy3
    df = df + jac_5_3*er_3_4*dwcdy4
    df = df + jac_5_3*er_3_5*dwcdy5
    !mm=4
    df = df + jac_5_4*er_4_1*dwcdy1
    df = df + jac_5_4*er_4_2*dwcdy2
    df = df + jac_5_4*er_4_3*dwcdy3
    df = df + jac_5_4*er_4_4*dwcdy4
    df = df + jac_5_4*er_4_5*dwcdy5
    !mm=5
    df = df + jac_5_5*er_5_1*dwcdy1
    df = df + jac_5_5*er_5_2*dwcdy2
    df = df + jac_5_5*er_5_3*dwcdy3
    df = df + jac_5_5*er_5_4*dwcdy4
    df = df + jac_5_5*er_5_5*dwcdy5
    fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) + df * detady_gpu(j)
   enddo  ! end of i-loop 
  enddo  ! end of k-loop
  !@cuf iercuda=cudaDeviceSynchronize()
!
 elseif (ilat==5) then  ! back side
 elseif (ilat==6) then  ! fore side
 endif
!
end subroutine bcrelax
