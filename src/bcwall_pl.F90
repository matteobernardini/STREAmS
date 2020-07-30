subroutine bcwall_pl(ilat)
!
! Apply Poinsot-Lele wall boundary conditions
!
 use mod_streams
 implicit none
!
 real(mykind) :: dwdy1,dwcdy1,dwdy2,dwcdy2,dwdy3,dwcdy3,dwdy4,dwcdy4,dwdy5,dwcdy5
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
!
 integer :: i,j,k,l,m,mm,ilat
 real(mykind) :: uu,vv,ww
 real(mykind) :: b1,b2,c,cc,ci,df,h,pp,qq,ri
!
 if (ilat==3) then ! lower side
! Redefining g(u)_y on the boundary
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    j         = 1 
    ri        =  1._mykind/w_gpu(i,j,k,1)
    uu        =  w_gpu(i,j,k,2) * ri
    vv        =  w_gpu(i,j,k,3) * ri
    ww        =  w_gpu(i,j,k,4) * ri
    qq        =  0.5_mykind * (uu*uu  + vv*vv + ww*ww)
    pp        =  gm1 * (w_gpu(i,j,k,5) - w_gpu(i,j,k,1) * qq)
    h         =  (w_gpu(i,j,k,5)  + pp) * ri
    cc        =  gamma * pp * ri
    c         =  sqrt(cc)
    ci        =  1._mykind/c 
!   left eigenvectors matrix (at Roe's state)
    b2        =   gm1/cc
    b1        =   b2 * qq
    el_1_1    =   0.5_mykind * (b1     + vv * ci)
    el_1_2    =  -0.5_mykind * (b2 * uu         )
    el_1_3    =  -0.5_mykind * (b2 * vv +     ci)
    el_1_4    =  -0.5_mykind * (b2 * ww         )
    el_1_5    =   0.5_mykind * b2
    el_2_1    =   1._mykind - b1
    el_2_2    =   b2 * uu
    el_2_3    =   b2 * vv
    el_2_4    =   b2 * ww
    el_2_5    =  -b2
    el_3_1    =   0.5_mykind * (b1     - vv * ci)
    el_3_2    =  -0.5_mykind * (b2 * uu         )
    el_3_3    =  -0.5_mykind * (b2 * vv -     ci)
    el_3_4    =  -0.5_mykind * (b2 * ww         )
    el_3_5    =   0.5_mykind * b2
    el_4_1    =  -ww
    el_4_2    =   0._mykind
    el_4_3    =   0._mykind
    el_4_4    =   1._mykind
    el_4_5    =   0._mykind
    el_5_1    =   uu
    el_5_2    =  -1._mykind
    el_5_3    =   0._mykind
    el_5_4    =   0._mykind
    el_5_5    =   0._mykind
!   right eigenvectors matrix (at Roe's state)
    er_1_1    =  1._mykind
    er_2_1    =  uu
    er_3_1    =  vv - c
    er_4_1    =  ww      
    er_5_1    =  h  - vv * c
    er_1_2    =  1._mykind
    er_2_2    =  uu
    er_3_2    =  vv
    er_4_2    =  ww
    er_5_2    =  qq
    er_1_3    =  1._mykind
    er_2_3    =  uu
    er_3_3    =  vv + c
    er_4_3    =  ww   
    er_5_3    =  h  + vv * c
    er_1_4    =  0._mykind
    er_2_4    =  0._mykind
    er_3_4    =  0._mykind
    er_4_4    =  1._mykind
    er_5_4    =  ww
    er_1_5    =  0._mykind
    er_2_5    = -1._mykind
    er_3_5    =  0._mykind
    er_4_5    =  0._mykind
    er_5_5    = -uu
!   Eigenvalues (killing the positive ones)
    ev1       =  min(vv-c,0._mykind)
    ev2       =  min(vv  ,0._mykind) 
    ev3       =  min(vv+c,0._mykind)
    ev4       =  ev2
    ev5       =  ev2
!   Derivatives of conservative variables
    dwdy1 = (-1.5_mykind*w_gpu(i,1,k,1)+2._mykind*w_gpu(i,2,k,1)-0.5_mykind*w_gpu(i,3,k,1))
    dwdy2 = (-1.5_mykind*w_gpu(i,1,k,2)+2._mykind*w_gpu(i,2,k,2)-0.5_mykind*w_gpu(i,3,k,2))
    dwdy3 = (-1.5_mykind*w_gpu(i,1,k,3)+2._mykind*w_gpu(i,2,k,3)-0.5_mykind*w_gpu(i,3,k,3))
    dwdy4 = (-1.5_mykind*w_gpu(i,1,k,4)+2._mykind*w_gpu(i,2,k,4)-0.5_mykind*w_gpu(i,3,k,4))
    dwdy5 = (-1.5_mykind*w_gpu(i,1,k,5)+2._mykind*w_gpu(i,2,k,5)-0.5_mykind*w_gpu(i,3,k,5))
!   Projection on the eigendirections
    dwcdy1 = 0._mykind
    dwcdy1 = dwcdy1 + el_1_1 * dwdy1 
    dwcdy1 = dwcdy1 + el_1_2 * dwdy2 
    dwcdy1 = dwcdy1 + el_1_3 * dwdy3 
    dwcdy1 = dwcdy1 + el_1_4 * dwdy4 
    dwcdy1 = dwcdy1 + el_1_5 * dwdy5 
    dwcdy1 = dwcdy1 * ev1 ! is the L of PL
!   dwcdy2 = 0._mykind
!   dwcdy2 = dwcdy2 + el_2_1 * dwdy1 
!   dwcdy2 = dwcdy2 + el_2_2 * dwdy2 
!   dwcdy2 = dwcdy2 + el_2_3 * dwdy3 
!   dwcdy2 = dwcdy2 + el_2_4 * dwdy4 
!   dwcdy2 = dwcdy2 + el_2_5 * dwdy5 
!   dwcdy2 = dwcdy2 * ev2 ! is the L of PL
!   dwcdy3 = 0._mykind
!   dwcdy3 = dwcdy3 + el_3_1 * dwdy1 
!   dwcdy3 = dwcdy3 + el_3_2 * dwdy2 
!   dwcdy3 = dwcdy3 + el_3_3 * dwdy3 
!   dwcdy3 = dwcdy3 + el_3_4 * dwdy4 
!   dwcdy3 = dwcdy3 + el_3_5 * dwdy5 
!   dwcdy3 = dwcdy3 * ev3 ! is the L of PL
!   dwcdy4 = 0._mykind
!   dwcdy4 = dwcdy4 + el_4_1 * dwdy1 
!   dwcdy4 = dwcdy4 + el_4_2 * dwdy2 
!   dwcdy4 = dwcdy4 + el_4_3 * dwdy3 
!   dwcdy4 = dwcdy4 + el_4_4 * dwdy4 
!   dwcdy4 = dwcdy4 + el_4_5 * dwdy5 
!   dwcdy4 = dwcdy4 * ev4 ! is the L of PL
!   dwcdy5 = 0._mykind
!   dwcdy5 = dwcdy5 + el_5_1 * dwdy1 
!   dwcdy5 = dwcdy5 + el_5_2 * dwdy2 
!   dwcdy5 = dwcdy5 + el_5_3 * dwdy3 
!   dwcdy5 = dwcdy5 + el_5_4 * dwdy4 
!   dwcdy5 = dwcdy5 + el_5_5 * dwdy5 
!   dwcdy5 = dwcdy5 * ev5 ! is the L of PL
!   wall boundary conditions
    dwcdy3 = dwcdy1 
    dwcdy2 = 0._mykind
    dwcdy4 = 0._mykind
    dwcdy5 = 0._mykind
!   Returning to conservative variables 
    df = 0._mykind
    df = df + er_1_1 * dwcdy1
    df = df + er_1_2 * dwcdy2
    df = df + er_1_3 * dwcdy3
    df = df + er_1_4 * dwcdy4
    df = df + er_1_5 * dwcdy5
    fl_gpu(i,j,k,1) = fl_gpu(i,j,k,1) + df * detady_gpu(j)
    df = 0._mykind
    df = df + er_2_1 * dwcdy1
    df = df + er_2_2 * dwcdy2
    df = df + er_2_3 * dwcdy3
    df = df + er_2_4 * dwcdy4
    df = df + er_2_5 * dwcdy5
    fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) + df * detady_gpu(j)
    df = 0._mykind
    df = df + er_3_1 * dwcdy1
    df = df + er_3_2 * dwcdy2
    df = df + er_3_3 * dwcdy3
    df = df + er_3_4 * dwcdy4
    df = df + er_3_5 * dwcdy5
    fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) + df * detady_gpu(j)
    df = 0._mykind
    df = df + er_4_1 * dwcdy1
    df = df + er_4_2 * dwcdy2
    df = df + er_4_3 * dwcdy3
    df = df + er_4_4 * dwcdy4
    df = df + er_4_5 * dwcdy5
    fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) + df * detady_gpu(j)
    df = 0._mykind
    df = df + er_5_1 * dwcdy1
    df = df + er_5_2 * dwcdy2
    df = df + er_5_3 * dwcdy3
    df = df + er_5_4 * dwcdy4
    df = df + er_5_5 * dwcdy5
    fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) + df * detady_gpu(j)
   enddo  ! end of i-loop 
  enddo  ! end of k-loop
 !@cuf iercuda=cudaDeviceSynchronize()
 elseif (ilat==4) then ! upper side
! Redefining g(u)_y on the boundary
  !$cuf kernel do(2) <<<*,*>>>
  do k=1,nz
   do i=1,nx
    j         = ny 
    ri        =  1._mykind/w_gpu(i,j,k,1)
    uu        =  w_gpu(i,j,k,2) * ri
    vv        =  w_gpu(i,j,k,3) * ri
    ww        =  w_gpu(i,j,k,4) * ri
    qq        =  0.5_mykind * (uu*uu  + vv*vv + ww*ww)
    pp        =  gm1 * (w_gpu(i,j,k,5) - w_gpu(i,j,k,1) * qq)
    h         =  (w_gpu(i,j,k,5)  + pp) * ri
    cc        =  gamma * pp * ri
    c         =  sqrt(cc)
    ci        =  1._mykind/c 
!   left eigenvectors matrix (at Roe's state)
    b2        =   gm1/cc
    b1        =   b2 * qq
    el_1_1    =   0.5_mykind * (b1     + vv * ci)
    el_1_2    =  -0.5_mykind * (b2 * uu         )
    el_1_3    =  -0.5_mykind * (b2 * vv +     ci)
    el_1_4    =  -0.5_mykind * (b2 * ww         )
    el_1_5    =   0.5_mykind * b2
    el_2_1    =   1._mykind - b1
    el_2_2    =   b2 * uu
    el_2_3    =   b2 * vv
    el_2_4    =   b2 * ww
    el_2_5    =  -b2
    el_3_1    =   0.5_mykind * (b1     - vv * ci)
    el_3_2    =  -0.5_mykind * (b2 * uu         )
    el_3_3    =  -0.5_mykind * (b2 * vv -     ci)
    el_3_4    =  -0.5_mykind * (b2 * ww         )
    el_3_5    =   0.5_mykind * b2
    el_4_1    =  -ww
    el_4_2    =   0._mykind
    el_4_3    =   0._mykind
    el_4_4    =   1._mykind
    el_4_5    =   0._mykind
    el_5_1    =   uu
    el_5_2    =  -1._mykind
    el_5_3    =   0._mykind
    el_5_4    =   0._mykind
    el_5_5    =   0._mykind
!   right eigenvectors matrix (at Roe's state)
    er_1_1    =  1._mykind
    er_2_1    =  uu
    er_3_1    =  vv - c
    er_4_1    =  ww      
    er_5_1    =  h  - vv * c
    er_1_2    =  1._mykind
    er_2_2    =  uu
    er_3_2    =  vv
    er_4_2    =  ww
    er_5_2    =  qq
    er_1_3    =  1._mykind
    er_2_3    =  uu
    er_3_3    =  vv + c
    er_4_3    =  ww   
    er_5_3    =  h  + vv * c
    er_1_4    =  0._mykind
    er_2_4    =  0._mykind
    er_3_4    =  0._mykind
    er_4_4    =  1._mykind
    er_5_4    =  ww
    er_1_5    =  0._mykind
    er_2_5    = -1._mykind
    er_3_5    =  0._mykind
    er_4_5    =  0._mykind
    er_5_5    = -uu
!   Eigenvalues (killing the negative ones)
    ev1       =  max(vv-c,0._mykind)
    ev2       =  max(vv  ,0._mykind) 
    ev3       =  max(vv+c,0._mykind)
    ev4       =  ev2
    ev5       =  ev2
!   Derivatives of conservative variables
    dwdy1 = (1.5_mykind*w_gpu(i,ny,k,1)-2._mykind*w_gpu(i,ny-1,k,1)+0.5_mykind*w_gpu(i,ny-2,k,1))
    dwdy2 = (1.5_mykind*w_gpu(i,ny,k,2)-2._mykind*w_gpu(i,ny-1,k,2)+0.5_mykind*w_gpu(i,ny-2,k,2))
    dwdy3 = (1.5_mykind*w_gpu(i,ny,k,3)-2._mykind*w_gpu(i,ny-1,k,3)+0.5_mykind*w_gpu(i,ny-2,k,3))
    dwdy4 = (1.5_mykind*w_gpu(i,ny,k,4)-2._mykind*w_gpu(i,ny-1,k,4)+0.5_mykind*w_gpu(i,ny-2,k,4))
    dwdy5 = (1.5_mykind*w_gpu(i,ny,k,5)-2._mykind*w_gpu(i,ny-1,k,5)+0.5_mykind*w_gpu(i,ny-2,k,5))
!   Projection on the eigendirections
!   dwcdy1 = 0._mykind
!   dwcdy1 = dwcdy1 + el_1_1 * dwdy1 
!   dwcdy1 = dwcdy1 + el_1_2 * dwdy2 
!   dwcdy1 = dwcdy1 + el_1_3 * dwdy3 
!   dwcdy1 = dwcdy1 + el_1_4 * dwdy4 
!   dwcdy1 = dwcdy1 + el_1_5 * dwdy5 
!   dwcdy1 = dwcdy1 * ev1 ! is the L of PL
!   dwcdy2 = 0._mykind
!   dwcdy2 = dwcdy2 + el_2_1 * dwdy1 
!   dwcdy2 = dwcdy2 + el_2_2 * dwdy2 
!   dwcdy2 = dwcdy2 + el_2_3 * dwdy3 
!   dwcdy2 = dwcdy2 + el_2_4 * dwdy4 
!   dwcdy2 = dwcdy2 + el_2_5 * dwdy5 
!   dwcdy2 = dwcdy2 * ev2 ! is the L of PL
    dwcdy3 = 0._mykind
    dwcdy3 = dwcdy3 + el_3_1 * dwdy1 
    dwcdy3 = dwcdy3 + el_3_2 * dwdy2 
    dwcdy3 = dwcdy3 + el_3_3 * dwdy3 
    dwcdy3 = dwcdy3 + el_3_4 * dwdy4 
    dwcdy3 = dwcdy3 + el_3_5 * dwdy5 
    dwcdy3 = dwcdy3 * ev3 ! is the L of PL
!   dwcdy4 = 0._mykind
!   dwcdy4 = dwcdy4 + el_4_1 * dwdy1 
!   dwcdy4 = dwcdy4 + el_4_2 * dwdy2 
!   dwcdy4 = dwcdy4 + el_4_3 * dwdy3 
!   dwcdy4 = dwcdy4 + el_4_4 * dwdy4 
!   dwcdy4 = dwcdy4 + el_4_5 * dwdy5 
!   dwcdy4 = dwcdy4 * ev4 ! is the L of PL
!   dwcdy5 = 0._mykind
!   dwcdy5 = dwcdy5 + el_5_1 * dwdy1 
!   dwcdy5 = dwcdy5 + el_5_2 * dwdy2 
!   dwcdy5 = dwcdy5 + el_5_3 * dwdy3 
!   dwcdy5 = dwcdy5 + el_5_4 * dwdy4 
!   dwcdy5 = dwcdy5 + el_5_5 * dwdy5 
!   dwcdy5 = dwcdy5 * ev5 ! is the L of PL
!   wall boundary conditions
    dwcdy1 = dwcdy3 
    dwcdy2 = 0._mykind
    dwcdy4 = 0._mykind
    dwcdy5 = 0._mykind
!   Returning to conservative variables 
    df = 0._mykind
    df = df + er_1_1 * dwcdy1
    df = df + er_1_2 * dwcdy2
    df = df + er_1_3 * dwcdy3
    df = df + er_1_4 * dwcdy4
    df = df + er_1_5 * dwcdy5
    fl_gpu(i,j,k,1) = fl_gpu(i,j,k,1) + df * detady_gpu(j)
    df = 0._mykind
    df = df + er_2_1 * dwcdy1
    df = df + er_2_2 * dwcdy2
    df = df + er_2_3 * dwcdy3
    df = df + er_2_4 * dwcdy4
    df = df + er_2_5 * dwcdy5
    fl_gpu(i,j,k,2) = fl_gpu(i,j,k,2) + df * detady_gpu(j)
    df = 0._mykind
    df = df + er_3_1 * dwcdy1
    df = df + er_3_2 * dwcdy2
    df = df + er_3_3 * dwcdy3
    df = df + er_3_4 * dwcdy4
    df = df + er_3_5 * dwcdy5
    fl_gpu(i,j,k,3) = fl_gpu(i,j,k,3) + df * detady_gpu(j)
    df = 0._mykind
    df = df + er_4_1 * dwcdy1
    df = df + er_4_2 * dwcdy2
    df = df + er_4_3 * dwcdy3
    df = df + er_4_4 * dwcdy4
    df = df + er_4_5 * dwcdy5
    fl_gpu(i,j,k,4) = fl_gpu(i,j,k,4) + df * detady_gpu(j)
    df = 0._mykind
    df = df + er_5_1 * dwcdy1
    df = df + er_5_2 * dwcdy2
    df = df + er_5_3 * dwcdy3
    df = df + er_5_4 * dwcdy4
    df = df + er_5_5 * dwcdy5
    fl_gpu(i,j,k,5) = fl_gpu(i,j,k,5) + df * detady_gpu(j)
   enddo  ! end of i-loop 
  enddo  ! end of k-loop
 !@cuf iercuda=cudaDeviceSynchronize()
 endif 
!
end subroutine bcwall_pl
