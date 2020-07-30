subroutine step
!
! Evaluation of the time step
!
 use mod_streams 
 implicit none
!
 integer :: i,j,k
 real(mykind) :: al,c,cc,dttemp,dtxi,dtyi,dtzi,evmax
 real(mykind) :: ri,rmu,rnu,tt,uu,vv,ww,rho,rhoe,qq,pp
 real(mykind), dimension(ny) :: evmax_mat_y_cpu
 
 evmax = 0._mykind
 !$cuf kernel do(1) <<<*,*>>> 
 do j=1,ny
  evmax_mat_y(j) = 0._mykind
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
 !$cuf kernel do(2) <<<*,*>>> 
 do j=1,ny
  do k=1,nz
   evmax_mat_yz(j,k) = 0._mykind
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()

 !$cuf kernel do(2) <<<*,*>>> 
 do k=1,nz
  do j=1,ny
   do i=1,nx
    rho  = w_gpu(i,j,k,1)
    ri   = 1._mykind/rho
    uu   = w_gpu(i,j,k,2)*ri
    vv   = w_gpu(i,j,k,3)*ri
    ww   = w_gpu(i,j,k,4)*ri
    rhoe = w_gpu(i,j,k,5)
    qq   = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
    pp   = gm1*(rhoe-rho*qq)
    tt   = pp*ri ! Note that temperature (i,j,k) is still the old one
    if (visc_type==1) then
     rmu  = sqgmr*tt**vtexp
    else
     rmu  = sqgmr*sqrt(tt)*(1._mykind+s2tinf)/(1._mykind+s2tinf/tt)
    endif
    rnu  = ri*rmu
    al   = rnu*ggmopr  ! molecular conductivity
    cc   = gamma*tt
    c    = sqrt(cc)
    dtxi = (abs(uu)+c)*dcsidx_gpu(i) 
    dtyi = (abs(vv)+c)*detady_gpu(j) 
    dtzi = (abs(ww)+c)*dzitdz_gpu(k) 
    evmax_mat_yz(j,k) = max(dtxi,evmax_mat_yz(j,k))
    evmax_mat_yz(j,k) = max(dtyi,evmax_mat_yz(j,k))
    evmax_mat_yz(j,k) = max(dtzi,evmax_mat_yz(j,k))
    dtxi  = rnu*dcsidxs_gpu(i)
    dtyi  = rnu*detadys_gpu(j)
    dtzi  = rnu*dzitdzs_gpu(k)
    evmax_mat_yz(j,k) = max(dtxi,evmax_mat_yz(j,k))
    evmax_mat_yz(j,k) = max(dtyi,evmax_mat_yz(j,k))
    evmax_mat_yz(j,k) = max(dtzi,evmax_mat_yz(j,k))
    dtxi  = al*dcsidxs_gpu(i)
    dtyi  = al*detadys_gpu(j)
    dtzi  = al*dzitdzs_gpu(k)
    evmax_mat_yz(j,k) = max(dtxi,evmax_mat_yz(j,k))
    evmax_mat_yz(j,k) = max(dtyi,evmax_mat_yz(j,k))
    evmax_mat_yz(j,k) = max(dtzi,evmax_mat_yz(j,k))
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()

 !$cuf kernel do(1) <<<*,*>>> 
 do j=1,ny
  do k=1,nz
   evmax_mat_y(j) = max(evmax_mat_yz(j,k), evmax_mat_y(j))
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 evmax_mat_y_cpu = evmax_mat_y
!
 evmax = maxval(evmax_mat_y_cpu)
!
 dttemp = 1._mykind/evmax 
!
!MPI global time step evaluation
!
 call mpi_allreduce(dttemp,dtmin,1,mpi_prec,mpi_min,mp_cart,iermpi)
!
end subroutine step
!
subroutine checkdt
!
! Evaluation of the time step
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,l
 real(mykind) :: al,c,cc,dtxi,dtyi,dtzi
 real(mykind) :: evmax_i_c,evmax_i_t,evmax_i_v
 real(mykind) :: evmax_j_c,evmax_j_t,evmax_j_v
 real(mykind) :: evmax_k_c,evmax_k_t,evmax_k_v
 real(mykind) :: ri,rmu,rnu,tt,uu,vv,ww,rho,rhoe,qq,pp
 real(mykind), dimension(9) :: vecdt,vecdtmin
!
 evmax_i_c = 0._mykind
 evmax_j_c = 0._mykind
 evmax_k_c = 0._mykind
 evmax_i_v = 0._mykind
 evmax_j_v = 0._mykind
 evmax_k_v = 0._mykind
 evmax_i_t = 0._mykind
 evmax_j_t = 0._mykind
 evmax_k_t = 0._mykind
 do k=1,nz
  do j=1,ny
   do i=1,nx
    rho  = w(1,i,j,k)
    ri   = 1._mykind/rho
    uu   = w(2,i,j,k)*ri
    vv   = w(3,i,j,k)*ri
    ww   = w(4,i,j,k)*ri
    rhoe = w(5,i,j,k)
    qq   = 0.5_mykind*(uu*uu+vv*vv+ww*ww)
    pp   = gm1*(rhoe-rho*qq)
    tt   = pp*ri ! Note that temperature (i,j,k) is still the old one
    if (visc_type==1) then
     rmu  = sqgmr*tt**vtexp
    else
     rmu  = sqgmr*sqrt(tt)*(1._mykind+s2tinf)/(1._mykind+s2tinf/tt)
    endif
    rnu  = ri*rmu
    al   = rnu*ggmopr  ! molecular conductivity
    cc   = gamma*tt
    c    = sqrt(cc)
    dtxi = (abs(uu)+c)*dcsidx(i) 
    dtyi = (abs(vv)+c)*detady(j) 
    dtzi = (abs(ww)+c)*dzitdz(k) 
    evmax_i_c = max(dtxi,evmax_i_c)
    evmax_j_c = max(dtyi,evmax_j_c)
    evmax_k_c = max(dtzi,evmax_k_c)
    dtxi  = rnu*dcsidxs(i)
    dtyi  = rnu*detadys(j)
    dtzi  = rnu*dzitdzs(k)
    evmax_i_v = max(dtxi,evmax_i_v)
    evmax_j_v = max(dtyi,evmax_j_v)
    evmax_k_v = max(dtzi,evmax_k_v)
    dtxi  = al*dcsidxs(i)
    dtyi  = al*detadys(j)
    dtzi  = al*dzitdzs(k)
    evmax_i_t = max(dtxi,evmax_i_t)
    evmax_j_t = max(dtyi,evmax_j_t)
    evmax_k_t = max(dtzi,evmax_k_t)
   enddo
  enddo
 enddo
!
 vecdt = 1000000._mykind
 vecdt(1) = 1._mykind/evmax_i_c
 vecdt(2) = 1._mykind/evmax_j_c
 vecdt(4) = 1._mykind/evmax_i_v
 vecdt(5) = 1._mykind/evmax_j_v
 vecdt(7) = 1._mykind/evmax_i_t
 vecdt(8) = 1._mykind/evmax_j_t
 if (ndim==3) then
  vecdt(3) = 1._mykind/evmax_k_c
  vecdt(6) = 1._mykind/evmax_k_v
  vecdt(9) = 1._mykind/evmax_k_t
 endif
!
! MPI global time step evaluation
!
 call mpi_allreduce(vecdt,vecdtmin,9,mpi_prec,mpi_min,mp_cart,iermpi)
!
 if (masterproc) then
  open(12,file='dt_limitations')
  write(12,100) (vecdtmin(l),l=1,9)
 100  format(20ES20.10)
  close(12)
 endif 
!
end subroutine checkdt
