subroutine compute_av
!
! Compute spanwise average
!
 use mod_streams
 implicit none
!
 real(mykind), dimension(nvmean,nx,ny) :: w_avz
 real(mykind) :: pp,pp2,rho,rho2,rhoe,rhou,rhov,rhow
 real(mykind) :: ri,rlam,rmu
 real(mykind) :: tt,tt2
 real(mykind) :: uu,vv,ww
 real(mykind) :: uu2,vv2,ww2
 integer :: i,j,k,npt
!
 w_avz = 0._mykind
!
 do k=1,nz
  do j=1,ny
   do i=1,nx
    rho  = w(1,i,j,k)
    rhou = w(2,i,j,k)
    rhov = w(3,i,j,k)
    rhow = w(4,i,j,k)
    rhoe = w(5,i,j,k)
    ri   = 1._mykind/rho
    uu   = rhou*ri
    vv   = rhov*ri
    ww   = rhow*ri
    pp   = rho*temperature(i,j,k)
    tt   = pp*ri
    if (visc_type==1) then
     rmu  = sqgmr*tt**vtexp
    else
     rmu  = sqgmr*sqrt(tt)*(1._mykind+s2tinf)/(1._mykind+s2tinf/tt)
    endif
    rho2 = rho*rho
    uu2  = uu*uu
    vv2  = vv*vv
    ww2  = ww*ww
    pp2  = pp*pp
    tt2  = tt*tt
!
    w_avz(1,i,j) = w_avz(1,i,j)+rho
    w_avz(2,i,j) = w_avz(2,i,j)+uu
    w_avz(3,i,j) = w_avz(3,i,j)+vv
    w_avz(4,i,j) = w_avz(4,i,j)+ww
    w_avz(5,i,j) = w_avz(5,i,j)+pp
    w_avz(6,i,j) = w_avz(6,i,j)+tt
    w_avz(7,i,j) = w_avz(7,i,j)+rho2
    w_avz(8,i,j) = w_avz(8,i,j)+uu2
    w_avz(9,i,j) = w_avz(9,i,j)+vv2
    w_avz(10,i,j) = w_avz(10,i,j)+ww2
    w_avz(11,i,j) = w_avz(11,i,j)+pp2
    w_avz(12,i,j) = w_avz(12,i,j)+tt2
    w_avz(13,i,j) = w_avz(13,i,j)+rhou
    w_avz(14,i,j) = w_avz(14,i,j)+rhov
    w_avz(15,i,j) = w_avz(15,i,j)+rhow
    w_avz(16,i,j) = w_avz(16,i,j)+rhou*uu
    w_avz(17,i,j) = w_avz(17,i,j)+rhov*vv
    w_avz(18,i,j) = w_avz(18,i,j)+rhow*ww
    w_avz(19,i,j) = w_avz(19,i,j)+rhou*vv
    w_avz(20,i,j) = w_avz(20,i,j)+rmu
!
   enddo
  enddo
 enddo
!
 w_avzg = 0._mykind
 npt    = nvmean*nx*ny
 call mpi_allreduce(w_avz,w_avzg,npt,mpi_prec,mpi_sum,mp_cartz,iermpi)
 w_avzg = w_avzg/nzmax
!
end subroutine compute_av
!
subroutine compute_av1d
!
! Compute spanwise average (in 1D)
!
 use mod_streams
 implicit none
!
 real(mykind), dimension(nvmean,ny) :: w_avxz
 real(mykind) :: pp,pp2,rho,rho2,rhoe,rhou,rhov,rhow
 real(mykind) :: ri,rlam,rmu
 real(mykind) :: tt,tt2
 real(mykind) :: uu,vv,ww
 real(mykind) :: uu2,vv2,ww2
 integer :: i,j,k,npt
!
 write(error_unit,*) "compute av from rank: ",nrank
!
 w_avxz = 0._mykind
!
 do k=1,nz
  do j=1,ny
   do i=1,nx
    rho  = w(1,i,j,k)
    rhou = w(2,i,j,k)
    rhov = w(3,i,j,k)
    rhow = w(4,i,j,k)
    rhoe = w(5,i,j,k)
    ri   = 1._mykind/rho
    uu   = rhou*ri
    vv   = rhov*ri
    ww   = rhow*ri
    pp   = rho*temperature(i,j,k)
    tt   = pp*ri
    if (visc_type==1) then
     rmu  = sqgmr*tt**vtexp
    else
     rmu  = sqgmr*sqrt(tt)*(1._mykind+s2tinf)/(1._mykind+s2tinf/tt)
    endif
    rho2 = rho*rho
    uu2  = uu*uu
    vv2  = vv*vv
    ww2  = ww*ww
    pp2  = pp*pp
    tt2  = tt*tt
!
    w_avxz(1,j) = w_avxz(1,j)+rho
    w_avxz(2,j) = w_avxz(2,j)+uu
    w_avxz(3,j) = w_avxz(3,j)+vv
    w_avxz(4,j) = w_avxz(4,j)+ww
    w_avxz(5,j) = w_avxz(5,j)+pp
    w_avxz(6,j) = w_avxz(6,j)+tt
    w_avxz(7,j) = w_avxz(7,j)+rho2
    w_avxz(8,j) = w_avxz(8,j)+uu2
    w_avxz(9,j) = w_avxz(9,j)+vv2
    w_avxz(10,j) = w_avxz(10,j)+ww2
    w_avxz(11,j) = w_avxz(11,j)+pp2
    w_avxz(12,j) = w_avxz(12,j)+tt2
    w_avxz(13,j) = w_avxz(13,j)+rhou
    w_avxz(14,j) = w_avxz(14,j)+rhov
    w_avxz(15,j) = w_avxz(15,j)+rhow
    w_avxz(16,j) = w_avxz(16,j)+rhou*uu
    w_avxz(17,j) = w_avxz(17,j)+rhov*vv
    w_avxz(18,j) = w_avxz(18,j)+rhow*ww
    w_avxz(19,j) = w_avxz(19,j)+rhou*vv
    w_avxz(20,j) = w_avxz(20,j)+rmu
!
   enddo
  enddo
 enddo
!
 w_avxzg = 0._mykind
 npt    = nvmean*ny
 call mpi_allreduce(w_avxz,w_avxzg,npt,mpi_prec,mpi_sum,mp_cart,iermpi)
 w_avxzg = w_avxzg/nzmax/nxmax
!
end subroutine compute_av1d
