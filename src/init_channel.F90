subroutine init_channel
!
! Apply wmean and perturbations for channel flow simulation
!
 use mod_streams
 implicit none
!
 integer :: i,j,k
 real(mykind), dimension(3) :: rr
 real(mykind) :: u0_02,u0_05
 real(mykind) :: rho,uu,vv,ww,ufluc,vfluc,wfluc,rhouu,rhovv,rhoww
!
 u0_02 = 0.02_mykind*u0
 u0_05 = 0.05_mykind*u0
!
 do k=1,nz
  do j=1,ny
   do i=1,nx
!
    rho = wmean(1,i,j)
    uu  = wmean(2,i,j)/rho
    vv  = wmean(3,i,j)/rho
    ww  = wmean(4,i,j)/rho

    call random_number(rr)
    rr = rr-0.5_mykind

    !rr = 0._mykind

    ufluc   = u0_02*rr(1)
    vfluc   = u0_02*rr(2)
    wfluc   = u0_02*rr(3)
    vfluc   = vfluc+u0_05*sin(0.5_mykind*pi*y(j))*cos(2*pi*z(k))
    wfluc   = wfluc+u0_05*sin(0.5_mykind*pi*y(j))*sin(2*pi*z(k))
    uu      = uu+ufluc
    vv      = vv+vfluc
    ww      = ww+wfluc
    rhouu   = rho*uu
    rhovv   = rho*vv
    rhoww   = rho*ww
    w(1,i,j,k) = rho
    w(2,i,j,k) = rhouu
    w(3,i,j,k) = rhovv
    w(4,i,j,k) = rhoww
    w(5,i,j,k) = p0*gm+0.5_mykind*(rhouu**2+rhovv**2+rhoww**2)/rho
   enddo
  enddo
 enddo
!
end subroutine init_channel
