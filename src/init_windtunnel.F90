subroutine init_windtunnel
!
! Apply wmean and perturbations for channel flow simulation
!
 use mod_streams
 implicit none
!
 integer :: i,j,k
!
 do k=1-ng,nz+ng
  do j=1-ng,ny+ng
   do i=1-ng,nx+ng
    w(1,i,j,k) = winf(1)
    w(2,i,j,k) = winf(2)
    w(3,i,j,k) = winf(3)
    w(4,i,j,k) = winf(4)
    w(5,i,j,k) = winf(5)
!   w(1,i,j,k) = rho0
!   w(2,i,j,k) = 0._mykind
!   w(3,i,j,k) = 0._mykind
!   w(4,i,j,k) = 0._mykind
!   w(5,i,j,k) = gm*p0
   enddo
  enddo
 enddo
!
end subroutine init_windtunnel
