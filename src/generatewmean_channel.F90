subroutine generatewmean_channel
!
! Generating field for channel initialization
!
 use mod_streams
 implicit none
!
 integer :: i,j
!
 wmean = 0._mykind
 do i=1,nx
  do j=1,ny
   wmean(1,i,j) = rho0
   wmean(2,i,j) = 1.5_mykind*rho0*u0*(1._mykind-y(j)**2)
   wmean(3,i,j) = 0._mykind
   wmean(4,i,j) = 0._mykind
   wmean(5,i,j) = p0*gm+0.5_mykind*(wmean(2,i,j)**2)/wmean(1,i,j)
  enddo
 enddo
!
end subroutine generatewmean_channel
