subroutine compute_coeff_xyz_midpi()
!
! Compute coefficients for mid-point interpolation on non equally-spaced grids
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,mm,l,ll
 real(mykind) :: xx,xxx
 real(mykind) :: yy,yyy
 real(mykind) :: zz,zzz
 real(mykind),dimension(:,:),allocatable ::amat
!
 allocate(amat(ivis,ivis))
!
 cx_midpi = 0._mykind
 cy_midpi = 0._mykind
 cz_midpi = 0._mykind
!
! x direction
!
 mm = ivis/2
 do i=0,nx
  xx = xh(i)
  do l=-mm+1,mm
   xxx = x(i+l)
   do ll=1,ivis
    amat(l+mm,ll) = (xxx-xx)**(ll-1)
   enddo
  enddo
  call invmat(amat,ivis)
  do l=1-mm,mm
   cx_midpi(i,l) = amat(1,mm+l)
  enddo
 enddo
!
! y direction
!
 mm = ivis/2
 do j=0,ny
  yy  = yh(j)
  yy = 0.5_mykind*(y(j)+y(j+1))
  do l=-mm+1,mm
   yyy = y(j+l)
   do ll=1,ivis
    amat(l+mm,ll) = (yyy-yy)**(ll-1)
   enddo
  enddo
  call invmat(amat,ivis)
  do l=1-mm,mm
   cy_midpi(j,l) = amat(1,mm+l)
  enddo
 enddo
!
! z direction
!
 mm = ivis/2
 do k=0,nz
  zz = zh(k)
  do l=-mm+1,mm
   zzz = z(k+l)
   do ll=1,ivis
    amat(l+mm,ll) = (zzz-zz)**(ll-1)
   enddo
  enddo
  call invmat(amat,ivis)
  do l=1-mm,mm
   cz_midpi(k,l) = amat(1,mm+l)
  enddo
 enddo
!
 deallocate(amat)
!
end subroutine compute_coeff_xyz_midpi
