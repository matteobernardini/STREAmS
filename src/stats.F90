subroutine stats2d
!
! Updating statistical information
!
 use mod_streams
 implicit none
!
 integer :: i,j,l
 real(mykind) :: fnsm1,fns
!
 if (masterproc) write(*,*) 'Computing statistics, itav =',itav+1
!
 call compute_av
!
 itav  = itav+1
 fnsm1 = real(itav-1,mykind)
 fns   = real(itav  ,mykind)
 do j=1,ny
  do i=1,nx
   do l=1,nvmean
    w_av(l,i,j) = (w_av(l,i,j)*fnsm1+w_avzg(l,i,j))/fns
   enddo
  enddo
 enddo
!
end subroutine stats2d
!
subroutine stats1d
!
! Updating statistical information
!
 use mod_streams
 implicit none
!
 integer :: j,l
 real(mykind) :: fnsm1,fns
!
 if (masterproc) write(*,*) 'Computing statistics, itav=',itav+1
!
 call compute_av1d
!
 itav  = itav+1
 fnsm1 = real(itav-1,mykind)
 fns   = real(itav  ,mykind)
 do j=1,ny
  do l=1,nvmean
   w_av_1d(l,j) = (w_av_1d(l,j)*fnsm1+w_avxzg(l,j))/fns
  enddo
 enddo
!
end subroutine stats1d
