subroutine bc(inr)
!
! Application of boundary conditions
!
 use mod_streams
 implicit none
!
 integer :: inr,ilat
!
! Table of values for parameter (ibc(ilat))
! ilat = 1,6
!
!     ____________           ibc = 1 -> 
!    /|         /|           ibc = 2 -> extrapolation 
!   / |  4     / |           ibc = 3 -> 
!  /  |    5  /  |  j        ibc = 4 -> nonreflecting
! /__________/   |  |        ibc = 5 -> wall (staggered)
! | 1 |______|_2 |  |____ i  ibc = 6 -> wall
! |  / 6     |  /  /         ibc = 7 -> oblique shock imposed
! | /    3   | /  /          ibc = 8 -> wall (PL type bc)
! |/         |/  k           ibc = 9 -> digital filtering for turbulent inflow
! /----------/
!
! inr = 0 -> steady-type    boundary conditions
! inr = 1 -> non-reflecting boundary conditions
!
  if (inr==0) then 
!
!  Steady-type BCs
!
!  'Physical' boundary conditions
!
   do ilat=1,2*ndim ! loop on all sides of the boundary (3D -> 6, 2D -> 4)
    if (ibc(ilat)==1) call bcfree(ilat)
    if (ibc(ilat)==2) call bcextr(ilat)
    if (ibc(ilat)==4) call bcextr(ilat)
    if (ibc(ilat)==5) call bcwall_staggered(ilat)
    if (ibc(ilat)==6) call bcwall(ilat)
    if (ibc(ilat)==7) call bcshk(ilat)
    if (ibc(ilat)==8) call bcwall(ilat)
    if (ibc(ilat)==9) then
     if (.not.dfupdated) then
      call bcdf(ilat)
      dfupdated = .true.
     endif
    endif
   enddo
!
  else
!
!  Unsteady-type BCs (update boundary fluxes)
!
   do ilat=1,2*ndim ! loop on all sides of the boundary (3D -> 6, 2D -> 4)
    if (ibc(ilat)==4) call bcrelax  (ilat)
    if (ibc(ilat)==7) call bcrelax  (ilat)
    if (ibc(ilat)==8) call bcwall_pl(ilat)
   enddo
  endif
!  
end subroutine bc
