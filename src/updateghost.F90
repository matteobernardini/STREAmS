subroutine updateghost
!
! Update ghost nodes
!
 use mod_streams
 implicit none
!
 call bcswap_prepare()
 call bc(0)
 call bcswap()
!
end subroutine updateghost
