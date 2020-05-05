subroutine updateghost
!
! Update ghost nodes
!
 use mod_streams
 implicit none
!
 call bc(0)
 call bcswap_prepare()
 call bcswap()
!
end subroutine updateghost
