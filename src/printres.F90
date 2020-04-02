subroutine printres
!
! Printing residual
!
 use mod_streams
 implicit none
!
 if (masterproc) then
  if (iflow==0) then
   write(* ,100) icyc,telaps,rtrms,dpdx,rhobulk,ubulk,tbulk
   write(20,100) icyc,telaps,rtrms,dpdx,rhobulk,ubulk,tbulk
  else
   write(* ,100) icyc,telaps,rtrms
   write(20,100) icyc,telaps,rtrms
  endif
 endif
!
 100 format(1I10,40ES20.10)
!
end subroutine printres
