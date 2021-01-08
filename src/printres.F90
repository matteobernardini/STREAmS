subroutine printres
!
! Printing residual
!
 use mod_streams
 implicit none
!
 real(mykind) :: sumq
!
 if (masterproc) then
  if (iflow==0) then
   sumq = heatflux+aeroheat+bulkcooling
   write(* ,100) icyc,telaps,rtrms,dpdx,rhobulk,ubulk,tbulk,heatflux,aeroheat,bulkcooling,sumq
   write(20,100) icyc,telaps,rtrms,dpdx,rhobulk,ubulk,tbulk,heatflux,aeroheat,bulkcooling,sumq
  else
   write(* ,100) icyc,telaps,rtrms
   write(20,100) icyc,telaps,rtrms
  endif
 endif
!
 100 format(1I10,40ES20.10)
!
end subroutine printres
