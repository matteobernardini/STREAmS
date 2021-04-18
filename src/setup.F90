subroutine setup
!
! Setup for the computation
!
 use mod_streams
 implicit none
!
!===================================================
 if (masterproc) write(error_unit,*) 'Allocation of variables'
 call allocate_vars()
!===================================================
 if (masterproc) write(error_unit,*) 'Reading input'
 call readinp()
!===================================================
 if (idiski==0) then
  if (masterproc) write(*,*) 'Generating mesh'
  call generategrid()
 else
  if (masterproc) write(*,*) 'Reading mesh'
  call readgrid()
 endif
 if (masterproc) write(*,*) 'Computing metrics'
 call computemetrics()
 if (enable_plot3d>0) call writegridplot3d()
!===================================================
 call constants()
!===================================================
 if (xrecyc>0._mykind) call recyc_prepare()
!===================================================
 if (masterproc) write(*,*) 'Initialize field'
 call init()
 call checkdt()
!===================================================
!
end subroutine setup
