subroutine setup
!
! Setup for the computation
!
 use mod_streams
 implicit none
!
!===================================================
 if (masterproc) write(*,*) 'Allocation of variables'
 call allocate_vars()
!===================================================
 if (masterproc) write(*,*) 'Reading input'
 call readinp()
!===================================================
 if (idiski==0) then
  if (masterproc) write(*,*) 'Generating mesh'
  call generategrid()
 else
  if (masterproc) write(*,*) 'Reading mesh'
  call readgrid()
 endif
 call computemetrics()
 if (idiski==0 .and. enable_plot3d > 0) call writegridplot3d()
!===================================================
 call constants()
!===================================================
 if (masterproc) write(*,*) 'Initialize field'
 call init()
 call checkdt()
!===================================================
!
end subroutine setup
