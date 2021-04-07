subroutine readinp
!
! Reading the input file
!
 use mod_streams
 implicit none
!
 integer :: l, i_skip
!
!
! iflow = 0 ==> Channel flow
!
!       rm:             Mach number based on bulk velocity and wall temperature
!       retauinflow:    Estimated friction Reynolds number (target)
!       trat:           Bulk temperature/wall temperature (if <=0., bulk temperature free to evolve)
!       pgradf:         Assigned pressure gradient (if 0. pressure gradient computed to maintain bulk velocity constant)
!
! iflow = 1,2 ==> Boundary layer,SBLI
!
!       rm:             Free-stream Mach number
!       retauinflow:    Estimated friction Reynolds number (target) at inflow
!       trat:           Wall-to-recovery-temperature ratio
!       deflec:         Deflection angle after shock (only for SBLI)
!
 open (unit=12,file='input.dat',form='formatted')
 do i_skip=1,33
  read (12,*)
 enddo
 read (12,*)
 read (12,*) tresduc, ximp, deflec, pgradf
 read (12,*)
 read (12,*)
 read (12,*) idiski, ncyc, cfl, nstep, nprint, io_type
 read (12,*)
 read (12,*)
 read (12,*) rm, retauinflow, trat, visc_type, tref_dimensional, dftscaling
 read (12,*)
 read (12,*)
 read (12,*) istat, nstat
 allocate( xstat(nstat))
 allocate(ixstat(nstat))
 allocate(igxstat(nstat))
 read (12,*)
 read (12,*)
 read (12,*) (xstat(l),l=1,nstat)
 read (12,*)
 read (12,*)
 read (12,*) dtsave, dtsave_restart, enable_plot3d, enable_vtk
 read (12,*)
 read (12,*)
 read (12,*) rand_start
 close(12)
!
 call check_input(2)
!
 !ivis   = 4
 !iorder = 4
 !iweno  = 2
!
 ibc   = 0
 ibcnr = 0
 select case (iflow)
 case (-1)
  ibc(1) = 1 ! Free stream
  ibc(2) = 4 ! Extrapolation + non reflecting treatment
  ibc(3) = 2 ! Extrapolation
  ibc(4) = 2
  ibc(5) = 2
  ibc(6) = 2
 case (0) ! Channel
  ibc(3) = 5 ! Staggered wall
  ibc(4) = 5
 case (1) ! BL
  ibc(1) = 9 ! Digital filtering
  ibc(2) = 4 ! Extrapolation + non reflecting treatment
  ibc(3) = 8 ! Wall + reflecting treatment
  ibc(4) = 4
 case (2) ! SBLI
  ibc(1) = 9
  ibc(2) = 4
  ibc(3) = 8
  ibc(4) = 7
 end select
!
 ndim = 3 ! default number of dimensions
 if (nzmax==1) ndim = 2
!
 do l=0,nsolmax
  tsol(l)         = l*dtsave
  tsol_restart(l) = l*dtsave_restart
 enddo
!
end subroutine readinp
