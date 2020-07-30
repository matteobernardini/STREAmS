subroutine readinp
!
! Reading the input file
!
 use mod_streams
 implicit none
!
 integer :: l, i_skip
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
 read (12,*) rm, retauinflow, trat, visc_type, s2tinf, dftscaling
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
  ibc(1) = 1
  ibc(2) = 4
  ibc(3) = 4
  ibc(4) = 4
! ibc(1) = 2
! ibc(2) = 2
! ibc(3) = 2
! ibc(4) = 2
 case (0)
  ibc(3) = 5
  ibc(4) = 5
 case (1)
  ibc(1) = 9
  ibc(2) = 4
  ibc(3) = 8
  ibc(4) = 4
 case (2)
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
