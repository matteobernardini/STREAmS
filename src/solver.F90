subroutine solver
!
! Solve the compressible NS equations
!
 use mod_streams
 implicit none
!
 real(mykind) :: elapsed,startTiming,endTiming
!
 integer :: i
 logical :: stop_streams
!
 icyc   = ncyc0
 telaps = telaps0
!
!Copy arrays from CPU to GPU or move alloc
 call copy_cpu_to_gpu()
!
 if (masterproc) write(*,*) 'Compute time step'
 dtmin = abs(cfl)
 if (cfl>0._mykind) call step()
 if (masterproc) write(*,*) 'Done'
!
 call updateghost() ! Needed here only for subsequent prims call
 call prims()
 if (tresduc>0._mykind.and.tresduc<1._mykind) then
  call sensor()
  call bcswapduc_prepare()
  call bcswapduc()
 endif
!
 open(20,file='output_streams.dat',position=stat_io)
 startTiming = mpi_wtime()
!
 stop_streams = .false.
 do i=1,ncyc
!
  icyc = icyc+1
!
  call rk() ! Third-order RK scheme
!
  if (io_type>0) call manage_solver()
!
  if (mod(i,nprint)==0) then
   call computeresidual()
   call printres()
  endif
!
  if (cfl>0._mykind) then
   if (mod(i,nstep)==0) call step() ! Compute the time step
  endif
!
  call mpi_barrier(mpi_comm_world,iermpi)
  inquire(file="stop.stop",exist=stop_streams)
  call mpi_barrier(mpi_comm_world,iermpi)
  if (stop_streams) exit
!
 enddo
!
 endTiming = mpi_wtime()
 elapsed = endTiming-startTiming
 if (ncyc>0) then
  if (masterproc) write(error_unit,*) 'Time-step time =', elapsed/ncyc
  if (masterproc) write(20,*) 'Time-step time =', elapsed/ncyc
 endif
!
end subroutine solver
