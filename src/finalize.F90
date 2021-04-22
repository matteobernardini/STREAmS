subroutine finalize
!
! Finalize the computation
!
 use mod_streams
 implicit none
 logical :: pending_wallp, pending_slicexy
 logical :: opened_wallp,  opened_slicexy
 real(mykind) :: elapsed,startTiming,endTiming
!
 if (io_type == 1) then
  if (xrecyc>0._mykind) call recyc
  call updateghost()
  call prims()
  call copy_gpu_to_cpu()
  call mpi_barrier(mpi_comm_world, iermpi)
  startTiming = mpi_wtime()
  call writerst_serial()
  call mpi_barrier(mpi_comm_world, iermpi)
  endTiming = mpi_wtime()
  elapsed = endTiming-startTiming
  if (masterproc) write(error_unit,*) 'I/O serial restart time =', elapsed
  if (iflow==-1) then
  elseif (iflow==0) then
   call writestatchann()
   call writestat1d()
  else 
   call writestatbl()
   call writestat2d_serial()
   call writedf_serial()
  endif
 elseif (io_type == 2) then
  if (xrecyc>0._mykind) call recyc
  call updateghost()
  call prims()
  call copy_gpu_to_cpu()
  startTiming = mpi_wtime()
  call writerst()
  endTiming = mpi_wtime()
  elapsed = endTiming-startTiming
  if (masterproc) write(error_unit,*) 'I/O MPI restart time =', elapsed
  if (iflow==-1) then
  elseif (iflow==0) then
   call writestatchann()
   call writestat1d()
  else 
   call writestatbl()
   call writestat2d()
   call writedf()
  endif
 endif
!
 if (io_type /= 0) then
  if (enable_plot3d>0) call writefield()
  if (enable_vtk>0) call writefieldvtk()
 endif
!
end subroutine finalize
