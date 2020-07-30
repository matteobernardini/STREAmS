subroutine solver
!
! Solve the NS equations
!
 use mod_streams
 implicit none
!
 real(mykind) :: elapsed,startTiming,endTiming
 integer :: i
 logical :: stop_streams
!
 icyc = ncyc0
 telaps = telaps0
!
 call copy_cpu_to_gpu()
!
 if (masterproc) write(*,*) 'Compute time step'
 if (cfl>0._mykind) then
  call step()
 else
  dtmin = abs(cfl)
 endif
 if (masterproc) write(*,*) 'Done'
 call updateghost()
 call prims()
 if (tresduc<1._mykind) then
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
  !if(io_type > 0) then
  ! call write_wallpressure
  !endif
 !
  if (io_type>0) then
!
   if (mod(icyc,istat)==0) then
    call updateghost()
    call prims()
    call copy_gpu_to_cpu()
    if (iflow==-1) then
    elseif (iflow==0) then
     call stats1d()
    else
     call stats2d()
    endif
    call reset_cpu_gpu()
   endif
!
   if (telaps>tsol(istore)) then
    call updateghost()
    call prims()
    call copy_gpu_to_cpu()
    if(enable_plot3d > 0) call writefield()
    if(enable_vtk > 0) call writefieldvtk()
    if (iflow>0) call writestatzbl()
    istore = istore+1
    call reset_cpu_gpu()
   endif
!
  endif
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
  if (io_type==1) then
!
   if (telaps>tsol_restart(istore_restart)) then
    call updateghost()
    call prims()
    call copy_gpu_to_cpu()
    call writerst_serial()
    if (iflow==-1) then
    elseif (iflow==0) then
     call writestat1d()
    else
     call writestat2d_serial()
     call writedf_serial()
    endif
    istore_restart = istore_restart+1
    call reset_cpu_gpu()
   endif
!
  elseif (io_type == 2) then
!
   if (telaps>tsol_restart(istore_restart)) then
    call updateghost()
    call prims()
    call copy_gpu_to_cpu()
    call writerst()
    if (iflow==-1) then
    elseif (iflow==0) then
     call writestat1d()
    else
     call writestat2d()
     call writedf()
    endif
    istore_restart = istore_restart+1
    call reset_cpu_gpu()
   endif
  endif
!
  inquire(file="stop.stop",exist=stop_streams)
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
