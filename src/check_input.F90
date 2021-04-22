subroutine check_input(step)
!
! Check input file to detect possible errors and warnings
!
 use mod_streams
 implicit none
!
 integer, intent(in) :: step
 logical :: file_exists, file2_exists
#ifdef USE_CUDA
 real(mykind) :: correction_factor, gpu_used_mem, gpu_total_mem_real
 integer(kind=cuda_count_kind) :: gpu_total_mem, gpu_free_mem
#endif
!
 if (step == 0) then
  inquire(file="input.dat", exist=file_exists)
  if (.not.file_exists) then
   call fail_input("input.dat file not found")
  endif
 endif
!
 if (step == 1) then
  if (iflow < -1 .or. iflow > 2) then
   call fail_input("iflow must be -1, 0, 1, or 2")
  endif 
  if (mod(iorder,2)/=0) then
   call fail_input("iorder must be even")
  endif
  if (mod(ivis,2)/=0) then
   call fail_input("ivis must be even")
  endif
  if (iorder>2*ng) then
   call fail_input("iorder cannot be greater than 2*ng")
  endif
  if (iorder>8) then
   call fail_input("iorder cannot be greater than 8")
  endif
  if (ivis>iorder) then
   call fail_input("ivis cannot be greater than iorder")
  endif
  if (iweno>iorder/2) then
   call fail_input("iweno cannot be greater than iorder/2")
  endif
  if (nblocks(1) <= 0 .or. nblocks(3) <= 0) then
   call fail_input("nblocks must be integer numbers greater than 0")
  endif
  if (nblocks(1)*nblocks(3) /= nproc) then
   call fail_input("nblocks(1)*nblocks(3) must be equal to number of MPI procs")
  endif
  if (mod(nxmax, nblocks(1)) /= 0) then
   call fail_input("nxmax must be multiple of nblocks(1)")
  endif
  if (mod(nzmax, nblocks(3)) /= 0) then
   call fail_input("nzmax must be multiple of nblocks(3)")
  endif
  if (iflow==0.and.jbgrid>0) then
   if (mod(nymax,2)/=0) call fail_input("nymax must be even if jbgrid > 0")
  endif
  if (iflow > 0) then
   inquire(file="database_bl.dat", exist=file_exists)
   if (.not.file_exists) then
    call fail_input("database_bl.dat file not found")
   endif
   if (mod(nxmax, nblocks(3)) /= 0) then
    call fail_input("if iflow > 0, nxmax must be multiple of nblocks(3)")
   endif
   if (mod(nymax, nblocks(1)) /= 0) then
    call fail_input("if iflow > 0, nymax must be multiple of nblocks(1)")
   endif
   if (rly <= rlywr) then
    call fail_input("if iflow > 0, rly must be greater than rlywr")
   endif
   if (nymax <= nymaxwr) then
    call fail_input("if iflow > 0, nymax must be greater than nymaxwr")
   endif
  endif 
#ifdef USE_CUDA
  if(masterproc) then
   write(error_unit,*) 'Checking GPU memory allocation'
   write(error_unit,'(A,(3(I0,2x)))') ' GPU accelerated version, grid size:',nxmax, nymax, nzmax
   write(error_unit,'(A,(3(I0,2x)))') ' GPU accelerated version, process grid size:',nx, ny, nz
   gpu_used_mem = 43._mykind      ! Number of 3D arrays on GPU
   correction_factor = 1.5_mykind ! Safety margin
   gpu_used_mem = gpu_used_mem+correction_factor 
   gpu_used_mem = gpu_used_mem*real((nx+2*ng),mykind)*real((ny+2*ng),mykind)*real((nz+2*ng),mykind)
   gpu_used_mem = gpu_used_mem*storage_size(1._mykind)/8._mykind/(1024._mykind**2)
   write(error_unit,*) 'Estimated Memory usage (MB) per MPI task: ',gpu_used_mem
   iercuda = cudaMemGetInfo(gpu_free_mem, gpu_total_mem)
   gpu_total_mem_real = real(gpu_total_mem,mykind)/(1024._mykind**2)
   write(error_unit, *) " Memory on GPU (MB): ", gpu_total_mem_real
   if (gpu_used_mem > 1.05_mykind*gpu_total_mem_real) then
    call fail_input("Requested memory greatly exceeds the total available. The &
     & simulation will surely crash. Aborting now!!!")
   elseif (gpu_used_mem > gpu_total_mem_real) then
    call warning_input("Requested memory exceeds the total available. The &
     & simulation will most probably run out of memory and crash!!!")
   endif
   print *," "
  endif
#endif
 endif
!
 if (step == 2) then
  if (iflow==0.or.iflow==1) then
   if (tresduc < 1._mykind) then
    call warning_input("iflow == 0 or 1 and tresduc < 1: you may consider to avoid &
     & the activation of shock capturing scheme by setting tresduc >= 1")
   endif
  endif
  if (idiski <= 1 .and. ncyc < istat ) then
   call warning_input("forcing istat=ncyc to compute final meaningful statistics")
   istat = ncyc
  endif
  if(io_type == 2) then
   if (idiski >= 1) then
    inquire(file="rst.bin", exist=file_exists)
    if (.not.file_exists) then
     inquire(file="rst.bak", exist=file2_exists)
     if (.not.file2_exists) then
      call fail_input("rst.bin or rst.bak file not found")
     endif
    endif
    inquire(file="finaltime.dat", exist=file_exists)
    if (.not.file_exists) then
     inquire(file="finaltime.bak", exist=file2_exists)
     if (.not.file2_exists) then
      call fail_input("finaltime.dat or finaltime.bak file not found")
     endif
    endif
    if (iflow > 0) then
     inquire(file="df.bin", exist=file_exists)
     if (.not.file_exists) then
      inquire(file="df.bak", exist=file2_exists)
      if (.not.file2_exists) then
       call fail_input("df.bin or df.bak file not found")
      endif
     endif
    endif
    if (idiski > 1) then
     inquire(file="stat.bin", exist=file_exists)
     if (.not.file_exists) then
      call fail_input("stat.bin file not found")
     endif
    endif
   endif
  endif
  if (iflow == 2) then
   if (ximp < 0._mykind .or. ximp > rlx) then
    call fail_input("ximp must be in the range [0:rlx]")
   endif
  endif
 endif
!
 if (step == 3) then
  if (iflow == 2) then
   if (xsh < 0._mykind .or. xsh > rlx) then
    call fail_input("xsh must be in the range [0:rlx]. Please fix ximp to ensure that")
   endif
   if (xrecyc>0._mykind) then
    if (xrecyc>xsh) call warning_input("Recycling station located downstream the &
     intersection between the oblique shock and the top domain")
   endif
  endif
 endif
!
endsubroutine check_input
!
subroutine fail_input(msg)
 use mod_streams, only : masterproc, mpi_comm_world, iermpi
 use, intrinsic :: iso_fortran_env, only : error_unit
 implicit none
 character(len=*) :: msg
 if (masterproc) then
  write(error_unit,*) "Input Error! ", msg
  call mpi_abort(mpi_comm_world, 15, iermpi)
 endif
endsubroutine fail_input

subroutine warning_input(msg)
 use mod_streams, only : masterproc, mpi_comm_world, iermpi
 use, intrinsic :: iso_fortran_env, only : error_unit
 implicit none
 character(len=*) :: msg
 if (masterproc) then
  write(error_unit,*) "Input Warning! ", msg
 endif
endsubroutine warning_input
