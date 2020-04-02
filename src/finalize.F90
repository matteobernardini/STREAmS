subroutine finalize
!
! Finalize the computation
!
 use mod_streams
 implicit none
!
 call updateghost()
 call prims()
 call copy_gpu_to_cpu()
 call writerst()
 if (iflow==-1) then
 elseif (iflow==0) then
  call writestatchann()
  call writestat1d()
 else 
  call writestatbl()
  call writestat2d()
  call writedf()
 endif
!
end subroutine finalize
