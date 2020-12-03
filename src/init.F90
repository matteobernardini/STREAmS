subroutine init
!
! Providing initial conditions
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,m
!
 dfupdated      = .true.
 ncyc0          = 0         ! number of cycles is set to zero
 istore         = 1         ! index of the first solution to store
 istore_restart = 1         ! index of the first solution to store restart
 telaps0        = 0._mykind ! time counter is set to zero
 stat_io        = 'rewind'  ! I/O status
!
! Freestream initialization (just to avoid 0 in double ghost nodes)
!
 do k=1-ng,nz+ng
  do j=1-ng,ny+ng
   do i=1-ng,nx+ng
    do m=1,nv
     w(m,i,j,k) = winf(m)
    enddo
   enddo
  enddo
 enddo
 ducros = .false.
 if (tresduc<=0._mykind) ducros = .true.
!
 if (idiski>=1) then
!
! Reading solution from file 
!
  if(io_type == 1) call readrst_serial
  if(io_type == 2) call readrst
  if (iflow==-1) then
  elseif (iflow==0) then
   if (idiski > 1) call readstat1d
   call generatewmean_channel
  else
   if(io_type == 1) then
    call readdf_serial
    if (idiski > 1) call readstat2d_serial
   elseif(io_type == 2) then
    call readdf
    if (idiski > 1) call readstat2d
   endif
   call generatewmean
   if (masterproc) call target_reystress
   call mpi_bcast(amat_df,9*ny,mpi_prec,0,mp_cart,iermpi)
  endif
!
  istore          = int(telaps0/dtsave) + 1
  istore_restart  = int(telaps0/dtsave_restart) + 1
  stat_io         = 'append'
!
  if (idiski==1) then
   w_av    = 0._mykind 
   w_av_1d = 0._mykind
   itav    = 0
  endif
!
 else ! start from scratch
!
! Setting the initial conditions
!
  w_av    = 0._mykind 
  w_av_1d = 0._mykind
  itav    = 0
!
  if (iflow==-1) then
   call init_windtunnel
  elseif (iflow==0) then
   call generatewmean_channel
   call init_channel
  else
   call generatewmean
   if (masterproc) call target_reystress
   call mpi_bcast(amat_df,9*ny,mpi_prec,0,mp_cart,iermpi)
   call initurb
  endif
!
 endif ! End of definition of the initial flowfield
!
end subroutine init
