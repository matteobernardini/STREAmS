subroutine writerst
!
! Writing rst.bin and finaltime.dat
!
 use mod_streams
 implicit none
!
 integer :: mpi_io_file
 integer :: filetype
 integer :: local_datatype
 integer, dimension(3) :: sizes     ! Dimensions of the total grid
 integer, dimension(3) :: subsizes  ! Dimensions of grid local to a procs
 integer, dimension(3) :: starts    ! Starting coordinates
 integer, dimension(3) :: memsizes
 integer :: size_real
 integer :: ntot,ntotxy,ntotyz
 integer (kind=mpi_offset_kind) :: offset
 integer :: i,j,k,l,m
 real(mykind) :: pp,rho,uu,vv,ww
!
 if (ncoords(3)==0) then
!
  open(unit=10,file='plotxy_'//chx//'.dat',form='formatted')
  write(10,*) 'zone i=',nx,', j=',ny
  k = nz/2+1
  do j=1,ny
   do i=1,nx
    rho = w(1,i,j,k)
    uu  = w(2,i,j,k)/rho
    vv  = w(3,i,j,k)/rho
    ww  = w(4,i,j,k)/rho
    pp  = rho*temperature(i,j,k)
    write(10,100) x(i),y(j),rho,uu,vv,ww,pp,temperature(i,j,k)
  100     format(20ES20.10)
    enddo
   enddo
   close(10)
 endif
!
 if (masterproc) write(*,*) 'Writing rst.bin'
!
 sizes(1) = nblocks(1)*nx
 sizes(2) = nblocks(2)*ny
 sizes(3) = nblocks(3)*nz
 subsizes(1) = nx
 subsizes(2) = ny
 subsizes(3) = nz
 starts(1) = 0 + ncoords(1)*subsizes(1)
 starts(2) = 0 + ncoords(2)*subsizes(2)
 starts(3) = 0 + ncoords(3)*subsizes(3)
 ntot = nx*ny*nz
!
 call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,iermpi)
 call mpi_type_commit(filetype,iermpi)
 call mpi_file_open(mp_cart,'rst.bin',mpi_mode_create+mpi_mode_wronly,mpi_info_null,mpi_io_file,iermpi)
 offset = 0
 do l=1,nv
  call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
  call mpi_file_write_all(mpi_io_file,w(l,1:nx,1:ny,1:nz),ntot,mpi_prec,istatus,iermpi)
  call mpi_type_size(mpi_prec,size_real,iermpi)
  do m=1,nblocks(1)*nblocks(2)*nblocks(3)
   offset = offset+size_real*ntot
  enddo
 enddo
!
 call mpi_file_close(mpi_io_file,iermpi)
 call mpi_type_free(filetype,iermpi)
!
 if (masterproc) then
  open(11,file='finaltime.dat')
  write(11,*) icyc,itav,telaps
  close(11)
 endif
!
end subroutine writerst

subroutine writerst_serial
!
! Writing rst.bin and finaltime.dat
!
 use mod_streams
 implicit none
!
 integer :: mpi_io_file
 integer :: filetype
 integer :: local_datatype
 integer, dimension(3) :: sizes     ! Dimensions of the total grid
 integer, dimension(3) :: subsizes  ! Dimensions of grid local to a procs
 integer, dimension(3) :: starts    ! Starting coordinates
 integer, dimension(3) :: memsizes
 integer :: size_real
 integer :: ntot,ntotxy,ntotyz
 integer (kind=mpi_offset_kind) :: offset
 integer :: i,j,k,l,m
 real(mykind) :: pp,rho,uu,vv,ww
!
 if (ncoords(3)==0) then
!
  open(unit=10,file='plotxy_'//chx//'.dat',form='formatted')
  write(10,*) 'zone i=',nx,', j=',ny
  k = nz/2+1
  do j=1,ny
   do i=1,nx
    rho = w(1,i,j,k)
    uu  = w(2,i,j,k)/rho
    vv  = w(3,i,j,k)/rho
    ww  = w(4,i,j,k)/rho
    pp  = rho*temperature(i,j,k)
    write(10,100) x(i),y(j),rho,uu,vv,ww,pp,temperature(i,j,k)
  100     format(20ES20.10)
    enddo
   enddo
   close(10)
 endif
!
 if (masterproc) write(*,*) 'Writing rst.bin'
!
 fl = w_order(1:nx,1:ny,1:nz,1:nv)
!
 open (11,file='rst1_'//chx//'_'//chz//'.bin',form='unformatted')
 !rewind (11)
 write(11) fl
 !write   (11) ((((w(m,i,j,k),m=1,nv), &
 !                           i=1,nx),     &
 !                           j=1,ny),     &
 !                           k=1,nz)
 close  (11)
!
 if (masterproc) then
  open(11,file='finaltime1.dat')
  write(11,*) icyc,itav,telaps
  close(11)
 endif
!
end subroutine writerst_serial
