subroutine readrst
!
! Reading file rst.bin and finaltime.dat for restart
!
 use mod_streams
 use mod_sys
 implicit none
!
 integer :: l,m
 integer :: mpi_io_file
 integer :: filetype
 integer :: local_datatype
 integer, dimension(3) :: sizes     ! Dimensions of the total grid
 integer, dimension(3) :: subsizes  ! Dimensions of grid local to a procs
 integer, dimension(3) :: starts    ! Starting coordinates
 integer, dimension(3) :: memsizes
 integer :: size_real
 integer :: ntot
 integer (kind=mpi_offset_kind) :: offset
 character(len=256) :: oldname, newname
 logical :: file1_exists, file2_exists
!
 if (masterproc) write(*,*) 'Reading restart files'
!
 inquire(file="finaltime.dat",exist=file1_exists)
 if (file1_exists) then
  open(11,file='finaltime.dat')
 else
  open(11,file='finaltime.bak')
 endif
 read(11,*) ncyc0,itav,telaps0
 close(11)
!
 if (file1_exists) then
  call mpi_barrier(mpi_comm_world, iermpi)
  if (masterproc) then
   oldname = c_char_"finaltime.dat"//c_null_char
   newname = c_char_"finaltime.bak"//c_null_char
   iermpi = rename_wrapper(oldname, newname)
   if (iermpi /= 0) write(error_unit,*) "Warning! Cannot rename file finaltime.dat to finaltime.bak"
  endif
  call mpi_barrier(mpi_comm_world, iermpi)
 endif
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
 inquire(file="rst.bin",exist=file2_exists)
 if (file2_exists) then
  call mpi_file_open(mp_cart,'rst.bin',mpi_mode_rdonly,mpi_info_null,mpi_io_file,iermpi)
 else
  call mpi_file_open(mp_cart,'rst.bak',mpi_mode_rdonly,mpi_info_null,mpi_io_file,iermpi)
 endif
 offset = 0
 do l=1,nv
  call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
  call mpi_file_read_all(mpi_io_file,w(l,1:nx,1:ny,1:nz),ntot,mpi_prec,istatus,iermpi)
  call mpi_type_size(mpi_prec,size_real,iermpi)
  do m=1,nblocks(1)*nblocks(2)*nblocks(3)
   offset = offset+size_real*ntot
  enddo
 enddo
 call mpi_file_close(mpi_io_file,iermpi)
 call mpi_type_free(filetype,iermpi)
 call mpi_barrier(mpi_comm_world, iermpi)
!
 if (file2_exists) then
  call mpi_barrier(mpi_comm_world, iermpi)
  if (masterproc) then
   oldname = c_char_"rst.bin"//c_null_char
   newname = c_char_"rst.bak"//c_null_char
   iermpi = rename_wrapper(oldname, newname)
   if (iermpi /= 0) write(error_unit,*) "Warning! Cannot rename file rst.bin to rst.bak"
  endif
  call mpi_barrier(mpi_comm_world, iermpi)
 endif
!
end subroutine readrst

subroutine readrst_serial
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
 integer :: i,j,k,l,m,iv
 real(mykind) :: pp,rho,uu,vv,ww
!
 if (masterproc) write(*,*) 'Reading rst.bin'
 write(*,*) 'START - Reading rst.bin from nrank = ',nrank
!
 open (11,file='rst0_'//chx//'_'//chz//'.bin',form='unformatted')
 !rewind (11)
 read   (11) fl
 !read   (11) ((((w(m,i,j,k),m=1,nv), &
 !                           i=1,nx),     &
 !                           j=1,ny),     &
 !                           k=1,nz)
 close  (11)
 write(*,*) 'END - Reading rst.bin from nrank = ',nrank
!
 do iv=1,5
  do k=1,nz
   do j=1,ny
    do i=1,nx
     w(iv,i,j,k) = fl(i,j,k,iv)
    enddo
   enddo
  enddo
 enddo
!
 write(*,*) 'COPY - Reading rst.bin from nrank = ',nrank
!
 write(*,*) 'STARTTIME - Reading finaltime0 from nrank = ',nrank
 open(11,file='finaltime0.dat')
 ! read(11,*) icyc,itav,telaps
 read(11,*) ncyc0,itav,telaps0
 close(11)
 write(*,*) 'ENDTIME - Reading finaltime0 from nrank = ',nrank
!
end subroutine readrst_serial
