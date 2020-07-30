subroutine readstat2d
!
! Reading file stat.bin with statistics
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
 integer :: ntotxy
 integer (kind=mpi_offset_kind) :: offset
 character(len=256) :: oldname, newname
!
 if (ncoords(3)==0) then
!
  sizes(1) = nblocks(1)*nx
  sizes(2) = nblocks(2)*ny
  sizes(3) = 1
  subsizes(1) = nx
  subsizes(2) = ny
  subsizes(3) = 1 
  starts(1) = 0 + ncoords(1)*subsizes(1)
  starts(2) = 0 + ncoords(2)*subsizes(2)
  starts(3) = 0
  ntotxy = nx*ny
!
  call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,iermpi)
  call mpi_type_commit(filetype,iermpi)
  call mpi_file_open(mp_cartx,'stat.bin',mpi_mode_rdonly,mpi_info_null,mpi_io_file,iermpi)
  offset = 0
  do l=1,nvmean
   call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
   call mpi_file_read_all(mpi_io_file,w_av(l,1:nx,1:ny),ntotxy,mpi_prec,istatus,iermpi)
   call mpi_type_size(mpi_prec,size_real,iermpi)
   do m=1,nblocks(1)*nblocks(2)
    offset = offset+size_real*ntotxy
   enddo
  enddo
  call mpi_file_close(mpi_io_file,iermpi)
  call mpi_type_free(filetype,iermpi)
!
 endif
!
 call mpi_bcast(w_av,nvmean*nx*ny,mpi_prec,0,mp_cartz,iermpi)
!
 if (masterproc) then
  oldname = c_char_"stat.bin"//c_null_char
  newname = c_char_"stat.bak"//c_null_char
  iermpi = rename_wrapper(oldname, newname)
  if (iermpi /= 0) write(error_unit,*) "Warning! Cannot rename file stat.bin to stat.bak"
 endif
 call mpi_barrier(mpi_comm_world, iermpi)
!
end subroutine readstat2d

subroutine readstat1d
!
! Reading file stat.bin with statistics
!
 use mod_streams
 use mod_sys
 implicit none
 character(len=256) :: oldname, newname
!
 if (masterproc) then
  open(12,file='stat.bin',form='unformatted')
  read(12) w_av_1d
  close(12)
 endif
!
 call mpi_bcast(w_av_1d,nvmean*ny,mpi_prec,0,mp_cart,iermpi)
!
 if (masterproc) then
  oldname = c_char_"stat.bin"//c_null_char
  newname = c_char_"stat.bak"//c_null_char
  iermpi = rename_wrapper(oldname, newname)
  if (iermpi /= 0) write(error_unit,*) "Warning! Cannot rename file stat.bin to stat.bak"
 endif
 call mpi_barrier(mpi_comm_world, iermpi)
!
end subroutine readstat1d
!
subroutine readstat2d_serial
!
! Reading file wavplot.dat
!
 use mod_streams
 implicit none
!
 integer :: i,j,l
 real(mykind) :: xx,yy
!
 if (ncoords(3)==0) then

  open(unit=12,file='wavplot0_'//chx//'_'//chy//'.dat',form='formatted')
  read(12,*)
  do j=1,ny
   do i=1,nx
    read(12,*) xx,yy,(w_av(l,i,j),l=1,nvmean)
   enddo
  enddo
  close(12)

! open(unit=12,file='wavplot_'//chx//'_'//chy//'.bin',form='unformatted')
! read(12) w_av
! close(12)

 endif
!
 call mpi_bcast(w_av,nvmean*nx*ny,mpi_prec,0,mp_cartz,iermpi)
!
end subroutine readstat2d_serial
