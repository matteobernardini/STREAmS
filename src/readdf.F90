subroutine readdf
!
! Reading file df.bin for digital filtering
!
 use mod_streams
 use mod_sys
 implicit none
!
 integer :: i,l,m
 integer :: mpi_io_file
 integer :: filetype
 integer :: local_datatype
 integer, dimension(3) :: sizes     ! Dimensions of the total grid
 integer, dimension(3) :: subsizes  ! Dimensions of grid local to a procs
 integer, dimension(3) :: starts    ! Starting coordinates
 integer, dimension(3) :: memsizes
 integer :: size_real
 integer :: ntotyz
 integer (kind=mpi_offset_kind) :: offset
 character(len=256) :: oldname, newname
 logical :: file_exists
!
 inquire(file="df.bin",exist=file_exists)
!
 if (ncoords(1)==0) then
!
  sizes(1) = 1
  sizes(2) = nblocks(2)*ny
  sizes(3) = nblocks(3)*nz
  subsizes(1) = 1
  subsizes(2) = ny
  subsizes(3) = nz
  starts(1) = 0
  starts(2) = 0 + ncoords(2)*subsizes(2)
  starts(3) = 0 + ncoords(3)*subsizes(3)
  ntotyz = ny*nz
!
  call mpi_type_create_subarray(3,sizes,subsizes,starts,mpi_order_fortran,mpi_prec,filetype,iermpi)
  call mpi_type_commit(filetype,iermpi)
  if (file_exists) then
   call mpi_file_open(mp_cartz,'df.bin',mpi_mode_rdonly,mpi_info_null,mpi_io_file,iermpi)
  else
   call mpi_file_open(mp_cartz,'df.bak',mpi_mode_rdonly,mpi_info_null,mpi_io_file,iermpi)
  endif
  offset = 0
  do i=1-ng,0
   do l=1,nv
    call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
    call mpi_file_read_all(mpi_io_file,w(l,i,1:ny,1:nz),ntotyz,mpi_prec,istatus,iermpi)
    call mpi_type_size(mpi_prec,size_real,iermpi)
    do m=1,nblocks(2)*nblocks(3)
     offset = offset+size_real*ntotyz
    enddo
   enddo
  enddo
  do l=1,3
   call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
   call mpi_file_read_all(mpi_io_file,vf_df(l,1:ny,1:nz),ntotyz,mpi_prec,istatus,iermpi)
   call mpi_type_size(mpi_prec,size_real,iermpi)
   do m=1,nblocks(2)*nblocks(3)
    offset = offset+size_real*ntotyz
   enddo
  enddo
  call mpi_file_close(mpi_io_file,iermpi)
  call mpi_type_free(filetype,iermpi)
!
 endif
!
 call mpi_barrier(mpi_comm_world, iermpi)
!
 if (file_exists) then
  call mpi_barrier(mpi_comm_world, iermpi)
  if (masterproc) then
   oldname = c_char_"df.bin"//c_null_char
   newname = c_char_"df.bak"//c_null_char
   iermpi = rename_wrapper(oldname, newname)
   if (iermpi /= 0) write(error_unit,*) "Warning! Cannot rename file df.bin to df.bak"
  endif
  call mpi_barrier(mpi_comm_world, iermpi)
 endif
!
end subroutine readdf

!
subroutine readdf_serial
!
! Reading df.bin
!
 use mod_streams
 implicit none
!
 integer :: i,l,m,j,k
 integer :: mpi_io_file
 integer :: filetype
 integer :: local_datatype
 integer, dimension(3) :: sizes     ! Dimensions of the total grid
 integer, dimension(3) :: subsizes  ! Dimensions of grid local to a procs
 integer, dimension(3) :: starts    ! Starting coordinates
 integer, dimension(3) :: memsizes
 integer :: size_real
 integer :: ntotyz
 integer (KIND=MPI_OFFSET_KIND) :: offset
!
 if (ncoords(1)==0) then
  open(12,file='vf0_df_'//chz//'.bin',form='unformatted')
  do k=1,nz
   do j=1,ny
    read(12) (vf_df(m,j,k),m=1,3)
   enddo
  enddo
  read(12) w(:,1-ng:0,:,:)
  close(12)
 endif
!
endsubroutine readdf_serial
