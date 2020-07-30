subroutine writedf
!
! Writing df.bin
!
 use mod_streams
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
 integer (KIND=MPI_OFFSET_KIND) :: offset
!
 if (masterproc) write(*,*) 'Writing df.bin'
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
  call mpi_file_open(mp_cartz,'df.bin',mpi_mode_create+mpi_mode_wronly,mpi_info_null,mpi_io_file,iermpi)
  offset = 0
  do i=1-ng,0
   do l=1,nv
    call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
    call mpi_file_write_all(mpi_io_file,w(l,i,1:ny,1:nz),ntotyz,mpi_prec,istatus,iermpi)
    call mpi_type_size(mpi_prec,size_real,iermpi)
    do m=1,nblocks(2)*nblocks(3)
     offset = offset+size_real*ntotyz
    enddo
   enddo
  enddo
  do l=1,3
   call mpi_file_set_view(mpi_io_file,offset,mpi_prec,filetype,"native",mpi_info_null,iermpi)
   call mpi_file_write_all(mpi_io_file,vf_df(l,1:ny,1:nz),ntotyz,mpi_prec,istatus,iermpi)
   call mpi_type_size(mpi_prec,size_real,iermpi)
   do m=1,nblocks(2)*nblocks(3)
    offset = offset+size_real*ntotyz
   enddo
  enddo
!
  call mpi_file_close(mpi_io_file,iermpi)
  call mpi_type_free(filetype,iermpi)
!
 endif
!
end subroutine writedf
!
subroutine writedf_serial
!
! Writing df.bin
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
  open(12,file='vf1_df_'//chz//'.bin',form='unformatted')
  do k=1,nz
   do j=1,ny
    write(12) (vf_df(m,j,k),m=1,3)
   enddo
  enddo
  write(12) w(:,1-ng:0,:,:)
  close(12)
 endif
!
endsubroutine writedf_serial
