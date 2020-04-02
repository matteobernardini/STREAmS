subroutine writegridplot3d
!
! Writing grid in plot3D format
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,l,m
!
 integer :: mpi_io_file
 integer :: filetype
 integer :: local_datatype
 integer,dimension(3) :: sizes     ! Dimensions of the total grid
 integer,dimension(3) :: subsizes  ! Dimensions of grid local to a procs
 integer,dimension(3) :: starts    ! Starting coordinates
 integer,dimension(3) :: memsizes
 integer :: size_real, size_integer
 integer :: ntot
 integer (kind=mpi_offset_kind) :: offset
 real(mykind), dimension(3,nx,ny,nz) :: grid3d
!
 character(4) :: nastore
!
 if (masterproc) write(*,*) 'Storing grid in plot3D format'
!
 if (masterproc) then
  open(unit=123, file='plot3dgrid.xyz', access="stream", form="unformatted")
   write(123) nxmax,nymax,nzmax
   flush(123)
  close(123)
 endif
 call mpi_barrier(mpi_comm_world,iermpi)
!
 do k=1,nz
  do j=1,ny
   do i=1,nx
     grid3d(1,i,j,k) = x(i)
     grid3d(2,i,j,k) = y(j)
     grid3d(3,i,j,k) = z(k)
   enddo
  enddo
 enddo
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
 call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,iermpi)
 call MPI_TYPE_COMMIT(filetype,iermpi)
 call MPI_FILE_OPEN(mp_cart,'plot3dgrid.xyz',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,iermpi)
 size_integer = storage_size(nx)/8
 offset = 3*size_integer
 do l=1,3
 call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
 call MPI_FILE_WRITE_ALL(mpi_io_file,grid3d(l,1:nx,1:ny,1:nz),ntot,mpi_prec,istatus,iermpi)
 call MPI_TYPE_SIZE(mpi_prec,size_real,iermpi)
  do m=1,nblocks(1)*nblocks(2)*nblocks(3)
   offset = offset+size_real*ntot
  enddo
 enddo
! 
 call MPI_FILE_CLOSE(mpi_io_file,iermpi)
 call MPI_TYPE_FREE(filetype,iermpi)
!
end subroutine writegridplot3d
