subroutine writefield
!
! Writing field (MPI I/O)
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
 integer :: size_real,size_integer
 integer :: ntot
 integer (kind=mpi_offset_kind) :: offset
!
 character(4) :: nastore
!
 size_real = storage_size(re)/8
 size_integer = storage_size(nx)/8
!
 if (masterproc) write(*,*) 'Storing sol', istore,'at time', telaps
 write(nastore,1004) istore
 1004 format(I4.4)
!
 if (masterproc) then
  open(unit=123, file='field_'//nastore//'.q', access="stream", form="unformatted")
   write(123) nxmax,nymax,nzmax
   write(123) rm, 0._mykind, re, telaps
   flush(123)
  close(123)
 endif
 call mpi_barrier(mp_cart,iermpi)
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
 call MPI_FILE_OPEN(mp_cart,'field_'//nastore//'.q',MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,iermpi)
 offset = 3*size_integer+4*size_real
 do l=1,nv
  call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
  select case (l)
  case(1)
   call MPI_FILE_WRITE_ALL(mpi_io_file,w(l,1:nx,1:ny,1:nz),ntot,mpi_prec,istatus,iermpi)
  case(2:4)
   call MPI_FILE_WRITE_ALL(mpi_io_file,w(l,1:nx,1:ny,1:nz)/w(1,1:nx,1:ny,1:nz),ntot,mpi_prec,istatus,iermpi)
  case(5)
   call MPI_FILE_WRITE_ALL(mpi_io_file,w(1,1:nx,1:ny,1:nz)*temperature(1:nx,1:ny,1:nz),ntot,mpi_prec,istatus,iermpi)
  end select
! if (l<nv) then
!  call MPI_FILE_WRITE_ALL(mpi_io_file,w(l,1:nx,1:ny,1:nz),ntot,mpi_prec,istatus,iermpi)
! else
!  call MPI_FILE_WRITE_ALL(mpi_io_file,w(1,1:nx,1:ny,1:nz)*temperature(1:nx,1:ny,1:nz),ntot,mpi_prec,istatus,iermpi)
! endif
  do m=1,nblocks(1)*nblocks(2)*nblocks(3)
   offset = offset+size_real*ntot
  enddo
 enddo
!
 call MPI_FILE_CLOSE(mpi_io_file,iermpi)
 call MPI_TYPE_FREE(filetype,iermpi)
!
end subroutine writefield
