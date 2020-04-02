subroutine writefieldvtk
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
 integer :: size_real
 integer :: ntot
 integer (KIND=MPI_OFFSET_KIND) :: offset
 integer (KIND=MPI_OFFSET_KIND) :: offset_x, offset_y, offset_z, delta_offset_w
!
 character(4) :: nastore
 character(128) :: filename
!
 integer, parameter :: int64_kind = selected_int_kind(2*range(1))
 integer(int64_kind) :: gridsize_64
 character(len=16) :: int2str
 character(len=32) :: int2str_o
 character(len=65536) :: xml_part
 character(len=7) :: vtk_float
!
 if (masterproc) print *, 'Storing VTK sol', istore,'at time', telaps
 write(nastore,1004) istore
 1004 format(I4.4)

 filename = 'field_'//nastore//'.vtr'

 call MPI_TYPE_SIZE(mpi_prec,size_real,iermpi)
 if(size_real == 4) then
  vtk_float = "Float32"
 elseif(size_real == 8) then
  vtk_float  = "Float64"
 else
  if(masterproc) write(*,*) "Error on VTK write! size_real must be either 4 or 8"
  call MPI_ABORT(MPI_COMM_WORLD, iermpi, iermpi)
 endif
 gridsize_64 = int(size_real,int64_kind)*int(nxmax,int64_kind)*int(nymax,int64_kind)*int(nzmax,int64_kind)

 if(storage_size(gridsize_64) /= 64) then
  if(masterproc) write(*,*) "Error on VTK write! Size of int64_kind integers is not 8 bytes!"
  call MPI_ABORT(MPI_COMM_WORLD, iermpi, iermpi)
 endif

 if (masterproc) then
  offset_x = 0 
  offset_y = size_real*nxmax + storage_size(gridsize_64)/8
  offset_z = offset_y + size_real*nymax + storage_size(gridsize_64)/8
  delta_offset_w = gridsize_64 + storage_size(gridsize_64)/8 ! the second part is because of the header of bytes before data

  open(unit=123, file=filename, access="stream", form="unformatted", status="replace")
   xml_part = ' &
   & <?xml version="1.0"?> &
   & <VTKFile type="RectilinearGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64"> &
   &  <RectilinearGrid WholeExtent="+1 +'//int2str(nxmax)//' +1 +'//int2str(nymax)//' +1 +'//int2str(nzmax)//'"> &
   &   <Piece Extent="+1 +'//int2str(nxmax)//' +1 +'//int2str(nymax)//' +1 +'//int2str(nzmax)//'"> &
   &    <Coordinates> &
   &     <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="X" format="appended" offset="'//int2str_o(offset_x)//'"/> &
   &     <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="Y" format="appended" offset="'//int2str_o(offset_y)//'"/> &
   &     <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="Z" format="appended" offset="'//int2str_o(offset_z)//'"/> &
   &    </Coordinates> &
   &   <PointData> '
   offset = offset_z + size_real*nzmax + storage_size(gridsize_64)/8
   do l=1,nv
    xml_part = trim(adjustl(xml_part)) // ' &
   & <DataArray type="'//vtk_float//'" NumberOfComponents="1" Name="w_'//int2str(l)//'" format="appended" &
   &  offset="'//int2str_o(offset)//'"/>'
    offset = offset + delta_offset_w
   enddo
   xml_part = trim(adjustl(xml_part)) // ' &
   &       </PointData> &
   &     </Piece> &
   &   </RectilinearGrid> &
   &   <AppendedData encoding="raw"> ' 
   write(123) trim(adjustl(xml_part))

   write(123) "_"
   !write(123) storage_size(gridsize_64)/8*int(nxmax,int64_kind) , xg(1:nxmax)
   !write(123) storage_size(gridsize_64)/8*int(nymax,int64_kind) , yg(1:nymax)
   !write(123) storage_size(gridsize_64)/8*int(nzmax,int64_kind) , zg(1:nzmax)
   write(123) size_real*int(nxmax,int64_kind) , xg(1:nxmax)
   write(123) size_real*int(nymax,int64_kind) , yg(1:nymax)
   write(123) size_real*int(nzmax,int64_kind) , zg(1:nzmax)
   flush(123)
  close(123)
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
 call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,mpi_prec,filetype,iermpi)
 call MPI_TYPE_COMMIT(filetype,iermpi)
!
 do l=1,nv
  call MPI_BARRIER(mp_cart, iermpi)
  if (masterproc) then
   open(unit=123, file=filename, access="stream", form="unformatted", position="append")
    write(123) gridsize_64
    flush(123)
   close(123)
  endif
  call MPI_BARRIER(mp_cart, iermpi)
  call MPI_FILE_OPEN(mp_cart,filename,MPI_MODE_RDWR,MPI_INFO_NULL,mpi_io_file,iermpi)
  call MPI_FILE_GET_SIZE(mpi_io_file, offset, iermpi)
  call MPI_BARRIER(mp_cart, iermpi)
  call MPI_FILE_SET_VIEW(mpi_io_file,offset,mpi_prec,filetype,"native",MPI_INFO_NULL,iermpi)
  call MPI_FILE_WRITE_ALL(mpi_io_file,w(l,1:nx,1:ny,1:nz),ntot,mpi_prec,istatus,iermpi)
  call MPI_FILE_CLOSE(mpi_io_file,iermpi)
 enddo
!
 call MPI_TYPE_FREE(filetype,iermpi)
!
 if (masterproc) then
  open(unit=123, file=filename, access="stream", position="append", form="unformatted")
   write(123) ' &
    &    </AppendedData> &
    &  </VTKFile>' 
  close(123)
 endif
end subroutine writefieldvtk
!
function int2str(int_num)
implicit none
integer :: int_num
character(len=16) :: int2str, ret_value
write(ret_value, "(I0)") int_num
int2str = ret_value
endfunction int2str

function int2str_o(int_num)
use mpi
implicit none
integer(KIND=MPI_OFFSET_KIND) :: int_num
character(len=32) :: int2str_o, ret_value
write(ret_value, "(I0)") int_num
int2str_o = ret_value
endfunction int2str_o

