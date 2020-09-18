subroutine startmpi
!
! Initialize MPI (and CUDA) environment
!
 use mod_streams
 implicit none
!
 logical :: reord
 logical :: remain_dims(ndims)
 integer :: i_skip, dims(2)
 !integer(kind=cuda_count_kind) :: req_stacksize_gpu
!
 call mpi_init(iermpi)
 call mpi_comm_rank(mpi_comm_world,nrank,iermpi)
 call mpi_comm_size(mpi_comm_world,nproc,iermpi)

#ifdef USE_CUDA
! 
!  Assign a different GPU to each MPI rank
!  Note: all the memory allocation should be dynamic, otherwise all the arrays will be allocated on device 0
!
 mydev=0
 call mpi_comm_split_type(mpi_comm_world,mpi_comm_type_shared,0,mpi_info_null,local_comm,iermpi)
 call mpi_comm_rank(local_comm,mydev,iermpi)
 iermpi = cudaSetDevice(mydev)
 write(*,*) "MPI rank",nrank,"using GPU",mydev
 iermpi = cudaStreamCreate(stream1)
 iermpi = cudaStreamCreate(stream2)
 !req_stacksize_gpu = 32000000
 !iermpi = cudaDeviceSetLimit( cudaLimitStackSize, req_stacksize_gpu )
 !write(*,*) "GPU stack enlarged"
 !possible limit arguments are cudaLimitStackSize, cudaLimitPrintfSize, and cudaLimitMallocHeapSize.
#endif
!
 allocate(ncoords(ndims))
 allocate(nblocks(ndims))
 allocate(pbc(ndims))
!
 call check_input(0)
!
 open (unit=12,file='input.dat',form='formatted')
 do i_skip = 1,15
  read (12,*)
 enddo
 read (12,*)
 read (12,*) iflow
 read (12,*)
 read (12,*)
 read (12,*) rlx,rly,rlz
 read (12,*)
 read (12,*)
 read (12,*) nxmax,nymax,nzmax
 read (12,*)
 read (12,*)
 read (12,*) nymaxwr,rlywr,dyp_target
 read (12,*)
 read (12,*)
 read (12,*) ng, ivis, iorder, iweno
 read (12,*)
 read (12,*)
 read (12,*) nblocks(1),nblocks(3)
 close(12)
!
 ngdf = iorder/2
!
 masterproc = .false.
 if (nrank==0) masterproc = .true.
!
 nblocks(2) = 1
!
 if(nblocks(1) <=0 .or. nblocks(3) <=0) then
     dims = [0,0]
     call MPI_Dims_create(nproc, 2, dims, iermpi)
     nblocks(1) = dims(1)
     nblocks(3) = dims(2)
     if (masterproc) write(error_unit,'(A,I0,A,I0)') 'Automatic MPI decomposition:', nblocks(1),' x ', nblocks(3)
 endif
!
 nx = nxmax/nblocks(1)
 ny = nymax/nblocks(2)
 nz = nzmax/nblocks(3)
!
 call check_input(1)
!
 pbc(1) = .false.
 pbc(2) = .false.
 pbc(3) = .true.
 if (iflow==-1) then ! Wind tunnel
  pbc(3) = .false.
 endif
 if (iflow==0) then ! Channel flow (periodicity also applied in x)
  pbc(1) = .true.
 endif
!
! Create 3D topology
!
 reord = .false.
 call mpi_cart_create(mpi_comm_world,ndims,nblocks,pbc,reord,mp_cart,iermpi)
 call mpi_cart_coords(mp_cart,nrank,ndims,ncoords,iermpi)
!
! Create 1D communicators
!
 remain_dims(1) = .true.
 remain_dims(2) = .false.
 remain_dims(3) = .false.
 call mpi_cart_sub(mp_cart,remain_dims,mp_cartx,iermpi)
 call mpi_comm_rank(mp_cartx,nrank_x,iermpi)
 call mpi_cart_shift(mp_cartx,0,1,ileftx,irightx,iermpi)
 remain_dims(2) = .true.
 remain_dims(1) = .false.
 remain_dims(3) = .false.
 call mpi_cart_sub(mp_cart,remain_dims,mp_carty,iermpi)
 call mpi_comm_rank(mp_carty,nrank_y,iermpi)
 call mpi_cart_shift(mp_carty,0,1,ilefty,irighty,iermpi)
 remain_dims(3) = .true.
 remain_dims(1) = .false.
 remain_dims(2) = .false.
 call mpi_cart_sub(mp_cart,remain_dims,mp_cartz,iermpi)
 call mpi_comm_rank(mp_cartz,nrank_z,iermpi)
 call mpi_cart_shift(mp_cartz,0,1,ileftz,irightz,iermpi)
!
end subroutine startmpi
