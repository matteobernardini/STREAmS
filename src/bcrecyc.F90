subroutine bcrecyc(ilat)
!
! Apply recycling-rescaling boundary condition
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,l,ilat,m,ntot
 integer :: j_inn,j_out
 real(mykind) :: weta,weta1,bdamp,disty_inn,disty_out
 real(mykind) :: rhof_inn,rhof_out,uf_inn,uf_out,vf_inn,vf_out,wf_inn,wf_out
 real(mykind) :: rho,rhou,rhov,rhow,uu,vv,ww,qq,pp,tt
 real(mykind) :: ufav,vfav,wfav,rhom,rhofluc,ufluc,vfluc,wfluc
 real(mykind) :: rhomean,uumean,vvmean,wwmean,tmean
!
 if (ilat/=1) call fail_input("Recycling implemented only for ilat = 1")
!
! Compute spanwise averages at the recycling station
 !$cuf kernel do(2) <<<*,*>>>
 do j=1,ny
  do i=1,ng
   do m=1,nv
    wrecycav_gpu(i,j,m) = 0._mykind
    do k=1,nz
     wrecycav_gpu(i,j,m) = wrecycav_gpu(i,j,m)+wrecyc_gpu(i,j,k,m)
    enddo
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
 ntot = ng*ny*nv
 call mpi_allreduce(MPI_IN_PLACE,wrecycav_gpu,ntot,mpi_prec,mpi_sum,mp_cartz,iermpi)
! 
! Remove average
 !$cuf kernel do(2) <<<*,*>>>
 do j=1,ny
  do i=1,ng
   ufav = wrecycav_gpu(i,j,2)/wrecycav_gpu(i,j,1)
   vfav = wrecycav_gpu(i,j,3)/wrecycav_gpu(i,j,1)
   wfav = wrecycav_gpu(i,j,4)/wrecycav_gpu(i,j,1)
   rhom = wrecycav_gpu(i,j,1)/nzmax
   do k=1,nz
    wrecyc_gpu(i,j,k,2) = wrecyc_gpu(i,j,k,2)/wrecyc_gpu(i,j,k,1)-ufav ! Velocity fluctuations
    wrecyc_gpu(i,j,k,3) = wrecyc_gpu(i,j,k,3)/wrecyc_gpu(i,j,k,1)-vfav
    wrecyc_gpu(i,j,k,4) = wrecyc_gpu(i,j,k,4)/wrecyc_gpu(i,j,k,1)-wfav
    wrecyc_gpu(i,j,k,1) = wrecyc_gpu(i,j,k,1)-rhom                     ! Density fluctuations
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
!
!$cuf kernel do(2) <<<*,*>>>
 do k=1,nz
  do j=1,ny
! 
   weta  = weta_inflow_gpu(j)
   weta1 = 1._mykind-weta
   bdamp = 0.5_mykind*(1._mykind-tanh(4._mykind*(eta_inflow_gpu(j)-2._mykind))) 
   j_inn = map_j_inn_gpu(j) 
   j_out = map_j_out_gpu(j) 
   disty_inn = (yplus_inflow_gpu(j)-yplus_recyc_gpu(j_inn))/(yplus_recyc_gpu(j_inn+1)-yplus_recyc_gpu(j_inn))    
   disty_out = (eta_inflow_gpu(j)-eta_recyc_gpu(j_out))/(eta_recyc_gpu(j_out+1)-eta_recyc_gpu(j_out))    
!   
   do i=1,ng
!
    if (j==1.or.j_inn>=ny.or.j_out>=ny) then
     rhofluc = 0._mykind
     ufluc   = 0._mykind
     vfluc   = 0._mykind
     wfluc   = 0._mykind
    else
     rhof_inn = wrecyc_gpu(i,j_inn,k,1)*(1._mykind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,1)*disty_inn
     rhof_out = wrecyc_gpu(i,j_out,k,1)*(1._mykind-disty_out)+wrecyc_gpu(i,j_out+1,k,1)*disty_out
     uf_inn   = wrecyc_gpu(i,j_inn,k,2)*(1._mykind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,2)*disty_inn
     uf_out   = wrecyc_gpu(i,j_out,k,2)*(1._mykind-disty_out)+wrecyc_gpu(i,j_out+1,k,2)*disty_out
     vf_inn   = wrecyc_gpu(i,j_inn,k,3)*(1._mykind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,3)*disty_inn
     vf_out   = wrecyc_gpu(i,j_out,k,3)*(1._mykind-disty_out)+wrecyc_gpu(i,j_out+1,k,3)*disty_out
     wf_inn   = wrecyc_gpu(i,j_inn,k,4)*(1._mykind-disty_inn)+wrecyc_gpu(i,j_inn+1,k,4)*disty_inn
     wf_out   = wrecyc_gpu(i,j_out,k,4)*(1._mykind-disty_out)+wrecyc_gpu(i,j_out+1,k,4)*disty_out
!
     rhofluc = rhof_inn*weta1+rhof_out*weta
     ufluc   =   uf_inn*weta1+  uf_out*weta
     vfluc   =   vf_inn*weta1+  vf_out*weta
     wfluc   =   wf_inn*weta1+  wf_out*weta
     rhofluc = rhofluc*bdamp
     ufluc   = ufluc  *bdamp*betarecyc 
     vfluc   = vfluc  *bdamp*betarecyc 
     wfluc   = wfluc  *bdamp*betarecyc 
    endif
!
!   rhofluc = wrecyc_gpu(i,j,k,1)
!   ufluc   = wrecyc_gpu(i,j,k,2)
!   vfluc   = wrecyc_gpu(i,j,k,3)
!   wfluc   = wrecyc_gpu(i,j,k,4)
!
    rhomean = wmean_gpu(1,1-i,j)
    uumean  = wmean_gpu(2,1-i,j)/rhomean
    vvmean  = wmean_gpu(3,1-i,j)/rhomean
    wwmean  = wmean_gpu(4,1-i,j)/rhomean
    tmean   = p0/rhomean
    rho     = rhomean + rhofluc
    uu      = uumean  + ufluc
    vv      = vvmean  + vfluc
    ww      = wwmean  + wfluc
    rhou    = rho*uu
    rhov    = rho*vv
    rhow    = rho*ww
    w_gpu(1-i,j,k,1) = rho
    w_gpu(1-i,j,k,2) = rhou 
    w_gpu(1-i,j,k,3) = rhov 
    w_gpu(1-i,j,k,4) = rhow 
    w_gpu(1-i,j,k,5) = p0*gm + 0.5_mykind*(rhou**2+rhov**2+rhow**2)/rho
   enddo
  enddo
 enddo
 !@cuf iercuda=cudaDeviceSynchronize()
! 
end subroutine bcrecyc
