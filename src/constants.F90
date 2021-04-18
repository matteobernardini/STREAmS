subroutine constants
!
! Define useful quantities
!
 use mod_streams
 implicit none
!
 integer :: ii,j,jj,kk,m
 real(mykind), dimension(3,ny) :: ylen ! Integral lengthscale in y
 real(mykind), dimension(3,ny) :: zlen ! Integral lengthscale in z
 real(mykind) :: reout
 real(mykind) :: trattmp,rmb
 real(mykind) :: sqrtt,sdivt,sdivt1,sqgmr2
 real(mykind) :: c1,c2,rmfac,rmcfac,rmufac
 real(mykind) :: alfa_channel,rm_infinity,tb_tw,tinf_tw
!      
 write(chx,1003) ncoords(1)
 write(chy,1003) ncoords(2)
 write(chz,1003) ncoords(3)
 1003 format(I3.3)
!
 if (visc_type==2) s2tinf = 110.4_mykind/tref_dimensional
!
 rho0 = 1._mykind
 u0   = sqrt(gamma)*rm
 v0   = 0._mykind
 w0   = 0._mykind
 s0   = 0._mykind
 p0   = 1._mykind
 t0   = 1._mykind
!
 if (iflow==-1) then ! wind tunnel
  reout = retauinflow 
 elseif (iflow==0) then
!
  rmfac = 0.5_mykind*gm1*rm*rm*rfac
  if (trat<=0._mykind) then
   tb_tw  = 1._mykind+rmfac
  else
   tb_tw  = trat
   target_tbulk = tb_tw
  endif
  rmb     = rm/sqrt(tb_tw) ! Mach bulk based on Tbulk
  if (masterproc) write(*,*) 'Mach bulk (defined with Tbulk) = ', rmb
  trattmp = 1._mykind/tb_tw/(1._mykind+0.5_mykind*gm1*rmb*rmb*rfac) ! Ratio between Twall and Trecovery based on Tbulk
  if (masterproc) write(*,*) 'Ratio between Twall and Trecovery based on Tbulk =', trattmp
  call get_reynolds_chan(ny/2,y(1:ny/2),yn(1:ny/2+1),retauinflow,rmb,trattmp,s2tinf,reout)
  if (masterproc) write(*,*) 'Re bulk (viscosity evaluated at Twall) = ', reout
!
! Connection to spatially developing compressible channel (see Pirozzoli, PRF 2019)
  if (.true.) then
   rm_infinity  = rm/sqrt(1._mykind-rmfac) ! Mach number at infinity
   rmfac = 0.5_mykind*gm1*rm_infinity*rm_infinity*rfac
   alfa_channel = (tb_tw*(1._mykind+rmfac)-1._mykind)/rmfac ! Alfa value (see PRF Pirozzoli)
   if (masterproc) write(*,*) 'Mach infinity = ', rm_infinity
   if (masterproc) write(*,*) 'Alfa channel = ', alfa_channel
   tinf_tw = 1._mykind/(1._mykind+0.5_mykind*gm1*rm_infinity*rm_infinity*rfac)
   if (visc_type==1) then
    rmufac  = tinf_tw**vtexp
   else
    sqrtt   =  sqrt(tinf_tw)
    sdivt   =  s2tinf/tinf_tw
    sdivt1  =  1._mykind+sdivt
    rmufac  = (1._mykind+s2tinf)*sqrtt/sdivt1  ! molecular viscosity
   endif
   if (masterproc) write(*,*) 'Re bulk (viscosity evaluated at t_infinity) = ', reout/rmufac
  endif
!
  if (visc_type==1) then
   rmufac  = tb_tw**vtexp
  else
   sqrtt   = sqrt(tb_tw)
   sdivt   = s2tinf/tb_tw
   sdivt1  = 1._mykind+sdivt
   rmufac  = (1._mykind+s2tinf)*sqrtt/sdivt1  ! molecular viscosity
  endif
!
  if (trat<=0._mykind) then
   if (masterproc) write(*,*) 'Estimated Re bulk (viscosity evaluated at tbulk) = ', reout/rmufac
  else
   if (masterproc) write(*,*) 'Re bulk (viscosity evaluated at tbulk) = ', reout/rmufac
  endif
!
 else
  call get_reynolds(ny,y(1:ny),retauinflow,rm,trat,s2tinf,reout)
  if (masterproc) write(*,*) 'Re = ', reout
 endif
!
 re = reout
 sqgmr = sqrt(gamma)*rm/re
!
! Adiabatic wall temperature and effective wall temperature
!
 taw   = t0*(1._mykind+0.5_mykind*gm1*rm*rm*rfac)
 twall = taw*trat
!
 if (iflow==-1) twall = t0
 if (iflow==0)  twall = t0
!
! winf is the vector of conservative variables at freestream conditions
!
 winf(1) = rho0
 winf(2) = rho0*u0
 winf(3) = rho0*v0
 winf(4) = rho0*w0
 winf(5) = p0*gm+0.5_mykind*(winf(2)**2)/winf(1)
!
 if (iflow==2) call osw() ! Oblique shock wave
!
! Coefficients for computation of convective terms
!
 dcoe(1,1) = 1._mykind/2._mykind
!
 dcoe(1,2) =  2._mykind/3._mykind
 dcoe(2,2) = -1._mykind/12._mykind
!
 dcoe(1,3) =  3._mykind/4._mykind
 dcoe(2,3) = -3._mykind/20._mykind
 dcoe(3,3) =  1._mykind/60._mykind
!
 dcoe(1,4) =  4._mykind/5._mykind
 dcoe(2,4) = -1._mykind/5._mykind
 dcoe(3,4) =  4._mykind/105._mykind
 dcoe(4,4) = -1._mykind/280._mykind
!
! Coefficients for computation of first derivatives (used in visflx and visflx_stag)
!
 select case (ivis/2)
 case (1)
  coeff_deriv1(1) = 0.5_mykind
 case (2)
  coeff_deriv1(1) = 2._mykind/3._mykind
  coeff_deriv1(2) = -1._mykind/12._mykind
 case (3)
  coeff_deriv1(1) = 0.75_mykind
  coeff_deriv1(2) = -0.15_mykind
  coeff_deriv1(3) = 1._mykind/60._mykind
 end select
!
! Coefficients for computation of laplacians (used in visflx)
!
 select case (ivis/2)
 case (1)
  coeff_clap(0) = -2._mykind
  coeff_clap(1) =  1._mykind
 case (2)
  coeff_clap(0) = -2.5_mykind
  coeff_clap(1) = 4._mykind/3._mykind
  coeff_clap(2) = -1._mykind/12._mykind
 case (3)
  coeff_clap(0) = -245._mykind/90._mykind
  coeff_clap(1) = 1.5_mykind
  coeff_clap(2) = -0.15_mykind
  coeff_clap(3) = 1._mykind/90._mykind
 endselect
!
! Coefficients for computation of first derivatives staggered (used in visflx_stag)
! Warning: they are not the coefficients for the computation of a cell-centred
! derivative reported in Lele; they are chosen to recover the laplacian
! coefficients on a uniform mesh
!
 select case (ivis/2)
 case (1)
  coeff_deriv1s(1) = 1._mykind
 case (2)
  coeff_deriv1s(1) =  5._mykind/4._mykind
  coeff_deriv1s(2) = -1._mykind/12._mykind
 case (3)
  coeff_deriv1s(1) = 49._mykind/36._mykind
  coeff_deriv1s(2) = -5._mykind/36._mykind
  coeff_deriv1s(3) =  1._mykind/90._mykind
! coeff_deriv1s(1) = 225._mykind/192._mykind
! coeff_deriv1s(2) = -25._mykind/384._mykind
! coeff_deriv1s(3) =   9._mykind/1920._mykind
 end select
!
! Coefficients for mid point interpolation (used in visflx_stag)
!
 select case (ivis/2)
 case (1)
  coeff_midpi(1) = 0.5_mykind
 case (2)
  coeff_midpi(1) =  9._mykind/16._mykind
  coeff_midpi(2) = -1._mykind/16._mykind
 case (3)
  coeff_midpi(1) =  75._mykind/128._mykind
  coeff_midpi(2) = -25._mykind/256._mykind
  coeff_midpi(3) =   3._mykind/256._mykind
 end select
!
 call compute_coeff_xyz_midpi()
!
!Coefficients for RK 3rd order time integration
!
 gamvec(1) =  8._mykind/15._mykind
 gamvec(2) =  5._mykind/12._mykind
 gamvec(3) =  3._mykind/ 4._mykind
 rhovec(1) =      0._mykind
 rhovec(2) = -17._mykind/60._mykind
 rhovec(3) =  -5._mykind/12._mykind
!
!Initialize random
!
 call init_random_seed(nrank, rand_start)
!
 if (iflow>0) call df_par() ! DF parameters (not needed for channel flow)
!
end subroutine constants    
