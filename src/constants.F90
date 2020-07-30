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
 real(mykind) :: tinf
 real(mykind) :: reout
!      
 write(chx,1003) ncoords(1)
 write(chy,1003) ncoords(2)
 write(chz,1003) ncoords(3)
 1003 format(I3.3)
!
 if (s2tinf==0._mykind) then
  tinf     = 300._mykind/(1._mykind+0.5*gm1*rm*rm)
  s2tinf   = 110.4_mykind/tinf
 endif
!
 if (iflow==-1) then ! wind tunnel
  reout = retauinflow 
 elseif (iflow==0) then
  call get_reynolds_channel(ny/2,y(1:ny/2),yn(1:ny/2+1),retauinflow,rm,trat,s2tinf,reout)
  if (masterproc) write(*,*) 'Re bulk (based on half channel) = ', reout
 else
  call get_reynolds(ny,y(1:ny),retauinflow,rm,trat,s2tinf,reout)
  if (masterproc) write(*,*) 'Re = ', reout
 endif
 re = reout
 sqgmr = sqrt(gamma)*rm/re
!
 rho0 = 1._mykind
 u0   = sqrt(gamma)*rm
 v0   = 0._mykind
 w0   = 0._mykind
 s0   = 0._mykind
 p0   = 1._mykind
 t0   = 1._mykind
!
! Adiabatic wall temperature
!
 taw  = t0*(1._mykind+0.5_mykind*gm1*rm*rm*rfac)
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
!Coefficients for convective terms
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
