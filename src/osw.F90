subroutine osw
!
! Implementation of shock relations
!
 use mod_streams
 implicit none
!
 integer :: ii
 real(mykind) :: sigma
 real(mykind) :: rmn,rrat,prat
 real(mykind) :: u0n,u0t,u1n,u1t,u1mod,u1,v1,w1,rho1,p1
!
 deflec = deflec*pi/180._mykind
 call shock_angle(rm,deflec,sigma)
!
 rmn    = rm*sin(sigma)
 rrat   = (gamma+1._mykind)*rmn**2/(2._mykind+(gamma-1._mykind)*rmn**2)
 prat   = 1._mykind+2._mykind*gamma/(gamma+1._mykind)*(rmn**2-1._mykind)
 u0n    = u0*sin(sigma)
 u0t    = u0*cos(sigma)
 u1n    = u0n/rrat
 u1t    = u0t
 u1mod  = sqrt(u1n**2+u1t**2)
 u1     =  u1mod*cos(deflec)
 v1     = -u1mod*sin(deflec)
 w1     = 0._mykind
 rho1   = rho0*rrat
 p1     = p0*prat
!
 winf1(1) = rho1
 winf1(2) = rho1*u1
 winf1(3) = rho1*v1
 winf1(4) = rho1*w1
 winf1(5) = p1*gm+0.5_mykind*(winf1(2)**2+winf1(3)**2+winf1(4)**2)/winf1(1)
!
 thetas = 0.5_mykind*pi-sigma
!
 xsh = ximp-rly*tan(thetas)
!
 call check_input(3)
!
 if (masterproc) write(*,*) 'xsh =', xsh
!
 call locateval(xg(1:nxmax),nxmax,xsh,ii)
 tanhsfac = dxg(ii)
!
end subroutine osw
