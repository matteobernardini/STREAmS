subroutine writestatbl
!
! Writing BL and SBLI statistics
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,l,j99,nkb,npt,ii
 real(mykind), dimension(nx,ny) :: ufav,vfav,wfav
 real(mykind), dimension(ny) :: uvd
 real(mykind) :: aa,alpha,bb,beta,cf,cfinc,delta99,deltav,dely,dstar
 real(mykind) :: dstarinc,dudyw,dx,dy,dyh,dz,fc,ff,ftheta,pdyn,pe,pp,pp2
 real(mykind) :: alfa,prmsw,retau,rethetainc,rethetawall
 real(mykind) :: rho,rho2,rhoe,rhop,rhou,rhov,rhow,ri,rmu,rmue,rmuw,rnuw
 real(mykind) :: shapef,shapefinc,tauw,theta,thetainc,tkel,tt,tt2
 real(mykind) :: tte,ttw,udel,uden,ue,unum,utau,uu,uu2,uup,vv,vv2,ww,ww2
 real(mykind) :: y99,yp,up,uvdp,urmsp,vrmsp,wrmsp,uvp,rhofac,prmsp
 real(mykind), dimension(nx) :: d99_vec,deltav_vec,utau_vec
!
 character(2) :: chstat
!
 if (ncoords(3)==0) then
!
  open(unit=10,file='wavplot_'//chx//'_'//chy//'.dat',form='formatted')
  write(10,*) 'zone i=',nx,', j=',ny
  do j=1,ny
   do i=1,nx
    ufav(i,j) = w_av(13,i,j)/w_av(1,i,j)
    vfav(i,j) = w_av(14,i,j)/w_av(1,i,j)
    wfav(i,j) = w_av(15,i,j)/w_av(1,i,j)
    write(10,100) x(i),y(j),(w_av(l,i,j),l=1,nvmean)
   enddo
  enddo
  close(10)
!
! Mean boundary layer properties
!
  udel = 0.99_mykind*u0
!
  open(10,file='cf_'//chx//'.dat',form='formatted')
  do i=1,nx
   dudyw = (-22._mykind*ufav(i,1)+36._mykind*ufav(i,2)-18._mykind*ufav(i,3)+ 4._mykind*ufav(i,4))/12._mykind
   dy    = (-22._mykind*   y(  1)+36._mykind*   y(  2)-18._mykind*   y(  3)+ 4._mykind*   y(  4))/12._mykind
   dudyw = dudyw/dy
   rhow  = w_av(1,i,1)
   ttw   = w_av(6,i,1)
   rmuw  = w_av(20,i,1)
   tauw  = rmuw*dudyw
   utau  = sqrt(abs(tauw)/rhow)
   rnuw  = rmuw/rhow
   deltav= rnuw/utau
   pdyn  = 0.5_mykind*u0**2
   cf    = tauw/pdyn
   j99   = 1
   do j=1,ny-1
    uu = ufav(i,j)
    if (uu>udel) then
     j99 = j-1
     exit
    endif
   enddo
   dely = y(j99+1)-y(j99)
   unum = udel-ufav(i,j99)
   uden = ufav(i,j99+1)-ufav(i,j99)
   delta99 = y(j99)+dely*(unum/uden) ! b._mykindl._mykind thickness
   retau = delta99/deltav
   prmsp = sqrt(abs(w_av(11,i,1)-w_av(5,i,1)**2))/tauw
   d99_vec(i)    = delta99
   deltav_vec(i) = deltav
   utau_vec(i)   = utau
!
!  Integral boundary layer thicknesses
!
   dstar     = 0.0_mykind
   theta     = 0.0_mykind
   dstarinc  = 0.0_mykind
   thetainc  = 0.0_mykind
   rhoe      = w_av(1,i,j99)
   pe        = w_av(5,i,j99)
   ue        = ufav(i,j99)
   do j=1,j99
    rho  = w_av(1,i,j)/rhoe
    rhop = w_av(1,i,j+1)/rhoe
    uu   = ufav(i,j)/ue
    uup  = ufav(i,j+1)/ue
    dy   = y(j+1)-y(j)
    dyh  = 0.5_mykind*dy
!   Trapezoidal rule
    dstar = dstar       + dyh*((1._mykind-rho*uu)+(1._mykind-rhop*uup))
    theta = theta       + dyh*((rho*uu*(1._mykind-uu))+(rhop*uup*(1._mykind-uup)))
    dstarinc = dstarinc + dyh*((1._mykind-uu)+(1._mykind-uup))
    thetainc = thetainc + dyh*((uu*(1._mykind-uu))+(uup*(1._mykind-uup)))
   enddo
   shapef    = dstar/theta ! Shape factor H
   shapefinc = dstarinc/thetainc ! Incompressible Shape factor H_i
!
!  VAN DRIEST II 
!
   tte  = pe/rhoe
   if (visc_type==1) then
    rmue = sqgmr*tte**vtexp
   else
    rmue = sqgmr*sqrt(tte)*(1._mykind+s2tinf)/(1._mykind+s2tinf/tte)
   endif
!  rfac = pr**(1._mykind/3._mykind)
   ff   = ttw/tte
   aa   = ((rfac*0.2_mykind*rm**2)/ff)**0.5_mykind
   bb   = (1+rfac*0.2_mykind*rm**2-ff)/ff
   alfa = 2*aa**2-bb
   alfa = alfa/(4*aa**2+bb**2)**0.5_mykind
   beta = bb/(4*aa**2+bb**2)**0.5_mykind
   fc   = rfac*0.2_mykind*rm**2
   fc   = fc/(asin(alfa)+asin(beta))**2
   ftheta = rmue/rmuw
   cfinc = cf*fc
   rethetainc = re*thetainc*ftheta
   rethetawall = re*theta*rmue/rmuw
!
   write(10,100) x(i),cf,retau,shapef,shapefinc,delta99,dstar,theta,utau,rethetainc,cfinc,rethetawall,prmsp
  enddo   
  close(10)
!
  do ii = 1,nstatloc
   i = ixstat(ii)
   write(chstat,1002) igxstat(ii)
!
   rhow = w_av(1,i,1)
   delta99 = d99_vec(i)
   deltav  = deltav_vec(i)
   utau    = utau_vec(i)
!
!  Van Driest velocity
!
   uvd(1) = 0._mykind
   do j=2,ny
    uvd(j) = uvd(j-1)+(ufav(i,j)-ufav(i,j-1))*0.5_mykind*(sqrt(w_av(1,i,j)/rhow)+sqrt(w_av(1,i,j-1)/rhow))
   enddo
!
   open(unit=15,file='stat'//chstat//'.prof')
   do j=1,ny
    y99    = y(j)/delta99
    yp     = y(j)/deltav
    up     = ufav(i,j)/utau
    uvdp   = uvd(j)/utau 
    urmsp  = sqrt(abs(w_av(16,i,j)/w_av(1,i,j)-ufav(i,j)*ufav(i,j)))/utau
    vrmsp  = sqrt(abs(w_av(17,i,j)/w_av(1,i,j)-vfav(i,j)*vfav(i,j)))/utau
    wrmsp  = sqrt(abs(w_av(18,i,j)/w_av(1,i,j)-wfav(i,j)*wfav(i,j)))/utau
    uvp    =         (w_av(19,i,j)/w_av(1,i,j)-ufav(i,j)*vfav(i,j)) /utau**2
    rhofac = sqrt(w_av(1,i,j)/rhow)
    prmsp  = sqrt(abs(w_av(11,i,j)-w_av(5,i,j)**2))/tauw
    write(15,100) y99,yp,up,uvdp,urmsp,vrmsp,wrmsp,uvp,rhofac,prmsp
   enddo 
   close(15)
!
   enddo
!
 endif
!
 100 format(200ES20.10)
 1002 format(I2.2)
!
end subroutine writestatbl
