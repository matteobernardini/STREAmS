subroutine writestatzbl
!
! Writing BL ans SBLI statistics (at selected time)
!
 use mod_streams
 implicit none
!
 integer :: i,j,k,l,j99,nkb,npt,ii
 real(mykind), dimension(nx,ny) :: ufav,vfav,wfav
 real(mykind) :: aa,alpha,bb,beta,cf,cfinc,delta99,deltav,dely,dstar
 real(mykind) :: dstarinc,dudyw,dx,dy,dyh,dz,fc,ff,ftheta,pdyn,pe,pp,pp2
 real(mykind) :: alfa,prmsw,retau,rethetainc,rethetawall
 real(mykind) :: rho,rho2,rhoe,rhop,rhou,rhov,rhow,ri,rmu,rmue,rmuw,rnuw
 real(mykind) :: shapef,shapefinc,tauw,theta,thetainc,tkel,tt,tt2
 real(mykind) :: tte,ttw,udel,uden,ue,unum,utau,uu,uu2,uup,vv,vv2,ww,ww2
 real(mykind) :: y99,yp,up,uvdp,urmsp,vrmsp,wrmsp,uvp,rhofac,prmsp
!
 character(4) :: chstore
!
 write(chstore,1004) istore
 1004 format(I4.4)
!
 call compute_av
!
 if (ncoords(3)==0) then
!
  do j=1,ny
   do i=1,nx
    ufav(i,j) = w_avzg(13,i,j)/w_avzg(1,i,j)
   enddo
  enddo
!
! Mean boundary layer properties
!
  udel = 0.99_mykind*u0
!
  open(10,file='cfz_'//chx//'_'//chstore//'.dat',form='formatted')
  do i=1,nx
   dudyw = (-22._mykind*ufav(i,1)+36._mykind*ufav(i,2)-18._mykind*ufav(i,3)+ 4._mykind*ufav(i,4))/12._mykind
   dy    = (-22._mykind*   y(  1)+36._mykind*   y(  2)-18._mykind*   y(  3)+ 4._mykind*   y(  4))/12._mykind
   dudyw = dudyw/dy
   rhow  = w_avzg(1,i,1)
   rmuw  = w_avzg(20,i,1)
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
   delta99 = y(j99)+dely*(unum/uden)
   retau = delta99/deltav
!
!  Integral boundary layer thicknesses
!
   dstar     = 0.0_mykind
   theta     = 0.0_mykind
   dstarinc  = 0.0_mykind
   thetainc  = 0.0_mykind
   rhoe      = w_avzg(1,i,j99)
   pe        = w_avzg(5,i,j99)
   ue        = ufav(i,j99)
   do j=1,j99
    rho  = w_avzg(1,i,j)/rhoe
    rhop = w_avzg(1,i,j+1)/rhoe
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
   shapef    = dstar/theta       ! Shape factor H
   shapefinc = dstarinc/thetainc ! Incompressible Shape factor H_i
   write(10,100) x(i),cf,retau,shapef,shapefinc
  enddo   
  close(10)
!
 endif
!
 100 format(200ES20.10)
!
end subroutine writestatzbl
