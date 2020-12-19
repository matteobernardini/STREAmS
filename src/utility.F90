subroutine locateval(xx,n,x,ii)
!
 use mod_streams, only: mykind
 implicit none
!
 integer, intent(in) :: n
 integer, intent(out) :: ii
 real(mykind), dimension(1:n), intent(in) :: xx
 real(mykind), intent(in) :: x
!
 if (x>xx(n)) then
  ii = n
 else
  ii = 0
  do while (xx(ii+1)<x)
   ii = ii+1
  enddo
 endif
!
end subroutine locateval
!
subroutine pol_int(x,y,n,xs,ys)
!
! Polynomial interpolation using Neville's algorithm  
! Order of accuracy of the interpolation is n-1
!
 use mod_streams, only: mykind
 implicit none
! 
 integer, intent(in) :: n
 real(mykind), dimension(n), intent(in) :: x,y
 real(mykind), intent(in) :: xs
 real(mykind), intent(out) :: ys
!
 integer :: i,m
 real(mykind), dimension(n)  :: v,vold
!
 v = y
!
 do m=2,n ! Tableu columns
  vold = v
  do i=1,n+1-m
   v(i) = (xs-x(m+i-1))*vold(i)+(x(i)-xs)*vold(i+1)
   v(i) = v(i)/(x(i)-x(i+m-1))
  enddo
 enddo
 ys = v(1)
!
end subroutine pol_int
!
subroutine init_random_seed(x, rand_start)
!
 implicit none
 integer, intent(in) :: x, rand_start
 integer :: i, n, clock
 integer, dimension(:), allocatable :: seed
 integer, dimension(8) :: values
! 
 call random_seed(size = n)
 allocate(seed(n))
 if (rand_start < 0) then
  !call system_clock(count=clock) ! not reproducible random
  call date_and_time(VALUES=values)
  clock = (values(6)+1)*(values(7)+1)*(values(8)+1)
  seed = clock + (x+1) * (/ (i - 1, i = 1, n) /)
 else
  seed = rand_start + (x+1) * (/ (i - 1, i = 1, n) /) ! reproducible 
 endif
 call random_seed(PUT = seed)
! 
 deallocate(seed)
!
end subroutine init_random_seed
!
!=================================================
!
subroutine get_reynolds(ny,eta,retau,rm,trat,s2tinf,reout)
!
 use mod_streams, only: mykind, tol_iter, visc_type, vtexp
 implicit none
!
 integer, intent(in) :: ny
 real(mykind), intent(in)  :: retau,rm,trat,s2tinf
 real(mykind), intent(out) :: reout
 real(mykind), dimension(ny), intent(in)  :: eta
 
 real(mykind), dimension(ny) :: uplus       ! Incompressible velocity profile
 real(mykind), dimension(ny) :: ucplus      ! Compressible velocity profile
 real(mykind), dimension(ny) :: ucplusold
 real(mykind), dimension(ny) :: density     ! Density profile
 real(mykind), dimension(ny) :: viscosity   ! Viscosity profile
 real(mykind), dimension(ny) :: temperature ! Temperature profile

 real(mykind) :: cf,dstar,gm1,gm1h,res,rfac,rhow,th,tr,tw,uci,ue,uu,yplus,du
 real(mykind) :: pr,alf,s,fuu
 real(mykind) :: gamma
 integer :: j
 real(mykind), external :: velmusker

 uplus = 0._mykind
 do j=2,ny
  yplus    = eta(j)*retau
  uplus(j) = velmusker(yplus,retau) ! Incompressible velocity profile (Musker, AIAA J 1979)
 enddo
!
 gamma = 1.4_mykind
 gm1   = gamma-1._mykind
 gm1h  = 0.5_mykind*gm1
 rfac  = 0.89_mykind                ! Recovery factor
 tr    = 1._mykind+gm1h*rfac*rm**2  ! Recovery temperature
 tw    = trat*tr                    ! Wall temperature
 rhow  = 1._mykind/tw
!alf   = 0.8259_mykind
 pr    = 0.72_mykind                ! Prandtl number
 s     = 1.1_mykind                 ! Reynolds analogy factor
 alf   = s*pr                       ! see Zhang
!
!Iterative loop to find the compressible profile
 ucplus = uplus
 do 
!
  ucplusold = ucplus
  do j=1,ny
   uu = ucplus(j)/ucplus(ny)
!  temperature(j) = tw+(tr-tw)*uu+(1._mykind-tr)*uu**2 ! Walz pag 399 Zhang JFM 2014
   fuu = alf*uu+(1-alf)*uu**2
   temperature(j) = tw+(tr-tw)*fuu+(1._mykind-tr)*uu**2 ! Duan & Martin (see Zhang)
   density(j) = 1._mykind/temperature(j)
   if (visc_type==1) then
    viscosity(j) = temperature(j)**vtexp ! Power-law
   else
    viscosity(j) = sqrt(temperature(j))*(1._mykind+s2tinf)/(1._mykind+s2tinf/temperature(j))
   endif
  enddo
  do j=2,ny
   du        = uplus(j)-uplus(j-1)
   uci       = 0.5_mykind*(sqrt(rhow/density(j))+sqrt(rhow/density(j-1)))
   ucplus(j) = ucplus(j-1)+uci*du
  enddo
  res = 0._mykind
  do j=1,ny
   res = res+abs(ucplus(j)-ucplusold(j))
  enddo
  res = res/ny
  if (res < tol_iter) exit
!  
 enddo ! End of iterative loop
!
 dstar = 0._mykind
 th    = 0._mykind
 ue    = ucplus(ny)
 cf    = 2._mykind*rhow/ue**2
!
 reout = retau*ue/rhow*viscosity(1) 
!
end subroutine get_reynolds
!
!=================================================
!
subroutine meanvelocity(ny,eta,u,rho,t,retau,rm,trat,th,cf)
!
 use mod_streams, only: mykind, tol_iter
 implicit none
!
 integer :: ny
 real(mykind), intent(in)  :: retau,rm,trat
 real(mykind), intent(out) :: th,cf
 real(mykind), dimension(ny), intent(in)  :: eta
 real(mykind), dimension(ny), intent(out) :: u,rho,t
!
 real(mykind), dimension(ny) :: uplus       ! Incompressible velocity profile
 real(mykind), dimension(ny) :: ucplus      ! Compressible velocity profile
 real(mykind), dimension(ny) :: ucplusold
 real(mykind), dimension(ny) :: density     ! Density profile
 real(mykind), dimension(ny) :: temperature ! Temperature profile
!
 real(mykind) :: dstar,dstari,du,dy,gamma,gm1,gm1h,res,rfac,rhow,thi,tr,tw,uci,ue,uu,uum,yplus
 real(mykind) :: pr,alf,s,fuu
 integer :: j
 real(mykind), external :: velmusker
!
 uplus = 0._mykind
 do j=2,ny
  yplus    = eta(j)*retau
  uplus(j) = velmusker(yplus,retau) ! Incompressible velocity profile(Musker, AIAA J 1979)
 enddo
!
 gamma = 1.4_mykind
 gm1   = gamma-1
 gm1h  = 0.5_mykind*gm1
 rfac  = 0.89_mykind                ! Recovery factor
 tr    = 1._mykind+gm1h*rfac*rm**2  ! Recovery temperature
 tw    = trat*tr                    ! Wall temperature
 rhow  = 1._mykind/tw
!alf   = 0.8259_mykind              ! Duan & Martin
 pr    = 0.72_mykind                ! Prandtl number
 s     = 1.1_mykind                 ! Reynolds analogy factor
 alf   = s*pr                       ! see Zhang
!
!Iterative loop to find the compressible profile
 ucplus = uplus
 do 
!
  ucplusold = ucplus
  do j=1,ny
   uu = ucplus(j)/ucplus(ny)
!  temperature(j) = tw+(tr-tw)*uu+(1._mykind-tr)*uu**2 ! Walz pag 399 Zhang JFM 2014
   fuu = alf*uu+(1-alf)*uu**2
   temperature(j) = tw+(tr-tw)*fuu+(1._mykind-tr)*uu**2 ! Duan & Martin (see Zhang)
   density(j) = 1._mykind/temperature(j)
  enddo
  do j=2,ny
   du        = uplus(j)-uplus(j-1)
   uci       = 0.5_mykind*(sqrt(rhow/density(j))+sqrt(rhow/density(j-1)))
   ucplus(j) = ucplus(j-1)+uci*du
  enddo
  res = 0._mykind
  do j=1,ny
   res = res+abs(ucplus(j)-ucplusold(j))
  enddo
  res = res/ny
  if (res < tol_iter) exit
!  
 enddo ! End of iterative loop
!
 dstar = 0._mykind
 th    = 0._mykind
 ue    = ucplus(ny)
 cf    = 2._mykind*rhow/ue**2
!
 do j=2,ny
  uu     = ucplus(j)  /ue
  uum    = ucplus(j-1)/ue
  dy     = eta(j)-eta(j-1)
  dstari = 0.5_mykind*((1._mykind-density(j)*uu)+(1._mykind-density(j-1)*uum))
  dstar  = dstar+dstari*dy
  thi    = 0.5_mykind*(density(j)*uu*(1._mykind-uu)+density(j-1)*uum*(1._mykind-uum))
  th     = th+thi*dy
 enddo
!
! du= (-22._mykind*ucplus(1)+36._mykind*ucplus(2)-18._mykind*ucplus(3)+ 4._mykind*ucplus(4))/12._mykind
! dy= (-22._mykind*eta(1)+36._mykind*eta(2)-18._mykind*eta(3)+ 4._mykind*eta(4))/12._mykind
! tauw = du/dy/ue/(retau*ue/rhow)
! print *, 2._mykind*tauw, cf
!
! m    = 2
! ue95 = 0.95_mykind*ue
! call locate(ucplus,ny,ue95,j95)
! jj95 = min(max(j95-(m-1)/2,1),ny+1-m)
! call polint(ucplus(jj95),eta(jj95),m,ue95,d95,d95err)
! ue99 = 0.99_mykind*ue
! call locate(ucplus,ny,ue99,j99)
! jj99 = min(max(j99-(m-1)/2,1),ny+1-m)
! call polint(ucplus(jj99),eta(jj99),m,ue99,d99,d99err)
!
! write(*,*) 'Effective friction Reynolds number (d99) = ', retau*d99
! write(*,*) 'Re_delta2 = ', retau*ue*th/rhow
!
 do j=1,ny
  u(j)   = ucplus(j)/ucplus(ny)
  rho(j) = density(j)
  t(j)   = temperature(j)
 enddo
!
end subroutine meanvelocity
!
!=================================================
!
 function velmusker(yplus,retau)
!
 use mod_streams, only: mykind
 implicit none
!
 real(mykind), intent(in) :: yplus,retau
 real(mykind) :: yp,eta,pi_wake
 real(mykind) :: velmusker
!
 yp        = yplus
 eta       = yp/retau
 eta       = min(1._mykind,eta)
 yp        = eta*retau
 pi_wake   = 0.434_mykind
 velmusker = 5.424_mykind*atan((2*yp-8.15_mykind)/16.7_mykind)+ &
  &log10((yp+10.6_mykind)**9.6_mykind/(yp**2-8.15_mykind*yp+86)**2)-3.51132976630723_mykind &
  & +2.44_mykind*(pi_wake*(6*eta**2-4*eta**3)+(eta**2*(1-eta)))
!
end function velmusker
!
subroutine gasdev_s(harvest)
!
 use mod_streams, only: mykind
 implicit none
!
 real(mykind), intent(out) :: harvest

! Local variables
 real(mykind)          :: rsq, v1, v2
 real(mykind), save    :: g
 logical, save :: gaus_stored = .false.
!
 if (gaus_stored) then
  harvest = g
  gaus_stored = .false.
 else
  do
   call random_number(v1)
   call random_number(v2)
   v1 = 2.0_mykind*v1 - 1.0_mykind
   v2 = 2.0_mykind*v2 - 1.0_mykind
   rsq = v1**2 + v2**2
   if (rsq > 0.0_mykind .and. rsq < 1.0_mykind) exit
  end do
  rsq = sqrt(-2.0_mykind*log(rsq)/rsq)
  harvest = v1*rsq
  g = v2*rsq
  gaus_stored = .true.
 end if
!
end subroutine gasdev_s
!
subroutine shock_angle(rm,delta,sigma)
!
 use mod_streams, only: mykind,pi,gamma,tol_iter
 implicit none
!
 real(mykind) :: rm,delta,sigma,diffs,rden,rm2,rnum,sigmao,sins
!
!For weak oblique shock
!
 sigma = 1._mykind*pi/180._mykind
!
 do
  sigmao = sigma
  rm2   = rm*rm
  rnum  = 1._mykind+0.5_mykind*tan(delta)*tan(sigma)*(rm2*(gamma+cos(2._mykind*sigma))+2._mykind)
  rden  = rm2 
  sins  = sqrt(rnum/rden)
  sigma = asin(sins)
  diffs = abs(sigma-sigmao)
  if (diffs < tol_iter) exit
 enddo
!
!For strong oblique shock
!
!sigma = 90._mykind*pi/180._mykind
!
!do
! sigmao = sigma
! sins  = sin(sigma)
! rm2   = rm*rm
! rnum  = rm2*sins*sins-1._mykind
! rden  = rm2*(gamma+cos(2._mykind*sigma))+2._mykind
! tans  = 2._mykind/tan(delta)*rnum/rden
! sigma = atan(tans)
! diffs = abs(sigma-sigmao)
! if (diffs<tol_iter) exit
!enddo
!
end subroutine shock_angle
!
subroutine get_reynolds_chan(ny,eta,yn,retau,rm,trat,s2tinf,reout)
!
 use mod_streams, only: mykind, tol_iter, visc_type, vtexp
 implicit none
!
 integer, intent(in) :: ny
 real(mykind), intent(in)  :: retau,rm,trat,s2tinf
 real(mykind), intent(out) :: reout
 real(mykind), dimension(ny), intent(in)  :: eta
 real(mykind), dimension(ny+1), intent(in) :: yn
! 
 real(mykind), dimension(ny) :: uplus       ! Incompressible velocity profile
 real(mykind), dimension(ny) :: ucplus      ! Compressible velocity profile
 real(mykind), dimension(ny) :: ucplusold
 real(mykind), dimension(ny) :: density     ! Density profile
 real(mykind), dimension(ny) :: viscosity   ! Viscosity profile
 real(mykind), dimension(ny) :: temperature ! Temperature profile

 real(mykind) :: gm1,gm1h,res,rfac,rhow,tr,tw,uci,ue,uu,yplus,du,dy,te
 real(mykind) :: pr,alf,s,fuu,rhosum,rhousum,ubulk
 real(mykind) :: gamma
 integer :: j,l
 real(mykind), external :: velmusker
!
 uplus = 0._mykind
 do j=2,ny
  yplus    = (1._mykind+eta(j))*retau
! uplus(j) = velmusker(yplus,retau) ! Incompressible velocity profile (Musker, AIAA J 1979)
  uplus(j) = 5.5_mykind+log(yplus)/0.4_mykind
 enddo
!
 gamma = 1.4_mykind
 gm1   = gamma-1._mykind
 gm1h  = 0.5_mykind*gm1
 rfac  = 0.89_mykind                ! Recovery factor
 te    = 1._mykind
 tr    = 1._mykind+gm1h*rfac*rm**2 ! Recovery temperature
 tw    = trat*tr                    ! Wall temperature
 rhow  = 1._mykind/tw
!alf   = 0.8259_mykind
 pr    = 0.72_mykind                ! Prandtl number
 s     = 1.1_mykind                 ! Reynolds analogy factor
 alf   = s*pr                       ! see Zhang
!
! Iterative loop to find the compressible profile
 ucplus = uplus
 do 
!
  ucplusold = ucplus
  do j=1,ny
   uu = ucplus(j)/ucplus(ny)
!  temperature(j) = tw+(tr-tw)*uu+(1._mykind-tr)*uu**2 ! Walz pag 399 Zhang JFM 2014
   fuu = alf*uu+(1-alf)*uu**2
   temperature(j) = tw+(tr-tw)*fuu+(te-tr)*uu**2 ! Duan & Martin (see Zhang)
   density(j) = 1._mykind/temperature(j)
   if (visc_type==1) then
    viscosity(j) = temperature(j)**vtexp ! Power-law
   else
    viscosity(j) = sqrt(temperature(j))*(1._mykind+s2tinf)/(1._mykind+s2tinf/temperature(j))
   endif
  enddo
  do j=2,ny
   du        = uplus(j)-uplus(j-1)
   uci       = 0.5_mykind*(sqrt(rhow/density(j))+sqrt(rhow/density(j-1)))
   ucplus(j) = ucplus(j-1)+uci*du
  enddo
  res = 0._mykind
  do j=1,ny
   res = res+abs(ucplus(j)-ucplusold(j))
  enddo
  res = res/ny
  if (res < tol_iter) exit
!   
 enddo ! End of iterative loop
!
 rhosum   = 0.
 rhousum  = 0.
 do j=1,ny
  dy      = yn(j+1)-yn(j)
  rhosum  = rhosum+density(j)*dy
  rhousum = rhousum+density(j)*ucplus(j)*dy
 enddo
 ubulk = rhousum/rhosum
!
 ue    = ucplus(ny)
!reout = retau*ue/rhow*viscosity(1) 
!reout = retau*ubulk*rhosum/rhow*viscosity(1) 
 reout = retau*ubulk*rhosum/rhow
!
end subroutine get_reynolds_chan
