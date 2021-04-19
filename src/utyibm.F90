subroutine bilinear(x0,x1,y0,y1,V,xx,yy,c)
!
 use mod_streams, only: mykind
 implicit none
!
 real(mykind), intent(in) :: x0,x1,y0,y1,xx,yy
 real(mykind), intent(out) :: c
 real(mykind), dimension(2,2), intent(in) :: V
 real(mykind) :: xd,yd,c1,c2
!
 xd = (xx-x0)/(x1-x0)
 yd = (yy-y0)/(y1-y0)
!
 c1 = V(1,1)*(1-xd)+V(2,1)*xd
 c2 = V(1,2)*(1-xd)+V(2,2)*xd
!
 c  = c1*(1-yd)+c2*yd
!
end
!
subroutine trilinear(x0,x1,y0,y1,z0,z1,V,xx,yy,zz,c)
!
 use mod_streams, only: mykind
 implicit none
!
 real(mykind), intent(in) :: x0,x1,y0,y1,z0,z1,xx,yy,zz
 real(mykind), intent(out) :: c
 real(mykind), dimension(2,2,2), intent(in) :: V
 real(mykind) :: xd,yd,zd,c1,c2,c11,c12,c21,c22
!
 xd = (xx-x0)/(x1-x0)
 yd = (yy-y0)/(y1-y0)
 zd = (zz-z0)/(z1-z0)
!
 c11 = V(1,1,1)*(1-xd)+V(2,1,1)*xd
 c12 = V(1,1,2)*(1-xd)+V(2,1,2)*xd
 c21 = V(1,2,1)*(1-xd)+V(2,2,1)*xd
 c22 = V(1,2,2)*(1-xd)+V(2,2,2)*xd
!
 c1  = c11*(1-yd)+c21*yd
 c2  = c21*(1-yd)+c22*yd
!
 c   = c1*(1-zd)+c2*zd
!
end

subroutine invert4x4(m,ierr)
use mod_streams, only: mykind
 implicit none
 real(mykind) , intent(inout) :: m(4,4)
 integer, intent(out) :: ierr
 integer :: i,j
 real(mykind) :: det, inv(4,4)

    inv(1,1) = m(2,2) * m(3,3) * m(4,4) - &
               m(2,2) * m(3,4) * m(4,3) - &
               m(3,2) * m(2,3) * m(4,4) + &
               m(3,2) * m(2,4) * m(4,3) + &
               m(4,2) * m(2,3) * m(3,4) - &
               m(4,2) * m(2,4) * m(3,3)

    inv(2,1) =-m(2,1) * m(3,3) * m(4,4) + &
               m(2,1) * m(3,4) * m(4,3) + &
               m(3,1) * m(2,3) * m(4,4) - &
               m(3,1) * m(2,4) * m(4,3) - &
               m(4,1) * m(2,3) * m(3,4) + &
               m(4,1) * m(2,4) * m(3,3)

    inv(3,1) = m(2,1) * m(3,2) * m(4,4) - &
               m(2,1) * m(3,4) * m(4,2) - &
               m(3,1) * m(2,2) * m(4,4) + &
               m(3,1) * m(2,4) * m(4,2) + &
               m(4,1) * m(2,2) * m(3,4) - &
               m(4,1) * m(2,4) * m(3,2)

    inv(4,1) =-m(2,1) * m(3,2) * m(4,3) + &
               m(2,1) * m(3,3) * m(4,2) + &
               m(3,1) * m(2,2) * m(4,3) - &
               m(3,1) * m(2,3) * m(4,2) - &
               m(4,1) * m(2,2) * m(3,3) + &
               m(4,1) * m(2,3) * m(3,2)

    inv(1,2) =-m(1,2) * m(3,3) * m(4,4) + &
               m(1,2) * m(3,4) * m(4,3) + &
               m(3,2) * m(1,3) * m(4,4) - &
               m(3,2) * m(1,4) * m(4,3) - &
               m(4,2) * m(1,3) * m(3,4) + &
               m(4,2) * m(1,4) * m(3,3)

    inv(2,2) = m(1,1) * m(3,3) * m(4,4) - &
               m(1,1) * m(3,4) * m(4,3) - &
               m(3,1) * m(1,3) * m(4,4) + &
               m(3,1) * m(1,4) * m(4,3) + &
               m(4,1) * m(1,3) * m(3,4) - &
               m(4,1) * m(1,4) * m(3,3)

    inv(3,2) =-m(1,1) * m(3,2) * m(4,4) + &
               m(1,1) * m(3,4) * m(4,2) + &
               m(3,1) * m(1,2) * m(4,4) - &
               m(3,1) * m(1,4) * m(4,2) - &
               m(4,1) * m(1,2) * m(3,4) + &
               m(4,1) * m(1,4) * m(3,2)

    inv(4,2) = m(1,1) * m(3,2) * m(4,3) - &
               m(1,1) * m(3,3) * m(4,2) - &
               m(3,1) * m(1,2) * m(4,3) + &
               m(3,1) * m(1,3) * m(4,2) + &
               m(4,1) * m(1,2) * m(3,3) - &
               m(4,1) * m(1,3) * m(3,2)

    inv(1,3) = m(1,2) * m(2,3) * m(4,4) - &
               m(1,2) * m(2,4) * m(4,3) - &
               m(2,2) * m(1,3) * m(4,4) + &
               m(2,2) * m(1,4) * m(4,3) + &
               m(4,2) * m(1,3) * m(2,4) - &
               m(4,2) * m(1,4) * m(2,3)

    inv(2,3) =-m(1,1) * m(2,3) * m(4,4) + &
               m(1,1) * m(2,4) * m(4,3) + &
               m(2,1) * m(1,3) * m(4,4) - &
               m(2,1) * m(1,4) * m(4,3) - &
               m(4,1) * m(1,3) * m(2,4) + &
               m(4,1) * m(1,4) * m(2,3)

    inv(3,3) = m(1,1) * m(2,2) * m(4,4) - &
               m(1,1) * m(2,4) * m(4,2) - &
               m(2,1) * m(1,2) * m(4,4) + &
               m(2,1) * m(1,4) * m(4,2) + &
               m(4,1) * m(1,2) * m(2,4) - &
               m(4,1) * m(1,4) * m(2,2)

    inv(4,3) =-m(1,1) * m(2,2) * m(4,3) + &
               m(1,1) * m(2,3) * m(4,2) + &
               m(2,1) * m(1,2) * m(4,3) - &
               m(2,1) * m(1,3) * m(4,2) - &
               m(4,1) * m(1,2) * m(2,3) + &
               m(4,1) * m(1,3) * m(2,2)

    inv(1,4) =-m(1,2) * m(2,3) * m(3,4) + &
               m(1,2) * m(2,4) * m(3,3) + &
               m(2,2) * m(1,3) * m(3,4) - &
               m(2,2) * m(1,4) * m(3,3) - &
               m(3,2) * m(1,3) * m(2,4) + &
               m(3,2) * m(1,4) * m(2,3)

    inv(2,4) = m(1,1) * m(2,3) * m(3,4) - &
               m(1,1) * m(2,4) * m(3,3) - &
               m(2,1) * m(1,3) * m(3,4) + &
               m(2,1) * m(1,4) * m(3,3) + &
               m(3,1) * m(1,3) * m(2,4) - &
               m(3,1) * m(1,4) * m(2,3)

    inv(3,4) =-m(1,1) * m(2,2) * m(3,4) + &
               m(1,1) * m(2,4) * m(3,2) + &
               m(2,1) * m(1,2) * m(3,4) - &
               m(2,1) * m(1,4) * m(3,2) - &
               m(3,1) * m(1,2) * m(2,4) + &
               m(3,1) * m(1,4) * m(2,2)

    inv(4,4) = m(1,1) * m(2,2) * m(3,3) - &
               m(1,1) * m(2,3) * m(3,2) - &
               m(2,1) * m(1,2) * m(3,3) + &
               m(2,1) * m(1,3) * m(3,2) + &
               m(3,1) * m(1,2) * m(2,3) - &
               m(3,1) * m(1,3) * m(2,2)

    det = m(1,1) * inv(1,1) + m(1,2) * inv(2,1) + &
          m(1,3) * inv(3,1) + m(1,4) * inv(4,1)

    ierr = 0
    if (det.eq.0._mykind) then
      ierr = 1
      write(*,*) 'Error in invert4x4'
      return
    end if

    det = 1._mykind/det

    do i = 1,4
      do j = 1,4
        m(i,j) = inv(i,j) * det
      enddo
    enddo

  return
  end subroutine invert4x4
!
      SUBROUTINE LUDCMP(A,N,NP,INDX,D)
      use mod_streams, only: mykind
      implicit none
      integer, parameter :: NMAX = 100
      real(mykind), parameter :: TINY = 1.0D-20
      integer, intent(in) :: N,NP
      integer, dimension(N), intent(out) :: INDX
      real(mykind), intent(out) :: D
      real(mykind), dimension(NP,NP), intent(inout) :: A
      real(mykind), dimension(NMAX) :: VV 
      integer :: I,J,K,IMAX
      real(mykind) :: SUM,DUM,AAMAX
!
      D=1._mykind
      DO 12 I=1,N
        AAMAX=0._mykind
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0._mykind) PAUSE 'Singular matrix.'
        VV(I)=1./AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0._mykind)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0._mykind)A(N,N)=TINY
      RETURN
      END
!
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      use mod_streams, only: mykind
      implicit none
      integer, intent(in) :: N,NP
      integer, dimension(N), intent(in) :: INDX
      real(mykind), dimension(NP,NP), intent(in) :: A
      real(mykind), dimension(N), intent(inout) :: B
      integer :: II,I,LL,J
      real(mykind) :: SUM
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
!
subroutine invmat(mat,n)
!
 use mod_streams, only: mykind
 implicit none
!
 integer, intent(in) :: n
 real(mykind), dimension(n,n), intent(inout) :: mat
 integer :: i,j
 integer, dimension(n) :: indx
 real(mykind), dimension(n,n) :: y
 real(mykind) :: d
!
 y = 0._mykind
 do i=1,n
  y(i,i) = 1._mykind
 enddo
 call ludcmp(mat,n,n,indx,d)
 do j=1,n
  call lubksb(mat,n,n,indx,y(1,j))
 enddo
 mat = y
!
end subroutine invmat
!
subroutine detmat(mat,n,det)
!
 use mod_streams, only: mykind
 implicit none
!
 integer, intent(in) :: n
 real(mykind), dimension(n,n), intent(inout) :: mat
 real(mykind), intent(out) :: det
 integer :: j
 integer, dimension(n) :: indx
 real(mykind) :: d
!
 call ludcmp(mat,n,n,indx,d)
 do j=1,n
  d = d*mat(j,j)
 enddo
 det = d
!
end subroutine detmat

FUNCTION GETDET(A,N)
! A function written in FORTRAN77 to calculate determinant of a square matrix

! Passed parameters:
! A = the matrix
! N = dimension of the square matrix
!      
! A modification of a code originally written by Ashwith J. Rego, available from
! http://www.dreamincode.net/code/snippet1273.htm
! Modified by Syeilendra Pramuditya, available from http://wp.me/p61TQ-zb
! Last modified on January 13, 2011 
!
 use mod_streams, only: mykind
 implicit none
!
 integer, intent(in) :: N
 real(mykind), dimension(n,n), intent(in) :: A
 real(mykind) :: getdet
 real(mykind), dimension(n,n) :: elem
 real(mykind) :: m,temp
 integer :: i,j,k,l
 logical :: DETEXISTS

      DO I=1,N
      DO J=1,N
      ELEM(I,J)=A(I,J)
      END DO
      END DO
        
      DETEXISTS = .TRUE.
      L = 1
      !CONVERT TO UPPER TRIANGULAR FORM
      DO K = 1, N-1
	  IF (ABS(ELEM(K,K)).LE.1.0D-20) THEN
	   DETEXISTS = .FALSE.
	   DO I = K+1, N
	    IF (ELEM(I,K).NE.0.0) THEN
      
	     DO J = 1, N
	      TEMP = ELEM(I,J)
	      ELEM(I,J)= ELEM(K,J)
	      ELEM(K,J) = TEMP
	     END DO
      
	     DETEXISTS = .TRUE.
	     L=-L
	     EXIT
      
	    END IF
      
	   END DO
	   IF (DETEXISTS .EQV. .FALSE.) THEN
	    GETDET = 0
	    RETURN
	   END IF
          END IF
          DO J = K+1, N
	   M = ELEM(J,K)/ELEM(K,K)
	   DO I = K+1, N
	    ELEM(J,I) = ELEM(J,I) - M*ELEM(K,I)
	   END DO
	  END DO
      END DO
	
      !CALCULATE DETERMINANT BY FINDING PRODUCT OF DIAGONAL ELEMENTS
      GETDET = L
      DO I = 1, N
	  GETDET = GETDET * ELEM(I,I)
      END DO
	
      END        
